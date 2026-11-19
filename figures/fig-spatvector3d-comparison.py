"""Builds the three-panel comparison figure of the ocean3d data model:

  A  SpatVector3D  a polygon with a curved top and bottom surface
  B  SpatEnvelope  A evaluated on a grid (vector3d_to_envelope()),
                   bounded by the seafloor
  C  SpatVoxel     B sliced at standard depths (envelope_to_voxel())

All three panels are drawn from one synthetic example, so B is derived from A
and C from B. Run from the project root:

  python3 figures/fig-spatvector3d-comparison.py \
    vignettes/figures/Comparison_SpatVector3D_SpatEnvelope_SpatVoxel.png
  python3 figures/fig-spatvector3d-comparison.py --no-arrows \
    vignettes/figures/Comparison_SpatVector3D_SpatEnvelope_SpatVoxel_no_arrows.png

--no-arrows leaves out the conversion-function arrows between the panels.
Needs numpy, matplotlib and the Noto Sans fonts.
"""
import argparse

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import (LightSource, LinearSegmentedColormap,  # noqa: E402
                               ListedColormap, Normalize)
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import (Ellipse, FancyArrowPatch,  # noqa: E402
                                FancyBboxPatch, Rectangle)
from matplotlib.path import Path  # noqa: E402
from mpl_toolkits.mplot3d import proj3d  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("out")
parser.add_argument("--no-arrows", action="store_true")
parser.add_argument("--dpi", type=int, default=200)
args = parser.parse_args()

FONT_DIR = "/usr/share/fonts/truetype/noto/"
for f in ("Regular", "Bold", "SemiBold"):
    font_manager.fontManager.addfont(f"{FONT_DIR}NotoSans-{f}.ttf")
plt.rcParams["font.family"] = "Noto Sans"

# ---- colours ----------------------------------------------------------------
reds = LinearSegmentedColormap.from_list(
    "depth", ["#fff4f3", "#f7b9b6", "#ec7a74", "#d23a35", "#9c0b10", "#6e0007"])
depth_norm = Normalize(0, 600)
greys = LinearSegmentedColormap.from_list(
    "seafloor", ["#ffffff", "#e3e5e7", "#b9bec3", "#8d959e", "#5c6875"])
seafloor_norm = Normalize(0, 800)
PRESENT, ABSENT = "#ef6f69", "#ececec"
PANELS = {  # badge fill, border
    "A": ("#e5dcf7", "#cdbeeb"),
    "B": ("#fed7d5", "#f6c9c6"),
    "C": ("#dbe7f8", "#bce2fa"),
}

# ---- synthetic example --------------------------------------------------------
# map units: x 0-12 (east), y 0-10 (north); depths in metres, positive down
ZMAX = 800
ELEV, AZIM = 12, -72  # low, mostly from the south, so depth layers show
cx, cy, ax_x, ax_y = 6.0, 5.0, 4.9, 3.9


def radius(th):
    return (1 + 0.10 * np.sin(3 * th + 0.4) + 0.07 * np.cos(2 * th - 0.8)
            + 0.04 * np.sin(5 * th))


def to_xy(r, th):
    R = radius(th) * r
    return cx + ax_x * R * np.cos(th), cy + ax_y * R * np.sin(th)


def d2(x, y):
    return ((x - cx) / ax_x) ** 2 + ((y - cy) / ax_y) ** 2


def top_depth(x, y):
    """Dome: shallowest (~25 m) in the middle, ~160 m at the edge."""
    return 25 + 135 * np.clip(d2(x, y), 0, 1.4) / 1.2


def bottom_depth(x, y):
    """Deeper to the north, with a bulge in the middle."""
    return 250 + 24 * (y - 1.0) + 180 * np.exp(-2.0 * d2(x, y))


def seafloor_depth(x, y):
    sf = (640 - 110 * np.exp(-((x - 2) ** 2 + (y - 7.5) ** 2) / 5)
          + 120 * np.exp(-((x - 9.5) ** 2 + (y - 1.5) ** 2) / 9)
          + 25 * np.sin(1.3 * x) * np.cos(0.9 * y))
    return np.clip(sf, 570, 790)


th = np.linspace(0, 2 * np.pi, 181)
bx, by = to_xy(np.ones_like(th), th)  # polygon boundary
verts_idx = np.arange(0, 180, 12)       # vertices shown on the footprint map

# template grid for B and C
NX, NY = 16, 12
xe, ye = np.linspace(0, 12, NX + 1), np.linspace(0, 10, NY + 1)
xc, yc = np.meshgrid((xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2)  # [ny, nx]
inside = Path(np.c_[bx, by]).contains_points(
    np.c_[xc.ravel(), yc.ravel()]).reshape(xc.shape)
seafloor = np.round(seafloor_depth(xc, yc), -1)
# vector3d_to_envelope(): top and bottom at each cell centre, capped by seafloor
depth_min = np.where(inside, np.round(top_depth(xc, yc)), np.nan)
depth_max = np.where(inside, np.minimum(np.round(bottom_depth(xc, yc)), seafloor),
                     np.nan)
# envelope_to_voxel() with its default bounds = "top": each level stands for
# the slab from it down to the next level, and a cell occupies the level when
# its envelope shares some water with that slab. The deepest level has no
# next level, so it counts only when the envelope reaches below it.
STD_DEPTHS = [50, 150, 250, 350, 500]
slab_lower = STD_DEPTHS
slab_upper = STD_DEPTHS[1:] + STD_DEPTHS[-1:]
presence = {}
for i, z in enumerate(STD_DEPTHS):
    if i == len(STD_DEPTHS) - 1:
        overlap = (depth_min <= slab_upper[i]) & (depth_max > slab_lower[i])
    else:
        overlap = (depth_min < slab_upper[i]) & (depth_max > slab_lower[i])
    presence[z] = inside & overlap

# ---- figure layout ---------------------------------------------------------------
PW, H = 7.9, 8.87  # panel slot width and figure height, inches
X0 = [0.0, 8.0, 16.0] if args.no_arrows else [0.0, 10.0, 19.74]
W = X0[-1] + PW
fig = plt.figure(figsize=(W, H))


def pfig(i, x, y):
    """Panel-relative (0-1 of the panel slot) to figure-fraction coordinates."""
    return ((X0[i] + x * PW) / W, y)


def padd(i, rect, **kw):
    x, y = pfig(i, rect[0], rect[1])
    return fig.add_axes([x, y, rect[2] * PW / W, rect[3]], **kw)


def dx(i, frac):
    """A panel-width fraction as a figure-width fraction."""
    return frac * PW / W


def panel_frame(i, letter, title, badge_w):
    fill, edge = PANELS[letter]
    fig.add_artist(FancyBboxPatch(
        ((X0[i] + 0.12) / W, 0.15 / H), (PW - 0.24) / W, (H - 0.19) / H,
        boxstyle="round,pad=0,rounding_size=0.008", transform=fig.transFigure,
        fc="none", ec=edge, lw=1.4))
    fig.add_artist(FancyBboxPatch(
        ((X0[i] + 0.25) / W, (H - 0.76) / H), badge_w / W, 0.59 / H,
        boxstyle="round,pad=0,rounding_size=0.008", transform=fig.transFigure,
        fc=fill, ec="none"))
    base = (H - 0.60) / H
    fig.text((X0[i] + 0.41) / W, base, letter, fontsize=36.7,
             fontweight="bold", va="baseline")
    fig.text((X0[i] + 1.14) / W, base, title, fontsize=27.4,
             fontweight="semibold", va="baseline")


def axes3d(i):
    ax = padd(i, [0.13, 0.19, 0.74, 0.76], projection="3d")
    ax.set_axis_off()
    ax.set_box_aspect((12, 10, 9.5))
    ax.view_init(elev=ELEV, azim=AZIM)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.set_zlim(-ZMAX, 0)
    return ax


CORNERS = [(0, 0), (12, 0), (12, 10), (0, 10)]


def frame_box(ax):
    """The whole box; draw it before the contents."""
    for z in (0, -ZMAX):
        xs, ys = zip(*(CORNERS + CORNERS[:1]))
        ax.plot(xs, ys, [z] * 5, color="0.15", lw=1.0)
    for (x, y) in CORNERS:
        ax.plot([x, x], [y, y], [0, -ZMAX], color="0.15", lw=1.0)


def front_edges(ax):
    """The box edges nearest the viewer; draw them after the contents."""
    for z in (0, -ZMAX):
        ax.plot([0, 12, 12], [0, 0, 10], [z] * 3, color="0.15", lw=1.0)
    ax.plot([12, 12], [0, 0], [0, -ZMAX], color="0.15", lw=1.0)


def screen(ax, x, y, z):
    x2, y2, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
    return fig.transFigure.inverted().transform(ax.transData.transform((x2, y2)))


def box_right(ax):
    return max(screen(ax, x, y, z)[0] for x, y in CORNERS for z in (0, -ZMAX))


def depth_axis(i, ax):
    for d in (0, 200, 400, 600, 800):
        ax.text(0, 0, -d, f"{d}   ", ha="right", va="center", fontsize=12)
    p0, p6 = screen(ax, 0, 0, -50), screen(ax, 0, 0, -650)
    fig.text(p0[0] - dx(i, 0.15), (p0[1] + p6[1]) / 2, "Depth (m)", rotation=90,
             ha="center", va="center", fontsize=14.5, fontweight="bold")
    fig.add_artist(FancyArrowPatch(
        (p0[0] - dx(i, 0.105), p0[1]), (p0[0] - dx(i, 0.105), p6[1]),
        transform=fig.transFigure, arrowstyle="-|>", mutation_scale=18,
        lw=2.2, color="black"))


def range_labels(i, ax, x, y, upper, lower, labels):
    """Leaders from (x, y) at two depths to a labelled double arrow."""
    xr = box_right(ax) + dx(i, 0.01)
    pts = []
    for z, lab in zip((upper, lower), labels):
        s0 = screen(ax, x, y, -z)
        fig.add_artist(Line2D([s0[0], xr], [s0[1], s0[1]], color="0.1", lw=1.0,
                              transform=fig.transFigure))
        fig.add_artist(Ellipse((xr, s0[1]), dx(i, 0.016), 0.0142,
                               transform=fig.transFigure, color="black"))
        fig.text(xr + dx(i, 0.018), s0[1], lab, va="center", fontsize=14)
        pts.append((xr, s0[1]))
    fig.add_artist(FancyArrowPatch(
        pts[0], pts[1], transform=fig.transFigure, arrowstyle="<|-|>",
        mutation_scale=14, lw=1.8, color="black", shrinkA=7, shrinkB=7))


def map_axes(i, left, width, title):
    a = padd(i, [left, 0.0495, width, 0.145])  # bottom edge shared by all panels
    a.set_xlim(0, 12)
    a.set_ylim(0, 10)
    a.set_aspect("equal")
    a.set_xticks([])
    a.set_yticks([])
    for s in a.spines.values():
        s.set_color("0.6")
        s.set_lw(0.6)
    a.set_title(title, fontsize=13.5, fontweight="bold", pad=8)
    return a


def colorbar(i, a, cmap, norm, ticks):
    pos = a.get_position()
    cax = fig.add_axes([pos.x1 + dx(i, 0.012), pos.y0, dx(i, 0.018), pos.height])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax,
                      ticks=ticks)
    cb.ax.invert_yaxis()
    cb.ax.tick_params(labelsize=12, length=3)
    cb.outline.set_linewidth(0.4)


def grid_map(a, values, cmap, norm, outline=True):
    a.pcolormesh(xe, ye, np.ma.masked_invalid(values), cmap=cmap, norm=norm,
                 edgecolors="none")
    a.vlines(xe, 0, 10, color="0.82", lw=0.4)
    a.hlines(ye, 0, 12, color="0.82", lw=0.4)
    if outline:  # edge of the footprint cells
        for j in range(NY):
            for k in range(NX):
                if not inside[j, k]:
                    continue
                if k == 0 or not inside[j, k - 1]:
                    a.plot([xe[k]] * 2, ye[j:j + 2], color="0.1", lw=1.0)
                if k == NX - 1 or not inside[j, k + 1]:
                    a.plot([xe[k + 1]] * 2, ye[j:j + 2], color="0.1", lw=1.0)
                if j == 0 or not inside[j - 1, k]:
                    a.plot(xe[k:k + 2], [ye[j]] * 2, color="0.1", lw=1.0)
                if j == NY - 1 or not inside[j + 1, k]:
                    a.plot(xe[k:k + 2], [ye[j + 1]] * 2, color="0.1", lw=1.0)


# ---- panel A: SpatVector3D ------------------------------------------------------------
panel_frame(0, "A", "SpatVector3D", 3.55)
ax = axes3d(0)
ax.computed_zorder = False
ls = LightSource(azdeg=315, altdeg=50)


def surf(Xs, Ys, D):
    fc = ls.shade_rgb(reds(depth_norm(D))[..., :3], -D, vert_exag=0.02,
                      blend_mode="soft", fraction=0.35)
    ax.plot_surface(Xs, Ys, -D, facecolors=fc, rstride=1, cstride=1,
                    linewidth=0, antialiased=False, shade=False)


frame_box(ax)
Rg, Tg = np.meshgrid(np.linspace(0, 1, 40), th)
X, Y = to_xy(Rg, Tg)
Tw, Th = np.meshgrid(np.linspace(0, 1, 30), th)
WX, WY = to_xy(np.ones_like(Th), Th)
WZ = top_depth(WX, WY) + Tw * (bottom_depth(WX, WY) - top_depth(WX, WY))
surf(X, Y, bottom_depth(X, Y))
surf(WX, WY, WZ)
surf(X, Y, top_depth(X, Y))
# only the south-facing half of each outline is in view from this angle
front = (th >= np.pi) & (th <= 2 * np.pi)
ax.plot(bx[front], by[front], -top_depth(bx, by)[front], color="#7a0a0a",
        lw=1.6)
ax.plot(bx[front], by[front], -bottom_depth(bx, by)[front], color="#5a0005",
        lw=1.6)
for rr in (0.35, 0.7):  # contour rings show the top surface's curvature
    rx, ry = to_xy(np.full_like(th, rr), th)
    ax.plot(rx[front], ry[front], -top_depth(rx, ry)[front], color="#b23a36",
            lw=0.8, alpha=0.6)
front_edges(ax)
fig.canvas.draw()
depth_axis(0, ax)
px, py = to_xy(1.0, -0.35)
range_labels(0, ax, px, py, top_depth(px, py), bottom_depth(px, py),
             ("top", "bottom"))

a = map_axes(0, 0.03, 0.205, "Footprint")
a.fill(bx, by, color="#f4c9c7")
a.plot(bx[verts_idx], by[verts_idx], "o", ms=5, color="0.1")
a.plot(bx, by, color="0.1", lw=1.5)
gx, gy = np.meshgrid(np.linspace(0, 12, 400), np.linspace(0, 10, 340))
g_in = Path(np.c_[bx, by]).contains_points(
    np.c_[gx.ravel(), gy.ravel()]).reshape(gx.shape)
for left, title, field in ((0.335, "top (m)", top_depth),
                           (0.66, "bottom (m)", bottom_depth)):
    a = map_axes(0, left, 0.205, title)
    f = np.where(g_in, field(gx, gy), np.nan)
    a.contourf(gx, gy, f, levels=np.linspace(0, 600, 25), cmap=reds,
               norm=depth_norm)
    a.contour(gx, gy, f, levels=np.arange(50, 600, 50), colors="k",
              linewidths=0.3, alpha=0.35)
    a.plot(bx, by, color="0.1", lw=1.5)
    colorbar(0, a, reds, depth_norm, [0, 200, 400, 600])

# ---- shared 3D pieces for B and C --------------------------------------------------------
DZ = 50
ze = np.arange(0, ZMAX + DZ, DZ)       # layer edges (depth)
zmid = (ze[:-1] + ze[1:]) / 2
sf_xy = seafloor.T                                     # [x, y]


SHADE = {"top": 1.0, "south": 0.86, "east": 0.72}  # light from above-left


def voxel_faces(ax, filled, colours, edge):
    """Add the faces of `filled` voxels ([x, y, layer]) that face the viewer
    (top, south, east) and aren't covered by a neighbour, as one collection
    so matplotlib depth-sorts the faces."""
    nx, ny, nz = filled.shape
    quads, fcs = [], []
    for i, j, k in zip(*np.nonzero(filled)):
        x0, x1, y0, y1 = xe[i], xe[i + 1], ye[j], ye[j + 1]
        z0, z1 = -ze[k], -ze[k + 1]
        c = np.asarray(colours[i, j, k][:3])
        if k == 0 or not filled[i, j, k - 1]:
            quads.append([(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0)])
            fcs.append(c * SHADE["top"])
        if j == 0 or not filled[i, j - 1, k]:
            quads.append([(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)])
            fcs.append(c * SHADE["south"])
        if i == nx - 1 or not filled[i + 1, j, k]:
            quads.append([(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)])
            fcs.append(c * SHADE["east"])
    ax.add_collection3d(Poly3DCollection(quads, facecolors=fcs, edgecolors=edge,
                                         linewidths=0.35))


def seafloor_voxels(ax):
    filled = zmid[None, None, :] > sf_xy[..., None]
    colours = np.broadcast_to(greys(seafloor_norm(sf_xy))[..., None, :],
                              filled.shape + (4,))
    voxel_faces(ax, filled, colours, (0.3, 0.3, 0.3, 0.3))


# ---- panel B: SpatEnvelope -------------------------------------------------------------
panel_frame(1, "B", "SpatEnvelope", 3.7)
ax = axes3d(1)
ax.computed_zorder = False  # draw order: box, seafloor, envelope, front edges
frame_box(ax)
seafloor_voxels(ax)
dmin_xy, dmax_xy = depth_min.T, depth_max.T
filled = ((zmid[None, None, :] >= dmin_xy[..., None])
          & (zmid[None, None, :] <= dmax_xy[..., None]))
fc = np.broadcast_to(reds(depth_norm(zmid))[None, None], filled.shape + (4,))
voxel_faces(ax, filled, fc, (0.55, 0.1, 0.1, 0.15))
front_edges(ax)
fig.canvas.draw()
depth_axis(1, ax)
# label the east face of the easternmost cell in the middle row
j = NY // 2
k = np.nonzero(inside[j])[0].max()
z_lo = zmid[(zmid >= depth_min[j, k]) & (zmid <= depth_max[j, k])]
range_labels(1, ax, xe[k + 1], yc[j, k], z_lo.min() - DZ / 2,
             z_lo.max() + DZ / 2, ("depth_min", "depth_max"))

for left, title, vals, cmap, norm, ticks in (
        (0.03, "depth_min (m)", depth_min, reds, depth_norm, [0, 200, 400, 600]),
        (0.335, "depth_max (m)", depth_max, reds, depth_norm, [0, 200, 400, 600]),
        (0.66, "Seafloor (m)", seafloor, greys, seafloor_norm,
         [0, 200, 400, 600, 800])):
    a = map_axes(1, left, 0.205, title)
    grid_map(a, vals, cmap, norm)
    colorbar(1, a, cmap, norm, ticks)

# ---- panel C: SpatVoxel -------------------------------------------------------------------
panel_frame(2, "C", "SpatVoxel", 3.15)
ax = axes3d(2)
ax.computed_zorder = False  # box, seafloor, planes deepest first, front edges
frame_box(ax)
seafloor_voxels(ax)
for z in sorted(STD_DEPTHS, reverse=True):
    quads, colours = [], []
    for jj in range(NY):
        for kk in range(NX):
            quads.append([(xe[kk], ye[jj], -z), (xe[kk + 1], ye[jj], -z),
                          (xe[kk + 1], ye[jj + 1], -z), (xe[kk], ye[jj + 1], -z)])
            colours.append(matplotlib.colors.to_rgba(PRESENT, 0.92)
                           if presence[z][jj, kk]
                           else matplotlib.colors.to_rgba(ABSENT, 0.35))
    ax.add_collection3d(Poly3DCollection(quads, facecolors=colours,
                                         edgecolors=(0.6, 0.6, 0.6, 0.5),
                                         linewidths=0.4))
front_edges(ax)
fig.canvas.draw()
xr = box_right(ax) + dx(2, 0.01)
for z in STD_DEPTHS:
    s0 = screen(ax, 12, 10, -z)
    fig.add_artist(Line2D([s0[0], xr], [s0[1], s0[1]], color="0.3", lw=0.8,
                          transform=fig.transFigure))
    fig.text(xr + dx(2, 0.012), s0[1], f"z = {z} m", va="center", fontsize=14)

presence_cmap = ListedColormap([ABSENT, PRESENT])
for n, z in enumerate(STD_DEPTHS):
    a = map_axes(2, 0.03 + n * 0.152, 0.138, f"z = {z} m")
    grid_map(a, presence[z].astype(float), presence_cmap, Normalize(0, 1),
             outline=False)
for n, (col, lab) in enumerate(((PRESENT, "Species\npresence"),
                                (ABSENT, "Absence"))):
    y = 0.155 - n * 0.06
    x = pfig(2, 0.815, y)[0]
    fig.add_artist(Rectangle((x, y - 0.012), dx(2, 0.03), 0.024,
                             transform=fig.transFigure, fc=col, ec="0.5", lw=0.6))
    fig.text(x + dx(2, 0.045), y, lab, va="center", fontsize=12)

# ---- conversion arrows ------------------------------------------------------------------------
if not args.no_arrows:
    for i, label in ((0, "vector3d_to_envelope()"), (1, "envelope_to_voxel()")):
        x_a, x_b = (X0[i] + PW + 0.24) / W, (X0[i + 1] - 0.18) / W
        y = 1 - 422 / 887
        fig.text((x_a + x_b) / 2, y + 0.04, label, ha="center", va="center",
                 fontsize=12.2, fontweight="bold")
        fig.add_artist(FancyArrowPatch((x_a, y), (x_b, y),
                                       transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=22,
                                       lw=2.2, color="black"))

fig.savefig(args.out, dpi=args.dpi, facecolor="white")
print(args.out, fig.get_size_inches() * args.dpi)
