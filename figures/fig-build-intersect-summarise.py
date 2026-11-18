"""Builds an overview of the package design in SPEC.md: an analysis is
build -> intersect -> summarise, and each function belongs to one step.

  1 Build      2D inputs (a range polygon with a depth range, observations
               with depths) become 3D objects: a species envelope and 3D
               buffers around the observations
  2 Intersect  where a species range and a fishing layer share 3D space;
               intersects_3d() per cell and mask() on a temperature voxel
  3 Summarise  volume, overlap, a temperature profile inside the range and
               per-buffer statistics

All panels use one synthetic example (16 x 12 cells of 10 km, depths in m),
so the numbers in panel 3 are computed from the shapes in panels 1 and 2.
Run from the project root:

  python3 figures/fig-build-intersect-summarise.py \
    vignettes/figures/Build_Intersect_Summarise.png

Needs numpy, matplotlib and the Noto Sans fonts.
"""
import argparse

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import (LinearSegmentedColormap, ListedColormap,  # noqa: E402
                               Normalize, to_rgba)
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("out")
parser.add_argument("--dpi", type=int, default=200)
args = parser.parse_args()

FONT_DIR = "/usr/share/fonts/truetype/noto/"
for f in ("Regular", "Bold", "SemiBold", "Italic"):
    font_manager.fontManager.addfont(f"{FONT_DIR}NotoSans-{f}.ttf")
plt.rcParams["font.family"] = "Noto Sans"
MONO = "DejaVu Sans Mono"

# ---- colours ------------------------------------------------------------------
SPECIES, FISHING, SHARED = "#ef6f69", "#5b8fd6", "#7b4fb0"
reds = LinearSegmentedColormap.from_list(
    "depth", ["#fff4f3", "#f7b9b6", "#ec7a74", "#d23a35", "#9c0b10", "#6e0007"])
depth_norm = Normalize(0, 600)
greys = LinearSegmentedColormap.from_list(
    "seafloor", ["#ffffff", "#e3e5e7", "#b9bec3", "#8d959e", "#5c6875"])
seafloor_norm = Normalize(0, 800)
thermal = LinearSegmentedColormap.from_list(
    "temp", ["#2c3e8f", "#4d8cc4", "#9fd0d8", "#f3e7a6", "#f2a35e", "#c8482f"])
temp_norm = Normalize(4, 18)
STEPS = {  # badge fill, border
    "1": ("#dcefe3", "#b9dcc6"),
    "2": ("#e9def5", "#d2c0ec"),
    "3": ("#fdebd0", "#f3d3a3"),
}
PLANNED = "0.5"

# ---- synthetic example ---------------------------------------------------------
ZMAX, CELL_KM = 800, 10
NX, NY = 16, 12
xe, ye = np.linspace(0, 12, NX + 1), np.linspace(0, 10, NY + 1)
xc, yc = np.meshgrid((xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2)  # [ny, nx]
DZ = 50
ze = np.arange(0, ZMAX + DZ, DZ)
zmid = (ze[:-1] + ze[1:]) / 2

# species range: a polygon with a top and bottom surface (as in the
# SpatVector3D figure), put on the grid and capped by the seafloor
cx, cy, ax_x, ax_y = 5.2, 5.6, 3.9, 3.0
th = np.linspace(0, 2 * np.pi, 181)


def radius(t):
    return 1 + 0.10 * np.sin(3 * t + 0.4) + 0.07 * np.cos(2 * t - 0.8)


bx = cx + ax_x * radius(th) * np.cos(th)
by = cy + ax_y * radius(th) * np.sin(th)
d2 = ((xc - cx) / ax_x) ** 2 + ((yc - cy) / ax_y) ** 2
seafloor = np.round(np.clip(
    640 - 110 * np.exp(-((xc - 2) ** 2 + (yc - 7.5) ** 2) / 5)
    + 120 * np.exp(-((xc - 9.5) ** 2 + (yc - 1.5) ** 2) / 9)
    + 25 * np.sin(1.3 * xc) * np.cos(0.9 * yc), 570, 790), -1)
sp_in = Path(np.c_[bx, by]).contains_points(
    np.c_[xc.ravel(), yc.ravel()]).reshape(xc.shape)
sp_min = np.where(sp_in, np.round(25 + 135 * np.clip(d2, 0, 1.4) / 1.2), np.nan)
sp_max = np.where(sp_in, np.minimum(
    np.round(250 + 24 * (yc - 1) + 180 * np.exp(-2 * d2)), seafloor), np.nan)

# fishing layer: a midwater fishery at 200-400 m over an ellipse to the east
fi_in = ((xc - 8.6) / 2.9) ** 2 + ((yc - 3.6) / 2.2) ** 2 <= 1
fi_min = np.where(fi_in, 200.0, np.nan)
fi_max = np.where(fi_in, 400.0, np.nan)

# intersect_3d(): the depth interval both share in each cell
sh_min, sh_max = np.maximum(sp_min, fi_min), np.minimum(sp_max, fi_max)  # NaN unless both
shared = sh_max > sh_min
intersects = np.where(sp_in & fi_in, shared, np.where(sp_in | fi_in, False, np.nan))


def volume_km3(dmin, dmax):
    """Sum of cell area x (depth_max - depth_min), as volume() does."""
    return np.nansum(np.where(dmax > dmin, dmax - dmin, 0)) / 1000 * CELL_KM ** 2


v_species = volume_km3(sp_min, sp_max)
v_fishing = volume_km3(fi_min, fi_max)
v_shared = volume_km3(sh_min, sh_max)

# a temperature voxel on the same grid, masked by the species range
temp = (3.5 + 14 * np.exp(-zmid[:, None, None] / 170)
        + 0.35 * xc[None] - 0.3 * yc[None])                  # [z, ny, nx]
in_range = ((zmid[:, None, None] >= sp_min[None])
            & (zmid[:, None, None] <= sp_max[None]))
masked = np.where(in_range, temp, np.nan)
level_n = in_range.sum(axis=(1, 2))
profile_mean = np.nanmean(np.where(level_n[:, None, None] > 0, masked, 0), axis=(1, 2))
profile_sd = np.nanstd(np.where(level_n[:, None, None] > 0, masked, 0), axis=(1, 2))
has_level = level_n > 0

# observations with depths, each buffered (radius 1 cell, +/- 50 m)
obs = [(1.6, 1.4, 120), (10.6, 8.6, 260), (11.0, 1.2, 450)]
BUF_R, BUF_DZ = 0.9, 50
buffers = []
for k, (ox, oy, oz) in enumerate(obs, start=1):
    near = (xc - ox) ** 2 + (yc - oy) ** 2 <= (BUF_R + 0.45) ** 2
    lev = np.abs(zmid - oz) <= BUF_DZ
    vals = temp[lev][:, near]
    buffers.append((k, vals.mean(), vals.std(), vals.size))

# ---- figure layout -----------------------------------------------------------------
PW, H, GAP = 7.9, 10.6, 1.0
X0 = [0.0, PW + GAP, 2 * (PW + GAP)]
W = X0[-1] + PW
fig = plt.figure(figsize=(W, H))


def pfig(i, x, y):
    return ((X0[i] + x * PW) / W, y)


def padd(i, rect, **kw):
    x, y = pfig(i, rect[0], rect[1])
    return fig.add_axes([x, y, rect[2] * PW / W, rect[3]], **kw)


def ptext(i, x, y, s, **kw):
    return fig.text(*pfig(i, x, y), s, **kw)


def panel_frame(i, number, title, question, badge_w):
    fill, edge = STEPS[number]
    fig.add_artist(FancyBboxPatch(
        ((X0[i] + 0.12) / W, 0.15 / H), (PW - 0.24) / W, (H - 0.19) / H,
        boxstyle="round,pad=0,rounding_size=0.008", transform=fig.transFigure,
        fc="none", ec=edge, lw=1.4))
    fig.add_artist(FancyBboxPatch(
        ((X0[i] + 0.25) / W, (H - 0.76) / H), badge_w / W, 0.59 / H,
        boxstyle="round,pad=0,rounding_size=0.008", transform=fig.transFigure,
        fc=fill, ec="none"))
    base = (H - 0.60) / H
    fig.text((X0[i] + 0.41) / W, base, number, fontsize=36.7, fontweight="bold",
             va="baseline")
    fig.text((X0[i] + 1.0) / W, base, title, fontsize=27.4, fontweight="semibold",
             va="baseline")
    fig.text((X0[i] + 0.3) / W, (H - 1.15) / H, question, fontsize=16,
             style="italic", color="0.25", va="center")


def function_list(i, now, planned, top=0.235):
    """The step's functions: existing ones in black, planned ones in grey."""
    fig.add_artist(Rectangle(pfig(i, 0.04, 0.03), 0.92 * PW / W, top - 0.03,
                             transform=fig.transFigure, fc="#fafafa", ec="0.85",
                             lw=0.8))
    y = top - 0.03
    for heading, names, colour in (("Exists now", now, "black"),
                                   ("Planned", planned, PLANNED)):
        if not names:
            continue
        ptext(i, 0.07, y, heading, fontsize=13, fontweight="bold", color=colour,
              va="center")
        y -= 0.03
        for row in names:
            ptext(i, 0.09, y, row, fontsize=12.5, family=MONO, color=colour,
                  va="center")
            y -= 0.026
        y -= 0.008


def map_axes(i, rect, title):
    a = padd(i, rect)
    a.set_xlim(0, 12)
    a.set_ylim(0, 10)
    a.set_aspect("equal")
    a.set_xticks([])
    a.set_yticks([])
    for s in a.spines.values():
        s.set_color("0.6")
        s.set_lw(0.6)
    a.set_title(title, fontsize=13.5, fontweight="bold", pad=6)
    return a


def grid_lines(a):
    a.vlines(xe, 0, 10, color="0.82", lw=0.4)
    a.hlines(ye, 0, 12, color="0.82", lw=0.4)


def colorbar(i, a, cmap, norm, ticks, label=None):
    pos = a.get_position()
    cax = fig.add_axes([pos.x1 + 0.012 * PW / W, pos.y0, 0.018 * PW / W,
                        pos.height])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax,
                      ticks=ticks)
    cb.ax.tick_params(labelsize=12, length=3)
    cb.outline.set_linewidth(0.4)
    if label:
        cb.set_label(label, fontsize=12)
    return cb


# ---- 3D helpers -----------------------------------------------------------------------
CORNERS = [(0, 0), (12, 0), (12, 10), (0, 10)]
SHADE = {"top": 1.0, "south": 0.86, "east": 0.72}


def axes3d(i, rect, elev=18, azim=-64):
    ax = padd(i, rect, projection="3d")
    ax.set_axis_off()
    ax.set_box_aspect((12, 10, 9.5))
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.set_zlim(-ZMAX, 0)
    ax.computed_zorder = False  # draw in the order added
    for z in (0, -ZMAX):
        xs, ys = zip(*(CORNERS + CORNERS[:1]))
        ax.plot(xs, ys, [z] * 5, color="0.15", lw=1.0)
    for (x, y) in CORNERS:
        ax.plot([x, x], [y, y], [0, -ZMAX], color="0.15", lw=1.0)
    return ax


def front_edges(ax):
    for z in (0, -ZMAX):
        ax.plot([0, 12, 12], [0, 0, 10], [z] * 3, color="0.15", lw=1.0)
    ax.plot([12, 12], [0, 0], [0, -ZMAX], color="0.15", lw=1.0)


def layers(dmin, dmax):
    """Voxels [x, y, layer] whose layer midpoint is inside [dmin, dmax]."""
    return ((zmid[None, None] >= dmin.T[..., None])
            & (zmid[None, None] <= dmax.T[..., None]))


def voxel_faces(ax, label, colours, edge, lw=0.35):
    """Faces of labelled voxels ([x, y, layer], 0 = empty) that face the viewer
    and border a voxel with a different label, as one depth-sorted collection.
    `colours` maps label -> RGBA, or is an array of RGBA per voxel."""
    nx, ny, nz = label.shape
    quads, fcs = [], []
    for i, j, k in zip(*np.nonzero(label)):
        lab = label[i, j, k]
        c = np.asarray(colours[lab] if isinstance(colours, dict)
                       else colours[i, j, k], dtype=float)
        x0, x1, y0, y1, z0, z1 = xe[i], xe[i + 1], ye[j], ye[j + 1], -ze[k], -ze[k + 1]
        if k == 0 or label[i, j, k - 1] != lab:
            quads.append([(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0)])
            fcs.append(np.r_[c[:3] * SHADE["top"], c[3]])
        if j == 0 or label[i, j - 1, k] != lab:
            quads.append([(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)])
            fcs.append(np.r_[c[:3] * SHADE["south"], c[3]])
        if i == nx - 1 or label[i + 1, j, k] != lab:
            quads.append([(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)])
            fcs.append(np.r_[c[:3] * SHADE["east"], c[3]])
    ax.add_collection3d(Poly3DCollection(quads, facecolors=fcs, edgecolors=edge,
                                         linewidths=lw))


def seafloor_voxels(ax):
    filled = zmid[None, None] > seafloor.T[..., None]
    colours = np.broadcast_to(greys(seafloor_norm(seafloor.T))[..., None, :],
                              filled.shape + (4,))
    voxel_faces(ax, filled.astype(int), colours, (0.3, 0.3, 0.3, 0.3))


def cylinder(ax, x, y, z, r, dz, colour):
    t = np.linspace(0, 2 * np.pi, 60)
    h = np.array([z - dz, z + dz])
    T, Hh = np.meshgrid(t, h)
    ax.plot_surface(x + r * np.cos(T), y + r * np.sin(T), -Hh, color=colour,
                    alpha=0.35, linewidth=0, shade=False)
    for zz in (z - dz, z + dz):
        ax.plot(x + r * np.cos(t), y + r * np.sin(t), -zz, color=colour, lw=1)


def depth_ticks(ax):
    for d in (0, 200, 400, 600, 800):
        ax.text(0, 0, -d, f"{d}   ", ha="right", va="center", fontsize=12)


def legend(i, x, y, items, dy=0.028):
    for n, (colour, lab) in enumerate(items):
        yy = y - n * dy
        x0 = pfig(i, x, yy)[0]
        fig.add_artist(Rectangle((x0, yy - 0.008), 0.025 * PW / W, 0.016,
                                 transform=fig.transFigure, fc=colour, ec="0.5",
                                 lw=0.6))
        ptext(i, x + 0.04, yy, lab, fontsize=12.5, va="center")


# ---- panel 1: Build -------------------------------------------------------------------------
panel_frame(0, "1", "Build", "What is the 3D shape or grid?", 2.45)
ptext(0, 0.06, 0.845, "2D inputs", fontsize=14, fontweight="bold")
a = map_axes(0, [0.06, 0.645, 0.24, 0.135], "polygon + depths")
a.fill(bx, by, color="#f4c9c7")
a.plot(bx, by, color="0.1", lw=1.3)
a = map_axes(0, [0.06, 0.445, 0.24, 0.135], "points + depths")
for (ox, oy, oz) in obs:
    a.plot(ox, oy, "o", ms=7, color="0.1")
    a.text(ox, oy - 1.6 if oy > 5 else oy + 1.0, f"{oz} m", ha="center",
           fontsize=11)
fig.add_artist(FancyArrowPatch(pfig(0, 0.33, 0.6), pfig(0, 0.4, 0.6),
                               transform=fig.transFigure, arrowstyle="-|>",
                               mutation_scale=20, lw=2, color="black"))
ptext(0, 0.43, 0.845, "3D objects", fontsize=14, fontweight="bold")
ax = axes3d(0, [0.4, 0.4, 0.58, 0.42])
seafloor_voxels(ax)
sp_vox = layers(sp_min, sp_max)
voxel_faces(ax, sp_vox.astype(int),
            np.broadcast_to(reds(depth_norm(zmid))[None, None], sp_vox.shape + (4,)),
            (0.55, 0.1, 0.1, 0.15))
for (ox, oy, oz) in obs:
    cylinder(ax, ox, oy, oz, BUF_R, BUF_DZ, "#3a6fb8")
    ax.scatter([ox], [oy], [-oz], color="black", s=14, depthshade=False)
front_edges(ax)
legend(0, 0.44, 0.39, [("#f4a19c", "species envelope (SpatEnvelope)"),
                       ("#a9c1e3", "3D buffers (SpatVector3D, planned)")])
function_list(
    0,
    ["as_voxel()  as_envelope()  vect_to_envelope()",
     "envelope_to_voxel()  voxel_to_envelope()",
     "occupied()  create_study_raster()",
     "woa_load_nc()  load_gebco_bathymetry()",
     "gfw_effort_to_raster()"],
    ["as_vector3d()  buffer_3d()",
     "vector3d_to_envelope()  track sampling"],
    top=0.32)

# ---- panel 2: Intersect -----------------------------------------------------------------------
panel_frame(1, "2", "Intersect", "Where do they meet?", 3.25)
ax = axes3d(1, [0.12, 0.47, 0.78, 0.40])
depth_ticks(ax)
seafloor_voxels(ax)
fi_vox = layers(fi_min, fi_max)
label = np.zeros(sp_vox.shape, int)
label[sp_vox] = 1
label[fi_vox] = 2
label[sp_vox & fi_vox] = 3
voxel_faces(ax, label, {1: to_rgba(SPECIES, 0.2), 2: to_rgba(FISHING, 0.4),
                        3: to_rgba(SHARED, 1.0)}, (0.2, 0.2, 0.2, 0.12))
front_edges(ax)
legend(1, 0.08, 0.515, [(SPECIES, "species range"), (FISHING, "fishing layer"),
                        (SHARED, "intersect_3d(): space they share")])
a = map_axes(1, [0.07, 0.235 + 0.03, 0.27, 0.12], "intersects_3d()")
a.pcolormesh(xe, ye, np.ma.masked_invalid(intersects.astype(float)),
             cmap=ListedColormap(["#d9d9d9", SHARED]), vmin=0, vmax=1)
grid_lines(a)
ptext(1, 0.36, 0.33, "TRUE", fontsize=12, color=SHARED, fontweight="bold",
      va="center")
ptext(1, 0.36, 0.305, "FALSE", fontsize=12, color="0.45", fontweight="bold",
      va="center")
ptext(1, 0.36, 0.28, "NA (blank)", fontsize=12, color="0.45", va="center")
lev100 = np.searchsorted(zmid, 100)
a = map_axes(1, [0.55, 0.265, 0.27, 0.12], "mask(temp, range)")
a.pcolormesh(xe, ye, np.ma.masked_invalid(masked[lev100]), cmap=thermal,
             norm=temp_norm)
grid_lines(a)
a.plot(bx, by, color="0.1", lw=1)
colorbar(1, a, thermal, temp_norm, [5, 10, 15], label="°C")
ptext(1, 0.69, 0.245, "at 100 m", fontsize=11.5, color="0.35", ha="center")
function_list(1, ["intersect_3d()  intersects_3d()  mask()"],
              ["the same three for SpatVector3D",
               "coverage_3d(): data cells in each shape,",
               "  with weights"],
              top=0.22)

# ---- panel 3: Summarise --------------------------------------------------------------------------
panel_frame(2, "3", "Summarise", "How much, and what values?", 3.6)


def card(i, x, y, w, h, call, value, note):
    fig.add_artist(FancyBboxPatch(pfig(i, x, y), w * PW / W, h,
                                  boxstyle="round,pad=0,rounding_size=0.006",
                                  transform=fig.transFigure, fc="white", ec="0.75",
                                  lw=1))
    ptext(i, x + w / 2, y + h - 0.03, call, fontsize=12.5, family=MONO,
          ha="center", va="center")
    ptext(i, x + w / 2, y + h / 2 - 0.008, value, fontsize=24, fontweight="bold",
          ha="center", va="center")
    ptext(i, x + w / 2, y + 0.022, note, fontsize=11.5, color="0.35",
          ha="center", va="center")


card(2, 0.06, 0.725, 0.42, 0.13, "volume(range)", f"{v_species:,.0f} km³",
     "species range")
card(2, 0.52, 0.725, 0.42, 0.13, "calc_volume_overlap()",
     f"{100 * v_shared / v_species:.0f} %",
     f"of the range is fished ({v_shared:,.0f} km³)")

a = padd(2, [0.15, 0.35, 0.36, 0.30])
zz, mm, ss = zmid[has_level], profile_mean[has_level], profile_sd[has_level]
a.fill_betweenx(zz, mm - ss, mm + ss, color=SPECIES, alpha=0.25, lw=0)
a.plot(mm, zz, "-o", color="#b8322c", ms=4, lw=1.8)
a.set_ylim(ZMAX * 0.75, 0)
a.set_xlim(3, 18)
a.set_xlabel("Temperature (°C)", fontsize=12.5)
a.set_ylabel("Depth (m)", fontsize=12.5)
a.tick_params(labelsize=11.5)
for s in ("top", "right"):
    a.spines[s].set_visible(False)
a.set_title("mean ± sd per depth\nin the masked range", fontsize=13.5,
            fontweight="bold", pad=8)

a = padd(2, [0.6, 0.35, 0.36, 0.30])
a.set_axis_off()
a.set_title("summarise_3d(), per buffer", fontsize=13.5, fontweight="bold",
            color=PLANNED, pad=8)
rows = [["id", "mean °C", "sd", "cells"]] + [
    [str(k), f"{m:.1f}", f"{s:.1f}", str(n)] for k, m, s, n in buffers]
tab = a.table(cellText=rows, loc="upper center", cellLoc="center",
              colWidths=[0.17, 0.33, 0.2, 0.25])
tab.auto_set_font_size(False)
tab.set_fontsize(12.5)
tab.scale(1, 2.0)
for (r, c), cell in tab.get_celld().items():
    cell.set_edgecolor("0.8")
    cell.set_text_props(color=PLANNED if r else "0.3",
                        fontweight="bold" if r == 0 else "normal")
    if r == 0:
        cell.set_facecolor("#f2f2f2")
function_list(2, ["volume()  calc_volume_overlap()"],
              ["summarise_3d(): statistics on intersect",
               "  results, per shape or group"],
              top=0.22)

# ---- arrows between the steps ---------------------------------------------------------------------
for i in (0, 1):
    y = 0.62
    fig.add_artist(FancyArrowPatch(((X0[i] + PW - 0.05) / W, y),
                                   ((X0[i + 1] + 0.05) / W, y),
                                   transform=fig.transFigure, arrowstyle="-|>",
                                   mutation_scale=26, lw=2.6, color="black"))

fig.savefig(args.out, dpi=args.dpi, facecolor="white")
print(args.out, fig.get_size_inches() * args.dpi,
      f"species {v_species:.0f} km3, fishing {v_fishing:.0f}, shared {v_shared:.0f}")
