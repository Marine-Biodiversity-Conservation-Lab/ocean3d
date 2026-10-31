# **Project Spec**

This document describes the overall design decisions for the `ocean3d` project. `ocean3d` is an R package that is designed to facilitate the analysis of marine habitats in 3D, enabling descriptions of interaction between different marine entities, including species, anthropogenic activities, and environmental variables. `ocean3d` is intended to empower users to conduct analyses in the 3D marine environment within the R ecosystem. 

The primary operation that `ocean3d` is focused on is the spatial querying of large 3D marine datasets. `ocean3d` enables asking questions like 
1. What range of environmental conditions is a species found in? 
2. How much does a given anthropogenic activity affect a species? 
3. What impact could protective measures have on a species?

while considering three-dimensional space. 

![Conceptual overview: a species' 2D range, depth range and bathymetry combine into a 3D distribution, which is then queried against environmental data](vignettes/figures/Figure1_conceptual3Ddistribution.png)

*Figure 1. A species' geographical range (A) and depth range, bounded by bathymetry, give its potential 3D distribution (B). Environmental data are then summarised within that space (C). From the Getting Started vignette.*

Technically, this means that `ocean3d` needs to be able to take 1D (points), 2D (polygons), and 3D (raster stack) and query how these intersect with one another. 

## Vector Representations in 3D space
For more background: https://en.wikipedia.org/wiki/Geometric_primitive
### 0D (Points / Vertices)
These are 0-dimensional, without length, width, or height dimensions. Points can be located in 3-dimensional space though, with 3 coordinate values. They are called vertices in the 3D graphics / modelling field. 

For `ocean3d`, these will be lon (longitude), lat (latitude), and depth (metres). 

### 1D (Lines / Edges)
Lines are two points connected, having a length but no width, therefore considered to be 1-dimensional. These are called edges in the 3D graphics / modelling field. 

`extract_to_track()` extracts environmental values along lines, such as tagging tracks or fishing vessel routes, by sampling points along each ordered track and extracting at them with `extract_to_point()` (#31). Lines are not yet a 3D shape of their own; see [Data model](#data-model).

### 2D (Polygons / Faces)
Polygons are constructed from points that exist on the same 2D plane, that are connected in order that creates a closed shape. With length and width, these are 2-dimensional shapes. In 3D graphics / modelling field, they are called faces. 

Currently, `ocean3d` doesn't handle faces that are natively 3D. These are cases where each point that is part of the polygon has a meaningful depth coordinate. Instead, we work with polygons that have depth values associated with the entire polygon. I.e. the polygon is represented on planes that aligned with depth levels. For example, IUCN Red List species ranges that have a species range represented as a 2D polygon with no depths, with separate depth range values provided. This is often called 2.5D representation. 

### 3D (Volumes)
Volumes are closed geometries that are formed by multiple polygons. 

Currently, `ocean3d` does not handle volumes. This might become useful to work with if we input 3D models, derived from techniques like photogrammetry. 

## Raster Representations in 3D space

### 2D (Raster) 

2D rasters are cellular representations with x, y dimensions. I.e. the whole study space is represented as a grid that has values in every cell. 

### 3D (Voxel)

![The same species distribution as a SpatEnvelope (depth_min and depth_max rasters) and as a SpatVoxel (presence at standard depth levels), converted with envelope_to_voxel()](vignettes/figures/Comparison_SpatEnvelope_SpatVoxel.png)

*Figure 2. The two raster representations of the same distribution. (A) A `SpatEnvelope` ([min-max 2.5D](#3d-min-max-25d)) stores `depth_min` and `depth_max` per cell, bounded by the seafloor. (B) `envelope_to_voxel()` converts it to a `SpatVoxel`, one layer per standard depth. From the Getting Started vignette.*

To extend rasters into 3D, we simply create multiple rasters for a set of standard depths. This is the convention from oceanography datasets, like Copernicus Marine and World Ocean Atlas (WOA). This raster stack is supported with netCDF file type and with the `terra` and `ncdf4` packages. This raster stack is called a **voxel model** in 3D graphics / modelling field. 

Multi-depth rasters used by `ocean3d` must encode depth in layer names using the format `{variable}_{field}_depth={value}` (e.g., `t_an_depth=0`, `t_an_depth=100` for temperature with annual means at 0 and 100 m depths). 

The convention for the depth bounds that a given raster layer represents depends on source netCDF. For example, WOA provides standard depth levels and explicit depth_bnds. For example, the 0 m layer has bounds from 0 to 2.5 m, the 5 m layer 2.5–7.5 m, the 10 m layer 7.5–12.5 m, etc. 

In these rasters, there is a continuous representation of values across x, y, and depth dimensions. For environmental variables, this could be temperature in degrees C at each coordinate and depth. For a species range, this could be probability of occurrence at each coordinate an depth. 

Data source utilities are responsible for converting other formats into this convention. Layer names are parsed by `depths()`, which utilities such as `extract_to_area()` use to determine which depth layers to select for a given depth range.

### 3D (Min-max 2.5D)

Colloquially referred to as the 2.5D model, this is often used for terrain models in GIS. Think like Digital Elevation Models, where the raster values correspond with the height of the landscape. 

We can represent spaces, like species ranges, by describing the space that they occupy with a minimum depth and maximum depth raster (Figure 2A). Maximum depth could be bound by the bathymetry raster, or by the maximum depth the species is found at. 

This representation is more efficient than the 3D Voxel representation, since it only needs two rasters for min and max respectively, compared to needing a raster for every standard depth. The downside is that the 2.5 representation doesn't model situations where there may not be presences between the min and max depths. 

## Data model

### Now: two models

`ocean3d` currently has two separate data models for the same kind of data (gridded marine variables by longitude, latitude and depth). Functions take one or the other, and there is no bridge between them.

| | netCDF model | `SpatVoxel` model |
|---|---|---|
| Functions | `extract_to_point()` and its `extract2d()` / `extract3d_*()` wrappers, `extract_to_buffer()`, `extract_to_track()` | `extract_to_area()`, `volume()`, `calc_volume_overlap()`, `intersect_3d()`, `intersects_3d()`, `mask()`, voxel profiles, `as_voxel()`, `as_envelope()`, `woa_load_nc()`, `load_gebco_bathymetry()`, `gfw_effort_to_raster()` |
| Input | netCDF paths, connections or file tables | `SpatVoxel` / `SpatEnvelope`, S4 classes extending `terra::SpatRaster` |
| Reads with | ncdf4, only the cells each observation needs | terra |
| Depth | the netCDF depth coordinate: units converted to metres, `positive` direction, bounds | layer names (`{variable}_depth={value}`), one set of depths for the whole grid |
| Time | CF time decoding (calendars, sub-daily steps, nearest or exact matching, series split across files) | none |
| Cell edges | CF bounds, or inferred from cell centres | terra's regular grid extent and resolution |

The netCDF model keeps 4D data (lon, lat, depth, time) and its CF metadata intact and reads only the cells it needs. The voxel model gets terra's raster operations (crop, mask, cell area, resampling), which the volume and overlap tools depend on. The cost of the split: a voxel from `woa_load_nc()` can't be passed to `extract_to_point()`, a netCDF file can't be passed to `extract_to_area()`, depth and grid handling is written twice with different rules, and the voxel tools have no time dimension.

### Planned: one model (#99)

This section is the plan from #99. None of it is implemented yet, and the first step (the stars trial below) may change it.

**Requirements.**

1. **Time-aware 3D analysis.** Time is a real dimension of the core object. Volume, overlap and profiles work through time (for example seasonal habitat volume, or overlap by month), not only on a static voxel.
2. **Lazy reads for data larger than memory.** Opening a dataset doesn't read its values, every operation reads only the slices it needs, and heavy processing runs per time step or per chunk, so peak memory follows the size of one slice rather than the whole dataset.
3. **CF-correct metadata.** Keep the current handling of time units and calendars, sub-daily steps, depth units and direction, bounds and longitude conventions (#71–#75, #88–#96).

No single package meets all three. terra stores depth and time as layer labels on a stack of 2D layers. stars has real dimensions and lazy `stars_proxy` objects, but only some operations stay lazy. ncdf4 reads exactly what is asked for but leaves everything above reading to the package.

**Direction: a hybrid.**

- **Core object:** one object that describes the dimensions (lon/lat with bounds, depth, time, CRS) without loading values, in the style of stars. It may be a `stars_proxy` itself if the trial supports that.
- **Reads:** lazy and exact, through stars proxies, or through the package's own ncdf4 `start` / `count` code where stars falls short. Remote Zarr stores (for example Copernicus Marine ARCO data) would be read through GDAL's multidimensional API (`stars::read_mdim()`).
- **Heavy spatial work:** per time step or per chunk. Read one slice, process it with terra, and combine the results across time.
- **`SpatVoxel` / `SpatEnvelope`:** become views of one time step of the core object, so the existing terra-based tools keep working on a slice.
- **Extraction:** point, buffer and track extraction keep reading only the cells they touch, and accept anything that can become the core object: a netCDF file, a `SpatVoxel`, or a Zarr store.

**First step: stars trial.** Before any code moves, compare stars with the current ncdf4 code:

1. On the synthetic netCDF fixtures from the tests and the #88–#96 cases: do decoded time, depth, bounds and values match the package's results?
2. On a large real file (for example a multi-year daily global Copernicus Marine product), measure peak memory and run time for opening as a proxy, point extraction across time, volume and overlap for one time step, and a polygon or buffer summary over a time window.
3. Optionally, the same dataset as a remote Zarr store.

The result decides whether stars proxies become the core object, or only describe dimensions while the package keeps its own lazy reads.

**Extraction on regions.** Point, buffer, area and track extraction are all "summarise the source over a 3D (or 4D) region". On top of the core object:

```r
# one horizontal geometry (point, line or polygon) + one depth range + one time range
region <- region_3d(geometry, depth = c(0, 50), time = c(t0, t1))

extract_to_region(source, region, statistic = "mean", weighting = "weighted")
```

A point is the zero-radius case, an area is a polygon, and a track is a buffered line. `extract_to_point()`, `extract_to_buffer()`, `extract_to_area()` and `extract_to_track()` stay as wrappers with their current arguments and results. The first version holds one geometry per region; collections of regions (one per observation or track segment), summarised in batches so each block of the source is read once, come later.

**Vector 3D shape class: `SpatSolid`.** The region is a new class alongside `SpatVoxel` and `SpatEnvelope`. `SpatEnvelope` is "footprint + top + bottom" on a grid; `SpatSolid` is the vector version:

| | `SpatVoxel` | `SpatEnvelope` | `SpatSolid` (planned) |
|---|---|---|---|
| Holds | values of a variable | a shape | a shape |
| Horizontal | raster cells | raster cells (non-`NA` = footprint) | one vector geometry: point, line or polygon |
| Depth | layer index, one set of depths for the grid | `[depth_min, depth_max]` per cell | top and bottom surfaces over the footprint |
| Time | none | none | optional interval |

- **Scope:** one depth interval per horizontal position (each vertical line enters the shape once and leaves once). This covers extruded shapes (cylinders, slabs, an MPA × a depth band), spheres, ellipsoids, cones, and shapes that follow the seafloor or mixed layer. Full 3D solids (overhangs, several depth bands in one column) are out of scope.
- **Surfaces:** `top` and `bottom` are each a constant (first version), a named shape with parameters (e.g. an ellipsoid buffer), a raster (e.g. relative to GEBCO bathymetry), or, as a fallback, a function of position. Every form reduces to "depth interval at these positions", so adding a shape doesn't change extraction code. `top <= bottom` everywhere in the footprint.
- **Extraction weight** per cell is the solid's volume inside that cell.
- **Conversion** to a `SpatEnvelope` evaluates top and bottom at each cell of a template grid, like `vect_to_envelope()`.
- `SpatVolume` extends to `c("SpatVoxel", "SpatEnvelope", "SpatSolid")`. Whether `SpatSolid` extends `terra::SpatVector` or something else depends on the core object chosen by the trial.

**Order of work.** Each step is its own issue and PR off `main`:

1. Document the current split: which functions take netCDF and which take `SpatVoxel` (README function overview and each function's `@param`).
2. Move the shared CF/netCDF helpers (time decoding, depth axis, bounds, dimension detection) out of `R/extract_to_points.R` into their own files.
3. Run the stars trial and write up the results with a recommendation for the core object.
4. Agree on the core data model and update this section.
5. Build the core object and lazy reader, then time-aware volume, overlap and profiles, then `SpatVoxel` as a time-step view.
6. Add `region_3d()` / `extract_to_region()` and reimplement the extraction functions as wrappers.
7. Add `SpatSolid` (one geometry, constant top and bottom, optional time) and its conversion to `SpatEnvelope`; then named shapes, raster and function surfaces.

The existing tests for both models, including the buffer and track tests, are the reference behaviour: they keep passing, and new tests check that the same synthetic data gives the same values from a netCDF file, from the core object and from the equivalent `SpatVoxel`. API changes (what `woa_load_nc()` and the voxel constructors return, time in volume and overlap) need NEWS entries and should settle before the rOpenSci and CRAN submissions (#12).

## `ocean3d` functionality

### Spatial extractions 
- **Point extraction**. `extract_to_point()` extracts the nearest value from netCDF files (ex. environmental data) to an input point based on its longitude, latitude, depth, and, when available, time.
- **Buffer extraction**. `extract_to_buffer()` summarises netCDF values within a geodesic circle × depth window × time window around each observation, with cell coverage, layer thickness and time weighting, or picks the value closest to a target.
- **Track extraction**. `extract_to_track()` samples points along ordered tracks (ex. tagging data), extracts netCDF values at them, and can summarise per track.
- **Area extraction**. `extract_to_area()` crops a `SpatVoxel` (ex. environmental data) to an input polygon and depth band. 
- **Volume extraction**. Extract values from multi-depth raster / voxel model (ex. environmental data) within a three-dimensional spatial volume defined by horizontal extent and depth range. 
  - This should be possible for both 2.5D (two-raster stack) and 3D (multi-depth raster stack). 

### Conversions between spatial representations
- **Create template voxel grid for study space**. Create the 3D voxel grid for the study space, that can be used as a template to convert other spatial representations to be compatible with it. 
- **Convert 2D polygons with depth data to 2.5D.** When separate depth ranges are available, use that to convert 2D polygons into 2.5D representations that are bound by bathymetry and depth ranges. Ex. IUCN Red List species ranges. 
- **Convert 2D rasters to 2.5D**. When depth ranges are available, use that to convert 2D raster data into 2.5D representations. Ex. working with Global Fishing Watch data on fishing effort by gear type.

### Spatial query
Three verbs cover every pairing of a 3D object (`SpatEnvelope`, `SpatVoxel`) with another 3D object or with a 2D one (a `SpatRaster` footprint, or polygons). A 2D object restricts the domain horizontally and leaves its depths alone.
- **`intersect_3d(x, y)`**: the 3D space the two share. Two envelopes intersect their depth intervals; two voxels keep the levels both occupy; an envelope is placed on a voxel's levels before the two are compared. The result is a domain, never field values.
- **`intersects_3d(x, y)`**: whether they share any 3D space, cell by cell, with terra's three answers: `TRUE`, `FALSE` (both present but disjoint, or only one present), `NA` (neither present).
- **`mask(x, mask)`**: `terra::mask()` made depth-aware. Keeps `x`'s values where the mask domain reaches that cell at that depth. A voxel masked by an envelope puts the envelope on the voxel's own levels, so the caller cannot misalign them.

`terra::intersect()` is a 2D test, so on a 3D object it is an error pointing to the verbs above. Dispatch names concrete classes only: S4 prefers a `SpatRaster` method over one on the `SpatVolume` union, so the union is a documentation device, not a dispatch target.

### Volume calculations 
- **Calculate the volume occupied by 2.5D representation**. Each rasterized range stores presence plus `depth_min/depth_max` per cell, clamped to the seafloor. Volume is the sum of `cell_area × (depth_max - depth_min)` over present cells. The overlap between two domains is the volume of `intersect_3d()`'s result, measured the same way.
- **Calculate the volume occupied by 3D representation**. Each raster layer has a height that it represents, which is defined by the height difference from the next raster layer. Ex. `t_an_depth=0`, `t_an_depth=100` would mean that the first `t_an_depth=0` layer has a height of 100m. 

### Utilities 
- **Download and prepare data.** Utilities for public APIs and open data incorporated into `ocean3d`. Includes WOA, Copernicus, IUCN Red List, Global Fishing Watch data. 
- **Downloads are cached, not manual.** Downloaded big data is cached in location accessible by package functions. Different projects read off the same copy of the data. Keeps analysis reproducible while compact. Downloaded data that is cached is documented by the package, so that the code can be handed to another user with clear information about what versions were used without having to send the big data alongside the code. 

## Conventions
- **Depth sign convention.** Depths are positive metres increasing downward, but GEBCO bathymetry is negative below sea level. The voxel constructor flips it and clamps land to 0. Document which convention any new argument uses.
- **`ncdf4` vs `terra`**. Point, buffer and track extraction read netCDF directly with `ncdf4`, reading only the cells each observation needs. Large raster operations (1000s of cells), like extracting areas and volumes, use `terra`. This split is planned to become one model; see [Data model](#data-model). Until then, say in each function's docs which input it takes.
- **Time.** Times are decoded from CF units and calendars and matched to observations in UTC. Only the netCDF model has a time dimension for now.

## Future directions

These come after the sections above are established, and are not grounded in existing project code.

### 3D species distribution modelling

Create 3D species distribution models from point observations, combining horizontal (X, Y) occurrence data with vertical (Z) depth information — extending traditional 2D SDMs by incorporating depth as an explicit dimension. 

- `create_3d_sdm(occurrences, bathymetry, env_rasters, depth_breaks)` — Build a 3D species distribution model from point observations with depth (X, Y, Z). Fits a model (e.g., MaxEnt, GLM) at each depth layer using environmental covariates extracted at that depth. Returns a multi-layer SpatRaster of predicted habitat suitability by depth.
- `predict_3d_habitat(model, env_rasters, depth_breaks, bathymetry)` — Generate 3D habitat suitability predictions from a fitted model. For each depth layer, extract environmental values and predict suitability. Mask cells where depth layer exceeds bathymetry. Returns multi-layer SpatRaster.
- **Handle spatial-autocorrelation in 3D**. Spatial autocorrelation in 2D is known challenge, tools in spatial modelling exist to handle. Need to investigate what is available / possible for spatial-autocorrelation in 3D. 
