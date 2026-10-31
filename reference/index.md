# Package index

## 3D classes

Raster-based 3D objects. A `SpatVoxel` holds values by depth level; a
`SpatEnvelope` holds a depth range per cell.

- [`SpatVoxel-class`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVoxel-class.md)
  : 3D voxel model: one layer per standard depth level
- [`SpatEnvelope-class`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatEnvelope-class.md)
  : 2.5D min-max envelope: exactly depth_min, depth_max
- [`SpatVolume-class`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/SpatVolume-class.md)
  : Union of the package's 3D domain representations
- [`depths()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/depths.md)
  : Depths of a voxel's layers

## Build

Turn inputs into 3D objects, and convert between them.

- [`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md)
  : Create a voxel object
- [`as_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_envelope.md)
  : Coerce a 2D footprint to a 2.5D min-max envelope
- [`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md)
  : Convert SpatVector or sf to SpatEnvelope
- [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md)
  : Convert Envelope 2.5D -\> Voxel 3D
- [`voxel_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_to_envelope.md)
  : Collapse Voxel 3D -\> Envelope 2.5D
- [`profile_flat()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_profiles.md)
  [`profile_equal()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/voxel_profiles.md)
  : Vertical profiles for voxel building
- [`occupied()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/occupied.md)
  : Reduce a variable voxel to a presence voxel
- [`create_study_raster()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/create_study_raster.md)
  : Create a study area raster grid

## Intersect

Find where 3D objects meet.

- [`intersect_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersect_3d.md)
  : The 3D space two objects share
- [`intersects_3d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/intersects_3d.md)
  : Do two objects share any 3D space?
- [`mask(`*`<SpatVoxel>`*`,`*`<SpatEnvelope>`*`)`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md)
  [`mask(`*`<SpatVoxel>`*`,`*`<SpatVoxel>`*`)`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md)
  [`mask(`*`<SpatEnvelope>`*`,`*`<SpatEnvelope>`*`)`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md)
  [`mask(`*`<SpatEnvelope>`*`,`*`<SpatVoxel>`*`)`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md)
  [`mask(`*`<SpatRaster>`*`,`*`<SpatEnvelope>`*`)`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md)
  [`mask(`*`<SpatRaster>`*`,`*`<SpatVoxel>`*`)`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/mask-3d.md)
  : Keep a raster's values inside a 3D domain

## Summarise

Measure 3D objects and their overlap.

- [`volume()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/volume.md)
  : Total 3D volume of a rasterized domain
- [`calc_volume_overlap()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/calc_volume_overlap.md)
  : Per-cell 3D volume overlap between two rasterized domains

## Extraction

Extract environmental values from netCDF files at points, in buffers,
along tracks, or from a `SpatVoxel` within an area. These run build,
intersect and summarise in one call.

- [`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md)
  : Extract values from netCDFs to one or more point observations
- [`extract2d()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract2d.md)
  : Extract values from a 2D netCDF variable
- [`extract3d_all()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract3d_all.md)
  : Extract nearest, surface and bottom values from a 3D netCDF variable
- [`extract3d_bottom()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract3d_bottom.md)
  : Extract the bottom available layer from a 3D netCDF variable
- [`extract3d_nearest()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract3d_nearest.md)
  : Extract the nearest valid depth layer from a 3D netCDF variable
- [`extract3d_surface()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract3d_surface.md)
  : Extract the surface layer from a 3D netCDF variable
- [`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md)
  : Summarise environmental values around observations
- [`extract_to_track()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_track.md)
  : Extract environmental values along tracks and lines
- [`extract_to_area()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_area.md)
  : Extract a 3D raster to an area and a depth band

## Data sources

Download, cache and load public datasets.

- [`load_gebco_bathymetry()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/load_gebco_bathymetry.md)
  : Load GEBCO bathymetry raster
- [`woa_download()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_download.md)
  : Download a WOA NetCDF file (with caching)
- [`woa_load_nc()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_load_nc.md)
  : Load a WOA NetCDF file
- [`woa_cache_dir()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_dir.md)
  : WOA cache directory
- [`woa_cache_clear()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_clear.md)
  : Clear the WOA cache
- [`gfw_effort_to_raster()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/gfw_effort_to_raster.md)
  : Rasterise a GFW effort tibble onto a target grid
- [`fetch_species_assessments()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/fetch_species_assessments.md)
  : Fetch species assessment data from IUCN Red List API
- [`fill_missing_depths()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/fill_missing_depths.md)
  : Fill missing depth values
