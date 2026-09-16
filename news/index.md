# Changelog

## ocean3d 0.1.1.9007

The package is renamed from `sharkabc3d` to `ocean3d`, to reflect that
the 3D workflows apply to any marine taxon or spatial layer, not just
sharks and rays.

### Breaking changes

- [`library(sharkabc3d)`](https://github.com/Marine-Biodiversity-Conservation-Lab/sharkabc3d)
  becomes
  [`library(ocean3d)`](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d).
  No function names, arguments or class names changed — only the package
  name.
- The WOA download cache moves with the package name, from
  `tools::R_user_dir("sharkabc3d", "cache")/woa` to
  `tools::R_user_dir("ocean3d", "cache")/woa`. Previously downloaded WOA
  files are not found at the new location; either re-download them or
  move the old `woa` directory across.
  [`woa_cache_dir()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/woa_cache_dir.md)
  reports the new path.
- Package `Title` and `Description` are broadened from “Shark and Ray
  Abiotic Covariates in 3 Dimensions” to “Three-Dimensional Marine
  Spatial Analysis”.
- The repository and pkgdown site move to
  `Marine-Biodiversity-Conservation-Lab/ocean3d`; GitHub redirects the
  old URLs, but `devtools::install_github()` calls should be updated.
