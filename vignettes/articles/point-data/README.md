# Point extraction tutorial data

These are website article resources, not installed example files or new
exported ocean3d functions. No remote services are needed to build the article.

## Source and attribution

- Product: Mediterranean Sea Physics Reanalysis, `MEDSEA_MULTIYEAR_PHY_006_004`.
- Dataset: `cmems_mod_med_phy-temp_my_4.2km_P1D-m_202511`.
- Producer: Centro Euro-Mediterraneo sui Cambiamenti Climatici (CMCC), Italy.
- Original: `cmems_mod_med_phy-temp_my_4.2km_P1D-m_1792410273421.nc`.
- Download date in source metadata: `2026-10-19T11:44:33.421Z`.
- Product DOI: <https://doi.org/10.48670/mds-00375>.
- Product: <https://data.marine.copernicus.eu/product/MEDSEA_MULTIYEAR_PHY_006_004/description>.
- Licence: <https://marine.copernicus.eu/user-corner/service-commitments-and-licence>.

E.U. Copernicus Marine Service Information; https://doi.org/10.48670/mds-00375.
Source temperatures retain that data licence. The source identifiers and
scientific references are preserved in each netCDF's global attributes.

## Files

- `temperature-points.nc` retains the original `thetao` and `bottomT` values
  on a small grid, all 50 original downloaded depth levels, and 19-21 September
  2026. It also contains an explicitly artificial, static `cell_id` field
  to demonstrate variable-specific date requirements. This field is an integer
  grid identifier, not an environmental observation or a Copernicus product.
- `temperature-2026-09-19.nc` and `temperature-2026-09-21.nc` contain `thetao`
  at a single timestamp each. They demonstrate selection between files,
  including when neither file has multiple time records.
- `point-display.rds` stores plot geometry and actual results and messages
  from the public R function for a finite set of menu choices. The widgets
  do not implement general selection algorithms in JavaScript. The example's
  file checksum is checked when rendering the article.
- `figures.R` and `explorer.js` draw the static and interactive figures.
  Dependencies are supplied by Plotly's local R htmlwidget assets, not a CDN.

The fictional fishing locations do not represent real fishing activity.
The cropped grid boundary does not represent a coastline or the full product
coverage. Depth boundaries used to display 3D prisms are inferred midpoints;
they do not change the point selection or perform interpolation.

## Regenerate

From the package root:

```r
pkgload::load_all()
source("data-raw/point-article.R")
prepare_point_article("path/to/the/original-downloaded-temperature.nc")
```

The original file is opened read-only. The preparation verifies that the
nearest, shallowest and deepest-valid values for the central observation
match those from the larger original file. Browser interactions change
the display of precomputed cases, not the source data or the extraction.
