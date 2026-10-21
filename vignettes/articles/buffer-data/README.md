# Real-data example for buffer extraction

These files accompany `extracting-and-weighting-buffer-data.Rmd`.
They are website article resources, not installed ocean3d example data.

## Source

- Product: Mediterranean Sea Physics Reanalysis, `MEDSEA_MULTIYEAR_PHY_006_004`.
- Dataset: `cmems_mod_med_phy-temp_my_4.2km_P1D-m_202511`.
- Variable: `thetao`, sea water potential temperature, in degrees Celsius.
- Producer: Centro Euro-Mediterraneo sui Cambiamenti Climatici (CMCC), Italy.
- System recorded in the source: MFS E3R1I.
- Original download timestamp: `2026-10-19T11:44:33.421Z`.
- Original filename: `cmems_mod_med_phy-temp_my_4.2km_P1D-m_1792410273421.nc`.
- Product DOI: <https://doi.org/10.48670/mds-00375>.
- Product page: <https://data.marine.copernicus.eu/product/MEDSEA_MULTIYEAR_PHY_006_004/description>.
- Licence: <https://marine.copernicus.eu/user-corner/service-commitments-and-licence>.
- Citation guidance: <https://help.marine.copernicus.eu/en/articles/4444611-citing-copernicus-marine-products-and-services>.

E.U. Copernicus Marine Service Information; https://doi.org/10.48670/mds-00375.
The source data retain the Copernicus Marine licence rather than acquiring
the package's software licence. The original netCDF's scientific references
and source identifiers are preserved in the subset metadata. The file includes
the interim-system reference DOI `10.25423/CMCC/MEDSEA_MULTIYEAR_PHY_006_004_E3R1I`.

## Files and preparation

- `temperature-subset.nc`: small rectangular subset of the original values,
  including neighbouring coordinate centres around the requested support.
  Temperature values are unchanged. No bounds variables have been added.
  Original-file and subset provenance are recorded as global attributes.
- `temperature-geometry.rds`: clipped display polygons, individual sample
  values and weights, and layer/time support for the one illustrated buffer.
  Polygon simplification (5 m tolerance) affects rendering only.
- `figures.R` and `explorer.js`: article-only drawing helpers. Their names
  are not exported package functions. Plotly is provided by R's htmlwidget
  dependencies rather than a remote CDN.

From the package root, load the current source and regenerate with:

```r
pkgload::load_all()
source("data-raw/buffer-article.R")
prepare_buffer_article("path/to/the/original-downloaded-temperature.nc")
```

The script does not download data or edit the source file. It independently
reconstructs figure support from the coordinates and compares the resulting
mean with the public extraction function. It also compares the subset with
the original, checks that all 450 selected values are valid, and records the
subset checksum in the RDS for verification when the article is rendered.

## Example

Observation: longitude 0 degrees, latitude 38 degrees, depth 60 m,
`2026-09-20 00:00:00 UTC`. Buffer: radius 10 km, depth 50-70 m,
time `2026-09-19` to `2026-09-21 UTC`. The selected support contains 30
horizontal cells, five depth layers and three timestamps (450 values).

Weighted mean: approximately 16.05743 degrees C.
Equal mean: approximately 15.98204 degrees C.

These values illustrate support weighting; they are not a scientific estimate
of an animal's experienced temperature. Inferred support does not establish
the actual physical averaging intervals of the source data. The source
metadata's daily-mean label is separate from that inference.
