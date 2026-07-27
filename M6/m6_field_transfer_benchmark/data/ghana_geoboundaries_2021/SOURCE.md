# Ghana administrative boundaries (2021)

These ESRI Shapefiles were generated in R from the unsimplified
`geoBoundaries` **gbOpen GHA ADM1** release:

- Boundary identifier: `GHA-ADM1-69750345`
- Administrative level: ADM1 (16 regions)
- Year represented: 2021
- Upstream source: OpenStreetMap / osm-boundaries.com
- Upstream boundary license: Creative Commons Attribution-ShareAlike 2.0
- geoBoundaries build: 12 December 2023
- Pinned upstream release commit: `9469f09`
- API metadata: <https://www.geoboundaries.org/api/current/gbOpen/GHA/ADM1/>
- Pinned archive: <https://github.com/wmgeolab/geoBoundaries/raw/9469f09/releaseData/gbOpen/GHA/ADM1/geoBoundaries-GHA-ADM1-all.zip>

`ghana_adm1_2021.*` contains the 16 source regional polygons.
`ghana_adm0_2021.*` is the national outline derived with `sf::st_union()` from
those same polygons so the exterior and regional boundaries are exactly aligned.

The local copies are unsimplified. They contain 79,310 ADM1 coordinate vertices
and 31,041 vertices after union to ADM0, compared with the 1:110m Natural Earth
outline previously used for Figure 1.

Recommended citation:

> Runfola, D. et al. (2020). geoBoundaries: A global database of political
> administrative boundaries. *PLOS ONE*, 15(4), e0231866.
> https://doi.org/10.1371/journal.pone.0231866

Refresh command (run from the repository root):

```r
adm1 <- geobounds::gb_get_adm1(
  "GHA", simplified = FALSE, release_type = "gbOpen"
)
```
