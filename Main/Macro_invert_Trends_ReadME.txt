# GRSM Macroinvertebrate Trends & Mapping

**Summary:** Cleans and aggregates GRSM macroinvertebrate data (with NCBI genus scores), links results to park hydrology, and produces time-series panels and maps, including per-stream trend estimates. In short: build stream/year (and summer) metrics, QC the sampling record, and visualize trends on both charts and maps.

## What this script does (condensed)
1. **Load inputs & paths** – Sets repo/Drive paths; reads macroinvertebrate tables and GRSM spatial layers (watersheds, streams, park border).
2. **Parse & join taxonomy** – Derives Genus, aggregates NCBI to genus level, cleans strings, and joins genus scores to specimens.
3. **QC the sampling record** – Flags duplicates and builds Year×Month timetables by sample code and by location/site for coverage checks.
4. **Aggregate metrics** – Collapses data from sample → day → site-year → stream-year (and stream-year-month) to compute richness, abundance, EPT richness/abundance, NCBI, and EPT proportion; tallies n_years_sampled.
5. **Attach geometry for mapping** – Links aggregated stream metrics to representative sampled coordinates; filters to sampled streams; creates label points.
6. **Visualize & analyze trends** – Generates time-series panels (annual and summer-only) with optional log scales and dashed LM fits; fits per-stream linear trends and maps slope symbols over GRSM streams.
