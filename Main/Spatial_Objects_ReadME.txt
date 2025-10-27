ReadMe for 1.Spatial_Objects

1. Define paths and boundaries – Loads the Great Smoky Mountains National Park (GRSM) boundary shapefile and unifies it into a single polygon for spatial clipping.

2. Assemble watershed data – Reads HUC-10 polygons from USGS Regions 03 and 06, merges them, and clips to the GRSM bounding box to retain only park-relevant watersheds.

3. Filter and export – Removes large outlying basins, keeping only in-park watersheds, and saves them as a GeoPackage.

4. Download and process streams – Retrieves NHDPlus HR flowlines for hydrologic units 0307 and 0601, merges and reprojects them, then clips to GRSM watersheds to isolate park streams.

5. Map sampled streams – Extracts named streams sampled for macroinvertebrate and fish surveys, creates midpoint labels, and plots maps showing sampled reaches within the park.

6. Download PRISM climate data – Automates acquisition of monthly 800 m precipitation and mean-temperature rasters (1980–2023) directly from NACSE.

7. Clip and aggregate climate rasters – Crops PRISM stacks to the park extent, masks to boundaries, and averages monthly data to produce annual precipitation and temperature layers.

8. Export and visualize – Writes annual mean rasters to GeoTIFF and plots representative 1980 layers to verify coverage and spatial patterns.