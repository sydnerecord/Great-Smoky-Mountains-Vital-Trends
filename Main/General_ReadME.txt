# Great-Smoky-Mountains-Vital-Trends
Code to analyze vital trends data for Great Smoky Mountains National Park

Description of Code

1. Spatial_Objects.R - This script builds USGS hydrology and climate layers from PRISM. It reads the park boundary, merges HUC10 watersheds (Regions 03 & 06), filters those within GRSM, and clips NHDPlus HR flowlines. Creates labeled maps for surveyed streams. Downloads 1980–2023 PRISM precipitation and temperature, masks to GRSM, computes annual means, and saves GeoTIFF stacks with quick preview plots

2. Macro_invert_trends.R – This code analyzes sampling effort (month of sampling) and generates plots of diversity indices for freshwater macroinvertebrates - abundance, richness, EPT richness and NCBI. Both scatterplot timelines and spatially explicit plots of trend slopes are generated.

3. Fish_Trends.R – This code uses three_pass fish data to analyze overall trends in fish metrics across streams - adult density, adult mass, etc. It also generates loops for species specific trends.

4. GRSM_fish_plotting.R. This builds a plotting function to use mixed model (with stream and/or species as a random effect) with Guassian or negative binomial fit. 