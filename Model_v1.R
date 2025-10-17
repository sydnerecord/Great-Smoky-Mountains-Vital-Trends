library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(lubridate)

#------- Load Locations & Data ----------

#---- Fish Data------

three_pass_locations <- read_csv("/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Fish/Three_Pass/Locations/GRSM_THREE_PASS.csv")
three_pass_loc_codes = three_pass_locations %>%
  dplyr::select(LAT, LON, WATERSHED, STATION_NAME, STREAMNAME) %>%
  distinct()

three_pass_data <- read_xlsx(
  "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Old/GRSM File Transfers/GRSM_Fish_3-Pass_Summary.xlsx",
  sheet = "Summary"
)
three_pass_data = three_pass_data %>%
  rename(STATION_NAME = Site) %>%
  left_join(three_pass_loc_codes, by = "STATION_NAME") %>%
  rename(LOC_NAME = STATION_NAME)
str(three_pass_data)
#write_csv(three_pass_data, "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Old/GRSM File Transfers/GRSM_Fish_3-Pass_Summary_locations.csv")

# Get response variables of interest
three_pass_summary <- three_pass_data %>%
  mutate(Year = year(Date)) %>%
  group_by(LOC_NAME, Year, WATERSHED, STREAMNAME, LAT, LON) %>%
  summarise(
    Richness = n_distinct(Species),
    ADTDens = mean(ADTDens, na.rm = TRUE),
    ADTBiom = mean(ADTBiom, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  dplyr::select(LOC_NAME, STREAMNAME, WATERSHED, LAT, LON,
                Year, Richness, ADTDens, ADTBiom)
str(three_pass_summary )

#---- Invert Data------
inverts_data <- read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export.csv')
invert_locations <- read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/Documents/Locations.csv')
invert_loc_codes <- invert_locations %>%
  mutate(
    station_code = STATION_NAME %>%
      str_remove("(I&M|I\\&M)$") %>%      # drop trailing "I&M"
      str_replace_all("[^A-Z0-9]", "")    # keep only A–Z, 0–9
  ) %>%
 dplyr::select(station_code, STATION_NAME, LOC_NAME, Watershed, StreamName, LAT, LON) %>%
  distinct()

 inverts_data = inverts_data %>%
  mutate(
    station_code2 = str_sub(Sample_Code, 1, 6)) %>% # a check
  left_join(invert_loc_codes, by = "LOC_NAME")
  
str(inverts_data)
str(invert_locations)
#save
#write_csv(inverts_data, '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export_locations.csv')


# Genus Richness
genus_rich <- inverts_data %>%
  mutate(
    Year  = year(coalesce(Start_Date, as.Date(EventDate))),
    Genus = sub("\\s.*", "", Lab_Scientific_Name)  # first token before space
  ) %>%
  group_by(LOC_NAME, Year) %>%
  summarise(Genus_richness = n_distinct(Genus, na.rm = TRUE), .groups = "drop")

#EPT
ept = read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/NCBI_EPT/EPT_site_summary.csv')

#NCBI
ncbi = read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/NCBI_EPT/NCBI_site_results.csv')

 #combine



# 1) make a common key to match EPT/NCBI to invert_loc_codes: "Location, Site" == LOC_NAME
ept_clean <- ept %>%
  mutate(LOC_NAME = paste0(Location, ", ", Site)) %>%
  dplyr::select(LOC_NAME, Year, EPT_rich)

ncbi_clean <- ncbi %>%
  mutate(LOC_NAME = paste0(Location, ", ", Site)) %>%
  dplyr::select(LOC_NAME, Year, NCBI)

# 2) combine EPT + NCBI (keep all site-years present in either table)
inv_metrics <- ept_clean %>%
  full_join(ncbi_clean, by = c("LOC_NAME", "Year"))

# 3) join genus richness + spatial columns
inverts_summary <- inv_metrics %>%
  left_join(genus_rich, by = c("LOC_NAME", "Year")) %>%
  left_join(
    invert_loc_codes %>%
      dplyr::select(LOC_NAME, STATION_NAME, Watershed, StreamName, LAT, LON),
    by = "LOC_NAME"
  ) %>%
  dplyr::select(
    LOC_NAME, STATION_NAME, Watershed, StreamName, LAT, LON,
    Year, EPT_rich, NCBI, Genus_richness
  )

inverts_summary
write_csv(inverts_summary, "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/NCBI_EPT/inverts_full_summary.csv")
# result: one row per site-year with spatial columns + EPT_rich + NCBI
inverts_summary

#-------------- Soils -----------
soil_noland  <- read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Soil_Quality/Soil_water_noland_divide_meta.csv') %>% janitor::clean_names()
soils <- read_csv("~/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Soil_Quality/Soils.csv")
soils_loc = read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Soil_Quality/soil_locations_all.csv')

soils_data = soils %>% 
  left_join(soils_loc, by = "LOC_NAME") %>%
  mutate(Year = year(Event_Date))

write_csv(soils_data, "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Soil_Quality/soil_data_and_locations.csv")


#--------- Vegetation -------
forest_locations    <- read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Forest_Health/Locations.csv')
forest_loc = forest_locations %>%
  dplyr::select(LOC_NAME, VS_WATERSHED, LAT, LON)

seedlings = read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Forest_Health/Seedlings.csv')
trees = read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Forest_Health/Trees.csv')

# Join location to data
seedlings_loc = seedlings %>% 
  left_join(forest_loc, by = "LOC_NAME") %>%
  mutate(Year = year(Event_Date))

trees_loc = trees %>% 
  left_join(forest_loc, by = "LOC_NAME") %>%
  mutate(Year = year(Event_Date))

write_csv(seedlings_loc, "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Forest_Health/Seedlings_with_locations.csv'")
write_csv(trees_loc, "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Forest_Health/trees_with_locations.csv'")

# Get abundance and Richness



# trees: richness = distinct species; abundance = count of live stems
trees_summary <- trees_loc %>%
  filter(is.na(Vigor_Desc) | Vigor_Desc != "Dead") %>%
  group_by(LOC_NAME, Year, VS_WATERSHED, LAT, LON) %>%
  summarise(
    Tree_richness  = n_distinct(SpeciesCode),
    Tree_abundance = dplyr::n(),
    .groups = "drop"
  )

# seedlings: richness = distinct species; abundance = sum of corrected stems
seedlings_summary <- seedlings_loc %>%
  group_by(LOC_NAME, Year, VS_WATERSHED, LAT, LON) %>%
  summarise(
    Seedling_richness  = n_distinct(SpeciesCode),
    Seedling_abundance = sum(Stem_Count * CORRECTION_FACTOR, na.rm = TRUE),
    .groups = "drop"
  )

# combine to one dataset per site-year with watershed, lat, lon
plant_diversity <- full_join(
  trees_summary,
  seedlings_summary,
  by = c("LOC_NAME", "Year", "VS_WATERSHED", "LAT", "LON")
) %>%
  mutate(
    Total_richness  = rowSums(dplyr::select(., dplyr::any_of(c("Tree_richness", "Seedling_richness"))), na.rm = TRUE),
    Total_abundance = rowSums(dplyr::select(., dplyr::any_of(c("Tree_abundance", "Seedling_abundance"))), na.rm = TRUE)
  ) %>%
  dplyr::select(
    LOC_NAME, VS_WATERSHED, LAT, LON, Year,
    Tree_richness, Tree_abundance,
    Seedling_richness, Seedling_abundance,
    Total_richness, Total_abundance
  ) %>%
  arrange(LOC_NAME, Year)

plant_diversity

#--------- Climate -------

temp =  rast('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Climate/PRISM_GRSM_precip_annual_mean_stack.tif')

precip = rast('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Climate/PRISM_GRSM_temp_annual_mean_stack.tif')


#---------------------------
#######################################################
#--------- Step 3 - Integrate Layers -----------------
#######################################################

# Add geometry and standarize
fish_sf <- three_pass_summary %>%
  filter(!is.na(LON) & !is.na(LAT)) %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4269)

inverts_sf <- inverts_summary %>%
  filter(!is.na(LON) & !is.na(LAT)) %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4269)

plants_sf <- plant_diversity %>%
  filter(!is.na(LON) & !is.na(LAT)) %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4269)

soils_sf <- soils_data %>%
  filter(!is.na(lon) & !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4269) #or LAT/LON depending on your column names

library(dplyr)
library(sf)
library(lubridate)

# --- 0) pick a projected CRS in meters (GRSM ~ UTM zone 17N)
crs_m <- 26917
fish_p    <- st_transform(fish_sf, crs_m)
inv_p     <- st_transform(inverts_sf, crs_m)
plant_p   <- st_transform(plants_sf, crs_m)
soils_p   <- st_transform(soils_sf, crs_m)

# --- 1) nearest neighbor indices + distances (meters)
# use fish as a baseline
idx_inv   <- st_nearest_feature(fish_p,   inv_p)
idx_plant <- st_nearest_feature(fish_p,   plant_p)
idx_soil  <- st_nearest_feature(fish_p,   soils_p)

dist_inv   <- as.numeric(st_distance(fish_p,   inv_p[idx_inv,   ], by_element = TRUE))
dist_plant <- as.numeric(st_distance(fish_p,   plant_p[idx_plant,], by_element = TRUE))
dist_soil  <- as.numeric(st_distance(fish_p,   soils_p[idx_soil, ], by_element = TRUE))

# --- 2) make a fish table (no geometry) with nearest-site IDs + distances
fish_nn <- fish_sf %>%
  st_drop_geometry() %>%
  mutate(
    nearest_invert = inverts_summary$LOC_NAME[idx_inv],
    nearest_plant  = plant_diversity$LOC_NAME[idx_plant],
    nearest_soil   = soils_data$LOC_NAME[idx_soil],
    dist_invert_m  = dist_inv,
    dist_plant_m   = dist_plant,
    dist_soil_m    = dist_soil,
    same_ws_invert = WATERSHED == inverts_summary$Watershed[idx_inv],
    same_ws_plant  = WATERSHED == plant_diversity$VS_WATERSHED[idx_plant]
  )

# --- 3) prep per-dataset minimal tables to join by nearest site + Year
inv_vars <- inverts_summary %>%
  select(LOC_NAME, Year, EPT_rich, NCBI, Genus_richness)

plant_vars <- plant_diversity %>%
  select(LOC_NAME, Year,
         Tree_richness, Tree_abundance,
         Seedling_richness, Seedling_abundance,
         Total_richness, Total_abundance)


# --- 4) join attributes onto fish by nearest site + Year
cross_taxon <- fish_nn %>%
  # inverts
  left_join(inv_vars,
            by = c("nearest_invert" = "LOC_NAME", "Year" = "Year")) %>%
  # plants
  left_join(plant_vars,
            by = c("nearest_plant" = "LOC_NAME", "Year" = "Year"),
            suffix = c("", "_plant")) %>%
  # soils
  left_join(soils_data,
            by = c("nearest_soil" = "LOC_NAME", "Year" = "Year"))

cross_taxon <- cross_taxon %>%
  rename(
    Fish_Richness = Richness,
    Fish_ADTDens  = ADTDens,
    Fish_ADTBiom  = ADTBiom,
    Invert_EPT_rich        = EPT_rich,
    Invert_NCBI            = NCBI,
    Invert_Genus_richness  = Genus_richness,
    Total_woody_richness     = Total_richness,
    Total_woody_abundance    = Total_abundance,
    Soil_pH   = SoilpH,
    Soil_CEC  = CEC,
    Soil_pOM  = pOM,
    Soil_Ca_ppm = Ca_ppm,
    Soil_Mg_ppm = Mg_ppm,
    Soil_K_ppm  = K_ppm
  ) %>%
  mutate(
    same_ws_invert = if_else(same_ws_invert, "Y", "N"),
    same_ws_plant  = if_else(same_ws_plant,  "Y", "N")
  )

# --- 5) optional filter: keep pairs within 500 m or same watershed (tune threshold as needed)
max_dist_m <- 800
cross_taxon_filtered <- cross_taxon %>%
  filter(dist_invert_m <= max_dist_m | same_ws_invert,
         dist_plant_m  <= max_dist_m | same_ws_plant)

#Check how close variables are to each other


# 1) stack points you want to cluster (use one point per record; fish is fine too)
# If you want all taxa together, bind with a type label.

# first use fish as a baseline
fish_pts   <- fish_sf %>% mutate(type = "fish")    # already EPSG:4269
invert_pts <- inverts_sf %>% mutate(type = "invert")
plant_pts  <- plants_sf %>% mutate(type = "plant")
soil_pts   <- soils_sf  %>% mutate(type = "soil")

all_pts <- rbind(fish_pts, invert_pts, plant_pts, soil_pts)

# 2) project to meters and compute within-radius groups
crs_m <- 26917
all_m <- st_transform(all_pts, crs_m)

r <- 2000  # try 2 km first (adjust after looking at coverage_tbl)
within_ix <- st_is_within_distance(all_m, all_m, dist = r)

# form cluster IDs: first point’s index in each connected component
cluster_id <- rep(NA_integer_, length(within_ix))
current <- 0L
unassigned <- seq_along(within_ix)

while(length(unassigned)){
  current <- current + 1L
  # BFS/DFS over neighbors
  to_visit <- unassigned[1]
  component <- integer(0)
  queue <- to_visit
  while(length(queue)){
    i <- queue[1]; queue <- queue[-1]
    if(i %in% component) next
    component <- c(component, i)
    nbs <- within_ix[[i]]
    queue <- c(queue, setdiff(nbs, component))
  }
  cluster_id[component] <- current
  unassigned <- setdiff(unassigned, component)
}

all_with_cluster <- all_pts %>%
  mutate(cluster_id = cluster_id)

# 3) compute cluster centroids and summaries (example stats)
cluster_summ <- all_with_cluster %>%
  st_drop_geometry() %>%
  group_by(cluster_id) %>%
  summarise(
    n_total   = n(),
    n_fish    = sum(type == "fish"),
    n_invert  = sum(type == "invert"),
    n_plant   = sum(type == "plant"),
    n_soil    = sum(type == "soil")
  ) %>%
  left_join(
    st_as_sf(all_with_cluster) %>%
      group_by(cluster_id) %>%
      summarise(geometry = st_centroid(st_union(geometry))),
    by = "cluster_id"
  ) %>%
  st_as_sf()

#inverts as baseline

rads <- c(800, 1200, 1600, 2000, 3200, 5000)

coverage_inv_tbl <- map_dfr(rads, function(r){
  df <- cross_taxon %>%
    mutate(
      has_fish  = (dist_invert_m <= r) | (same_ws_invert %in% c(TRUE, "Y")),  # reuse invert distances as fish→invert symmetric
      has_plant = (dist_plant_m  <= r) | (same_ws_plant  %in% c(TRUE, "Y")),
      has_soil  =  dist_soil_m   <= r
    )
  
  tibble(
    radius_m        = r,
    n_rows          = nrow(df),
    pct_has_fish    = mean(df$has_fish,  na.rm = TRUE),
    pct_has_plant   = mean(df$has_plant, na.rm = TRUE),
    pct_has_both_FP = mean(df$has_fish & df$has_plant, na.rm = TRUE),
    pct_has_all     = mean(df$has_fish & df$has_plant & df$has_soil, na.rm = TRUE),
    median_d_fish   = median(df$dist_invert_m, na.rm = TRUE),
    median_d_plant  = median(df$dist_plant_m,  na.rm = TRUE),
    median_d_soil   = median(df$dist_soil_m,   na.rm = TRUE)
  )
})
coverage_inv_tbl 



# use 3.2 km radius

library(sf)
library(dplyr)

# ---- 1. combine all taxon datasets into one sf ----
fish_pts   <- fish_sf   %>% mutate(taxon = "fish")
invert_pts <- inverts_sf %>% mutate(taxon = "invert")
plant_pts  <- plants_sf  %>% mutate(taxon = "plant")
soil_pts   <- soils_sf   %>% mutate(taxon = "soil")

all_pts <- bind_rows(fish_pts, invert_pts, plant_pts, soil_pts)

# ---- 2. project to meters (NAD83 / UTM zone 17N for Great Smokies) ----
all_m <- st_transform(all_pts, 26917)

# ---- 3. build clusters: all points within 3.2 km of one another ----
r_km <- 3.2
r_m  <- r_km * 1000

# adjacency list: for each point, which others are within 3.2 km
within_ix <- st_is_within_distance(all_m, all_m, dist = r_m)

# assign cluster IDs
cluster_id <- rep(NA_integer_, length(within_ix))
current <- 0L
unassigned <- seq_along(within_ix)

while (length(unassigned)) {
  current <- current + 1L
  to_visit <- unassigned[1]
  component <- integer(0)
  queue <- to_visit
  while (length(queue)) {
    i <- queue[1]; queue <- queue[-1]
    if (i %in% component) next
    component <- c(component, i)
    nbs <- within_ix[[i]]
    queue <- c(queue, setdiff(nbs, component))
  }
  cluster_id[component] <- current
  unassigned <- setdiff(unassigned, component)
}

all_with_cluster <- all_pts %>%
  mutate(cluster_id = cluster_id)

# ---- 4. summarise each cluster ----
cluster_summary <- all_with_cluster %>%
  st_drop_geometry() %>%
  group_by(cluster_id) %>%
  summarise(
    n_total   = n(),
    n_fish    = sum(taxon == "fish"),
    n_invert  = sum(taxon == "invert"),
    n_plant   = sum(taxon == "plant"),
    n_soil    = sum(taxon == "soil")
  )

# ---- 5. get centroid for each cluster ----
cluster_centroids <- all_with_cluster %>%
  group_by(cluster_id) %>%
  summarise(geometry = st_centroid(st_union(geometry))) %>%
  st_as_sf()

# ---- 6. join counts + centroids into one sf object ----
clusters_3_2km <- left_join(cluster_centroids, cluster_summary, by = "cluster_id")



# ------ combine variables within radius


# 1) Make per-taxon lookup tables: LOC_NAME -> cluster_id
fish_lookup  <- all_with_cluster %>% filter(taxon == "fish")   %>%
  st_drop_geometry() %>% distinct(LOC_NAME, cluster_id)

inv_lookup   <- all_with_cluster %>% filter(taxon == "invert") %>%
  st_drop_geometry() %>% distinct(LOC_NAME, cluster_id)

plant_lookup <- all_with_cluster %>% filter(taxon == "plant")  %>%
  st_drop_geometry() %>% distinct(LOC_NAME, cluster_id)

soil_lookup  <- all_with_cluster %>% filter(taxon == "soil")   %>%
  st_drop_geometry() %>% distinct(LOC_NAME, cluster_id)

# 2) Prepare per-taxon tables with clear names + cluster_id
# FISH (already per site-year)
fish_dat <- three_pass_summary %>%
  left_join(fish_lookup, by = "LOC_NAME") %>%
  rename(
    Fish_Richness = Richness,
    Fish_ADTDens  = ADTDens,
    Fish_ADTBiom  = ADTBiom
  ) %>%
  select(cluster_id, Year, Fish_Richness, Fish_ADTDens, Fish_ADTBiom)

# INVERTS (per site-year)
inv_dat <- inverts_summary %>%
  left_join(inv_lookup, by = "LOC_NAME") %>%
  rename(
    Invert_EPT_rich       = EPT_rich,
    Invert_NCBI           = NCBI,
    Invert_Genus_richness = Genus_richness
  ) %>%
  select(cluster_id, Year, Invert_EPT_rich, Invert_NCBI, Invert_Genus_richness)

# PLANTS (per site-year)
plant_dat <- plant_diversity %>%
  left_join(plant_lookup, by = "LOC_NAME") %>%
  rename(
    Plant_Tree_richness      = Tree_richness,
    Plant_Tree_abundance     = Tree_abundance,
    Plant_Seedling_richness  = Seedling_richness,
    Plant_Seedling_abundance = Seedling_abundance
    # If you also want totals, add:
    # ,Plant_Total_richness     = Total_richness
    # ,Plant_Total_abundance    = Total_abundance
  ) %>%
  select(cluster_id, Year,
         Plant_Tree_richness, Plant_Tree_abundance,
         Plant_Seedling_richness, Plant_Seedling_abundance
         # ,Plant_Total_richness, Plant_Total_abundance
  )

# SOILS (may have multiple samples per site-year; reduce to site-year first)
soils_site_year <- soils_data %>%
  mutate(Year = year(Event_Date)) %>%
  group_by(LOC_NAME, Year) %>%
  summarise(
    Soil_pH   = mean(SoilpH, na.rm = TRUE),
    Soil_CEC  = mean(CEC,    na.rm = TRUE),
    Soil_pOM  = mean(pOM,    na.rm = TRUE),
    Soil_Ca_ppm = mean(Ca_ppm, na.rm = TRUE),
    Soil_Mg_ppm = mean(Mg_ppm, na.rm = TRUE),
    Soil_K_ppm  = mean(K_ppm,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(soil_lookup, by = "LOC_NAME") %>%
  select(cluster_id, Year, Soil_pH, Soil_CEC, Soil_pOM, Soil_Ca_ppm, Soil_Mg_ppm, Soil_K_ppm)

# 3) Aggregate to cluster-year
# Choices:
#   - Richness indices: MEAN across sites in cluster-year (could also use median)
#   - Densities/biomass: MEAN
#   - Abundances (counts): SUM (or mean if they’re densities/standardized)
cluster_fish <- fish_dat %>%
  group_by(cluster_id, Year) %>%
  summarise(
    Fish_Richness = mean(Fish_Richness, na.rm = TRUE),
    Fish_ADTDens  = mean(Fish_ADTDens,  na.rm = TRUE),
    Fish_ADTBiom  = mean(Fish_ADTBiom,  na.rm = TRUE),
    .groups = "drop"
  )

cluster_inv <- inv_dat %>%
  group_by(cluster_id, Year) %>%
  summarise(
    Invert_EPT_rich       = mean(Invert_EPT_rich,       na.rm = TRUE),
    Invert_NCBI           = mean(Invert_NCBI,           na.rm = TRUE),
    Invert_Genus_richness = mean(Invert_Genus_richness, na.rm = TRUE),
    .groups = "drop"
  )

cluster_plant <- plant_dat %>%
  group_by(cluster_id, Year) %>%
  summarise(
    Plant_Tree_richness      = mean(Plant_Tree_richness,      na.rm = TRUE),
    Plant_Seedling_richness  = mean(Plant_Seedling_richness,  na.rm = TRUE),
    Plant_Tree_abundance     = sum(Plant_Tree_abundance,      na.rm = TRUE),
    Plant_Seedling_abundance = sum(Plant_Seedling_abundance,  na.rm = TRUE),
    # If you kept totals:
    # Plant_Total_richness     = mean(Plant_Total_richness,     na.rm = TRUE),
    # Plant_Total_abundance    = sum(Plant_Total_abundance,     na.rm = TRUE),
    .groups = "drop"
  )

cluster_soil <- soils_site_year %>%
  group_by(cluster_id, Year) %>%
  summarise(
    Soil_pH   = mean(Soil_pH,   na.rm = TRUE),
    Soil_CEC  = mean(Soil_CEC,  na.rm = TRUE),
    Soil_pOM  = mean(Soil_pOM,  na.rm = TRUE),
    Soil_Ca_ppm = mean(Soil_Ca_ppm, na.rm = TRUE),
    Soil_Mg_ppm = mean(Soil_Mg_ppm, na.rm = TRUE),
    Soil_K_ppm  = mean(Soil_K_ppm,  na.rm = TRUE),
    .groups = "drop"
  )

# 4) Combine all cluster-year tables and attach climate
cluster_year <- cluster_fish %>%
  full_join(cluster_inv,   by = c("cluster_id","Year")) %>%
  full_join(cluster_plant, by = c("cluster_id","Year")) %>%
  full_join(cluster_soil,  by = c("cluster_id","Year")) 


  # add climate at the centroid (cluster_id, Year)
  left_join(clusters_climate, by = c("cluster_id","Year")) %>%
  arrange(cluster_id, Year)

# 5) Bring centroid geometry if you want an sf table (geometry repeats by year)
# Keep geometry by joining with the sf object on the left 
  cluster_year_sf <- cluster_centroids_sf %>%
    right_join(cluster_year, by = "cluster_id")


# Add Climate
  
library(dplyr)
library(sf)
library(terra)
library(stringr)

# 0) Identify stacks (yours are named opposite)
prism_ppt   <- temp     # ppt_* layers
prism_tmean <- precip   # tmean_* layers

# 1) Make sure cluster centroids are in 4269 to match PRISM
clusters_ll <- st_transform(cluster_year_sf, 4269)

# 2) Parse the years available in each stack from layer names
ppt_years   <- as.integer(str_extract(names(prism_ppt),   "\\d{4}"))
tmean_years <- as.integer(str_extract(names(prism_tmean), "\\d{4}"))

# 3) Prepare a terra SpatVector of cluster points (same order as clusters_ll)
clusters_v <- terra::vect(clusters_ll)

# 4) Loop over years: extract the single matching ppt & tmean layer for all rows of that year
yrs <- sort(unique(clusters_ll$Year))
out_list <- list()

for (y in yrs) {
  i_ppt   <- which(ppt_years   == y)
  i_tmean <- which(tmean_years == y)
  if (!length(i_ppt) || !length(i_tmean)) {
    # No matching climate layer for this year → fill with NA for those rows
    out_list[[length(out_list)+1]] <- clusters_ll %>%
      st_drop_geometry() %>%
      filter(Year == y) %>%
      transmute(cluster_id, Year, ppt_mm = NA_real_, tmean_C = NA_real_)
    next
  }
  
  # rows for this year (indices relative to clusters_ll)
  idx_rows <- which(clusters_ll$Year == y)
  
  # extract at those rows for the year's layers
  vals_ppt   <- terra::extract(prism_ppt[[i_ppt[1]]],   clusters_v[idx_rows, ])
  vals_tmean <- terra::extract(prism_tmean[[i_tmean[1]]], clusters_v[idx_rows, ])
  
  # build a small table for this year
  out_list[[length(out_list)+1]] <- clusters_ll %>%
    st_drop_geometry() %>%
    slice(idx_rows) %>%
    transmute(
      cluster_id,
      Year = y,
      ppt_mm  = as.numeric(vals_ppt[[2]]),
      tmean_C = as.numeric(vals_tmean[[2]])
    )
}

clusters_climate <- dplyr::bind_rows(out_list)

# 5) Attach climate back to your cluster-year sf (keeps geometry)
cluster_year_sf <- clusters_ll %>%
  left_join(clusters_climate, by = c("cluster_id","Year"))

# (Optional) clearer names
cluster_year_sf <- cluster_year_sf %>%
  rename(Precip_mm = ppt_mm, Temp_C = tmean_C)

# (Optional) also store LAT/LON as numeric columns for export
cluster_year_sf <- cluster_year_sf %>%
  mutate(
    LON = st_coordinates(geometry)[,1],
    LAT = st_coordinates(geometry)[,2]
  )

# Plain dataframe for modeling / CSV:
cluster_year <- cluster_year_sf %>% st_drop_geometry()

# Example save:
# Diagnostics on Cluster



# 0) sanity: how many points got clustered and did any miss a cluster_id?
n_all_pts <- nrow(all_pts)
n_all_clustered <- nrow(all_with_cluster)
n_na_cluster <- sum(is.na(all_with_cluster$cluster_id))

cluster_basic <- tibble(
  n_all_pts_binded   = n_all_pts,
  n_all_pts_clustered= n_all_clustered,
  n_with_na_cluster  = n_na_cluster,
  n_clusters         = dplyr::n_distinct(all_with_cluster$cluster_id, na.rm = TRUE)
)


# some plots

#---------------------------------------------
# 0) Choose variables to track by theme
#---------------------------------------------
var_groups <- list(
  fish   = c("Fish_Richness","Fish_ADTDens","Fish_ADTBiom"),
  invert = c("Invert_EPT_rich","Invert_NCBI","Invert_Genus_richness"),
  plant  = c("Plant_Tree_richness","Plant_Seedling_richness",
             "Plant_Tree_abundance","Plant_Seedling_abundance"),
  soil   = c("Soil_pH","Soil_CEC","Soil_pOM","Soil_Ca_ppm","Soil_Mg_ppm","Soil_K_ppm"),
  clim   = c("Precip_mm","Temp_C")
)

# small helper to check presence and drop all-NA variables
vars_present <- function(df, vars){
  keep <- vars[vars %in% names(df)]
  if (length(keep) == 0) stop("None of the requested variables are present.")
  # drop vars that are entirely NA
  keep[vapply(keep, function(v) any(!is.na(df[[v]])), logical(1))]
}

#---------------------------------------------
# 1) Park-wide yearly trend: median ± IQR
#---------------------------------------------
plot_yearly_trend <- function(df, vars, title = NULL){
  keep <- vars_present(df, vars)
  df_long <- df %>%
    dplyr::select(cluster_id, Year, dplyr::all_of(keep)) %>%
    pivot_longer(cols = all_of(keep), names_to = "variable", values_to = "value")
  
  yearly_summary <- df_long %>%
    group_by(Year, variable) %>%
    summarise(
      n_clusters = n_distinct(cluster_id[!is.na(value)]),
      med = median(value, na.rm = TRUE),
      p25 = quantile(value, 0.25, na.rm = TRUE),
      p75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(yearly_summary, aes(Year, med)) +
    geom_ribbon(aes(ymin = p25, ymax = p75), alpha = 0.25) +
    geom_line(linewidth = 1) +
    facet_wrap(~ variable, scales = "free_y", ncol = 3) +
    labs(x = NULL, y = NULL, title = title, subtitle = "Median (line) with IQR (ribbon)")
}

#---------------------------------------------
# 2) Spaghetti by cluster + bold park-wide median
#   (sample clusters if you want fewer lines)
#---------------------------------------------
plot_spaghetti <- function(df, vars, title = NULL, sample_n = NULL, seed = 1){
  keep <- vars_present(df, vars)
  
  df_use <- df %>%
    dplyr::select(cluster_id, Year, dplyr::all_of(keep))
  
  if (!is.null(sample_n)) {
    set.seed(seed)
    keep_ids <- df_use %>% distinct(cluster_id) %>% slice_sample(n = sample_n) %>% pull(cluster_id)
    df_use <- df_use %>% filter(cluster_id %in% keep_ids)
  }
  
  df_long <- df_use %>%
    tidyr::pivot_longer(cols = dplyr::all_of(keep), names_to = "variable", values_to = "value")
  
  med_year <- df_long %>%
    dplyr::group_by(Year, variable) %>%
    dplyr::summarise(med = median(value, na.rm = TRUE), .groups = "drop")
  
  ggplot(df_long, aes(Year, value, group = cluster_id)) +
    geom_line(alpha = 0.25) +
    theme_plot() +
    # key fix: don't inherit group=cluster_id; group by variable for the median
    geom_line(
      data = med_year,
      aes(x = Year, y = med, group = variable),
      linewidth = 1.2,
      inherit.aes = FALSE
    ) +
    scale_y_log10() +
    facet_wrap(~ variable, scales = "free_y", ncol = 3) +
    labs(x = NULL, y = NULL, title = title, subtitle = "Thin = clusters; Thick = park-wide median")
}

#---------------------------------------------
# 3) Z-score heatmap (by variable): cluster × year
#   Good for spotting synchronous years or clusters.
#---------------------------------------------
plot_z_heatmap <- function(df, vars, title = NULL){
  keep <- vars_present(df, vars)
  
  df_long <- df %>%
    dplyr::select(cluster_id, Year, dplyr::all_of(keep)) %>%
    pivot_longer(cols = all_of(keep), names_to = "variable", values_to = "value")
  
  # z-score within each variable across all cluster-years
  df_z <- df_long %>%
    group_by(variable) %>%
    mutate(z = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>%
    ungroup()
  
  # order clusters by mean z in the most recent year (or choose a reference year)
  ref_year <- max(df_z$Year, na.rm = TRUE)
  clust_order <- df_z %>%
    filter(Year == ref_year) %>%
    group_by(cluster_id) %>%
    summarise(mz = mean(z, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mz)) %>%
    pull(cluster_id)
  
  df_z$cluster_id <- factor(df_z$cluster_id, levels = clust_order)
  
  ggplot(df_z, aes(Year, cluster_id, fill = z)) +
    geom_tile() +
    facet_wrap(~ variable, scales = "free", ncol = 3) +
    labs(x = NULL, y = "cluster_id (ordered at ref year)", title = title, subtitle = paste("Z-score per variable; ref year =", ref_year))
}

#---------------------------------------------
# 4) Examples — call what you need
#---------------------------------------------
# Park-wide medians ± IQR
plot_yearly_trend(cluster_year, var_groups$fish,  "Fish indices — yearly trend")
plot_yearly_trend(cluster_year, var_groups$invert,"Inverts — yearly trend")
plot_yearly_trend(cluster_year, var_groups$plant, "Plants — yearly trend")
plot_yearly_trend(cluster_year, var_groups$soil,  "Soils — yearly trend")
plot_yearly_trend(cluster_year, var_groups$clim,  "Climate — yearly trend")

# Spaghetti (sample 40 of your 271 clusters to keep it readable)
plot_spaghetti(cluster_year, var_groups$fish,  "Fish — spaghetti",  sample_n = 40)
plot_spaghetti(cluster_year, var_groups$invert,"Inverts — spaghetti",sample_n = 40)
plot_spaghetti(cluster_year, var_groups$plant,  "Plants — spaghetti",  sample_n = 40)
plot_spaghetti(cluster_year, var_groups$soil,   "Soils — spaghetti",   sample_n = 40)
plot_spaghetti(cluster_year, var_groups$clim,   "Climate — spaghetti", sample_n = 40)

# Z-score heatmap
plot_z_heatmap(cluster_year, var_groups$fish,  "Fish — z heatmap")
plot_z_heatmap(cluster_year, var_groups$plant, "Plants — z heatmap")
plot_z_heatmap(cluster_year, var_groups$invert,  "Invert— z heatmap")

#---------------------------------------------
# 5) Optional: smooth medians with GAM (lightweight)
#---------------------------------------------
plot_yearly_trend_smooth <- function(df, vars, title = NULL){
  keep <- vars_present(df, vars)
  df_long <- df %>%
    dplyr::select(cluster_id, Year, dplyr::all_of(keep)) %>%
    pivot_longer(cols = all_of(keep), names_to = "variable", values_to = "value")
  
  yearly_summary <- df_long %>%
    group_by(Year, variable) %>%
    summarise(
      med = median(value, na.rm = TRUE),
      p25 = quantile(value, 0.25, na.rm = TRUE),
      p75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(yearly_summary, aes(Year, med)) +
    geom_ribbon(aes(ymin = p25, ymax = p75), alpha = 0.2) +
    geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, k = 6)) +
    facet_wrap(~ variable, scales = "free_y", ncol = 3) +
    labs(x = NULL, y = NULL, title = title, subtitle = "GAM-smoothed median with IQR band")
}
# Example:
# plot_yearly_trend_smooth(cluster_year, c(var_groups$fish, var_groups$clim), "Fish + Climate (smoothed)")




# Correlations bdtween varfiables

#--- choose numeric variables only
vars_for_corr <- c(
  var_groups$fish,
  var_groups$invert,
  var_groups$plant,
  var_groups$soil,
  var_groups$clim
)

df_corr <- cluster_year %>%
  select(all_of(vars_for_corr)) %>%
  drop_na()

# Compute correlation matrix
corr_mat <- cor(df_corr, use = "pairwise.complete.obs", method = "spearman")

# Turn to long format for plotting
corr_long <- as.data.frame(as.table(corr_mat)) %>%
  rename(var1 = Var1, var2 = Var2, rho = Freq)

ggplot(corr_long, aes(var1, var2, fill = rho)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1,1)) +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Spearman correlations among cluster-year variables",
       fill = "ρ")

# 1) per-taxon site counts before clustering (unique LOC_NAME with coords)
count_sites_pre <- bind_rows(
  fish_sf   %>% mutate(taxon = "fish"),
  inverts_sf %>% mutate(taxon = "invert"),
  plants_sf  %>% mutate(taxon = "plant"),
  soils_sf   %>% mutate(taxon = "soil")
) %>%
  st_drop_geometry() %>%
  group_by(taxon) %>%
  summarise(
    n_rows          = n(),                        # site-year rows in the sf object
    n_unique_sites  = n_distinct(LOC_NAME),       # unique sites
    .groups = "drop"
  )

# 2) per-taxon “cluster lookup” coverage (did every LOC_NAME get a cluster_id?)
fish_lookup  <- all_with_cluster %>% filter(taxon == "fish")   %>% st_drop_geometry() %>% distinct(LOC_NAME, cluster_id)
inv_lookup   <- all_with_cluster %>% filter(taxon == "invert") %>% st_drop_geometry() %>% distinct(LOC_NAME, cluster_id)
plant_lookup <- all_with_cluster %>% filter(taxon == "plant")  %>% st_drop_geometry() %>% distinct(LOC_NAME, cluster_id)
soil_lookup  <- all_with_cluster %>% filter(taxon == "soil")   %>% st_drop_geometry() %>% distinct(LOC_NAME, cluster_id)

# unmatched LOC_NAMEs (present in data tables but missing from lookups)
fish_unmatched  <- three_pass_summary %>% distinct(LOC_NAME) %>% anti_join(fish_lookup,  by = "LOC_NAME")
inv_unmatched   <- inverts_summary   %>% distinct(LOC_NAME) %>% anti_join(inv_lookup,   by = "LOC_NAME")
plant_unmatched <- plant_diversity   %>% distinct(LOC_NAME) %>% anti_join(plant_lookup, by = "LOC_NAME")
soil_unmatched  <- soils_data        %>% distinct(LOC_NAME) %>% anti_join(soil_lookup,  by = "LOC_NAME")

lost_sites_tbl <- tibble(
  taxon          = c("fish","invert","plant","soil"),
  n_sites_pre    = c(
    n_distinct(three_pass_summary$LOC_NAME),
    n_distinct(inverts_summary$LOC_NAME),
    n_distinct(plant_diversity$LOC_NAME),
    n_distinct(soils_data$LOC_NAME)
  ),
  n_sites_with_cluster = c(
    n_distinct(fish_lookup$LOC_NAME),
    n_distinct(inv_lookup$LOC_NAME),
    n_distinct(plant_lookup$LOC_NAME),
    n_distinct(soil_lookup$LOC_NAME)
  ),
  n_sites_lost   = c(
    nrow(fish_unmatched),
    nrow(inv_unmatched),
    nrow(plant_unmatched),
    nrow(soil_unmatched)
  )
)

# 3) list the actual lost site IDs (if any)
lost_sites_list <- list(
  fish_lost_sites  = fish_unmatched$LOC_NAME,
  invert_lost_sites= inv_unmatched$LOC_NAME,
  plant_lost_sites = plant_unmatched$LOC_NAME,
  soil_lost_sites  = soil_unmatched$LOC_NAME
)

# 4) cluster diagnostics: size distribution and taxon composition
cluster_size_dist <- all_with_cluster %>%
  st_drop_geometry() %>%
  count(cluster_id, name = "points_in_cluster") %>%
  arrange(desc(points_in_cluster))

cluster_taxa_mix <- all_with_cluster %>%
  st_drop_geometry() %>%
  count(cluster_id, taxon) %>%
  tidyr::pivot_wider(names_from = taxon, values_from = n, values_fill = 0) %>%
  arrange(cluster_id)

# 5) coverage by cluster-year after aggregation (did any cluster-years drop out?)
cluster_year_counts <- cluster_year %>%
  group_by(cluster_id) %>%
  summarise(n_years = n_distinct(Year), .groups = "drop") %>%
  arrange(desc(n_years))

# 6) summary prints
cluster_basic
count_sites_pre
lost_sites_tbl
# If any n_sites_lost > 0, inspect:
lost_sites_list
cluster_size_dist %>% head(20)
cluster_taxa_mix %>% head(20)
cluster_year_counts %>% head(20)


