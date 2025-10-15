############################################################
## GRSM MAPS — ORGANIZED, LIGHTLY ANNOTATED (minimal edits)
## Order:
##   1) Load packages & data
##   2) Watersheds (outline + layer)
##   3) Streams (NPS FeatureServer) + quick plot with watershed
##   4) GRSM datasets (soils/inverts/fish/forest) + watershed plot
##   5) Add streams on top (busy plot)
##   6) NEON data + plot
############################################################

#############################
# 1) PACKAGES & RAW DATA
#############################
library(tidyverse)
library(sf)
library(janitor)
library(geojsonsf)
library(neonstore)
library(neonUtilities)
library(ggnewscale)
library(viridisLite)

# GRSM CSVs
three_pass_locations <- read_csv("/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Fish/Three_Pass/Locations/GRSM_THREE_PASS.csv")
inverts             <- read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export.csv')
invert_locations    <- read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/Documents/Locations.csv')
forest_locations    <- read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Forest_Health/Locations.csv')

# Noland Divide soil metadata (IRMA 705202)
soil_noland  <- read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Soil_Quality/Soil_water_noland_divide_meta.csv') %>% janitor::clean_names()
noland_station <- unique(soil_noland$station_id)

# Park-wide soils (attribute list for matching Veg&Soils plots)
soils <- read_csv("~/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Soil_Quality/Soils.csv")
soils_loc <- soils %>% pull(LOC_NAME) %>% unique()


########################################
# 2) WATERSHEDS (park outline & subbasins)
########################################
# Read watershed polygons (index [2] keeps the desired layer per your workflow)
watershed <- st_read('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Locations/GRSM_WATERSHEDS/GRSM_WATERSHEDS.shp')[2]
plot(watershed)  # quick base plot check

# Dissolve to a single outline (keeps a clean park border for maps)
watershed_outline <- watershed |>
  st_make_valid() |>
  st_union() |>
  st_cast("MULTIPOLYGON") |>
  st_as_sf()

# Quick outline check (base graphics)
plot(st_geometry(watershed_outline), col = NA, border = "black", lwd = 2, axes = TRUE)

#############################
# 3) STREAMS (NPS FeatureServer)
#############################
# Pull named flowlines directly as GeoJSON (keeps GNIS_NAME)
streams <- st_read(
  "https://services1.arcgis.com/fBc8EJBxQRMcHlei/arcgis/rest/services/GRSM_HYDROLOGY/FeatureServer/0/query?where=1=1&outFields=*&f=geojson",
  quiet = TRUE
)

# Minimal prep for colorizing by name
streams_named <- streams %>%
  filter(!is.na(GNIS_NAME), GNIS_NAME != "") %>%
  mutate(color_id = factor(GNIS_NAME))
n <- n_distinct(streams_named$color_id)
pal <- viridis(n, option = "turbo")  # wide palette for many names

# ---- Plot A: Watersheds + Streams only (clean baseline) ----
# Align CRS just in case
streams_named <- st_transform(streams_named, st_crs(watershed))
ggplot() +
  geom_sf(data = watershed_outline, fill = "grey98", color = "black", linewidth = 0.8) +
  geom_sf(data = watershed, fill = NA, color = "grey70", linewidth = 0.3) +
  geom_sf(data = streams_named, aes(color = color_id), linewidth = 0.35, alpha = 0.9) +
  geom_sf_text(data = streams_named, aes(label = GNIS_NAME, color = color_id),
               size = 1.5, check_overlap = TRUE) +
  scale_color_manual(values = pal, guide = "none") +
  coord_sf() +
  labs(title = "GRSM Watersheds + Streams (baseline)") +
  theme_minimal(base_size = 14)

##########################################
# 4) STANDARDIZE GRSM SITE DATA (points)
##########################################
# Noland Divide sites (IRMA 705202)
noland_meta <- read_csv("https://irma.nps.gov/DataStore/DownloadFile/705202?Reference=2304536") %>% clean_names()
noland_soil_coords <- noland_meta %>%
  transmute(
    station_id = location_id,
    latitude,
    longitude,
    datum = lat_lon_datum
  ) %>%
  distinct()

# Veg & Soils plots (ArcGIS REST → GeoJSON → sf)
endpoint <- "https://services1.arcgis.com/fBc8EJBxQRMcHlei/ArcGIS/rest/services/GRSM_VEG_SOIL_PLOTS_VS/FeatureServer/0/query"
u <- paste0(endpoint,
            "?where=1%3D1&outFields=LOC_NAME,PANEL,VS_WATERSHED,GRTS_SITE",
            "&returnGeometry=true&outSR=4326&f=geojson")
pts <- st_read(u, quiet = TRUE)

# Keep coordinates for sites that appear in your soils CSV
soil_coordinates <- pts %>%
  filter(LOC_NAME %in% soils_loc) %>%
  mutate(lon = st_coordinates(geometry)[,1],
         lat = st_coordinates(geometry)[,2]) %>%
  st_drop_geometry() %>%
  dplyr::select(LOC_NAME, lon, lat, PANEL, VS_WATERSHED, GRTS_SITE)

# Re-join to preserve order/dupes
soil_coords_full <- tibble(LOC_NAME = soils_loc) %>%
  left_join(soil_coordinates, by = "LOC_NAME")

# Quick structure checks (unchanged from your script)
str(invert_locations)
str(three_pass_locations)
str(forest_locations)
str(soil_coordinates)
str(noland_soil_coords)

# ---- Standardize to a common schema ----
soils_parkwide_sites <- soil_coordinates %>%
  transmute(dataset = "soils_parkwide",
            site_id = LOC_NAME,
            display_name = LOC_NAME,
            lon, lat,
            datum = "WGS84")  # outSR=4326

noland_sites <- noland_soil_coords %>%
  transmute(dataset = "soils_noland",
            site_id = station_id,
            display_name = station_id,
            lon = longitude, lat = latitude,
            datum = datum) %>%
  distinct()

invert_sites <- invert_locations %>%
  clean_names() %>%
  transmute(dataset = "inverts",
            site_id = station_name,
            display_name = loc_name,
            lon = lon, lat = lat,
            datum = datum)

fish_sites <- three_pass_locations %>%
  clean_names() %>%
  transmute(dataset = "fish_threepass",
            site_id = station_name,
            display_name = paste0(park_pref_name, ", ", section),
            lon = lon, lat = lat,
            datum = "NAD83")

forest_sites <- forest_locations %>%
  clean_names() %>%
  transmute(dataset = "forests",
            site_id = loc_name,
            display_name = loc_name,
            lon = lon, lat = lat,
            datum = datum)

# ---- Transform datums to WGS84 (EPSG:4326) ----
to_wgs84 <- function(df) {
  wgs   <- df %>% filter(datum %in% c("WGS84","WGS 84","WGS_84"))
  nad   <- df %>% filter(datum %in% c("NAD83","NAD 83","NAD_83"))
  other <- df %>% filter(!datum %in% c("WGS84","WGS 84","WGS_84","NAD83","NAD 83","NAD_83"))
  
  wgs_out <- if (nrow(wgs)) {
    st_as_sf(wgs, coords = c("lon","lat"), crs = 4326, remove = FALSE) %>% st_drop_geometry()
  } else wgs
  
  nad_out <- if (nrow(nad)) {
    st_as_sf(nad, coords = c("lon","lat"), crs = 4269, remove = FALSE) %>%
      st_transform(4326) %>%
      mutate(lon = st_coordinates(geometry)[,1],
             lat = st_coordinates(geometry)[,2]) %>%
      st_drop_geometry()
  } else nad
  
  bind_rows(wgs_out, nad_out, other)
}

soils_parkwide_sites <- to_wgs84(soils_parkwide_sites)
noland_sites         <- to_wgs84(noland_sites)
invert_sites         <- to_wgs84(invert_sites)
fish_sites           <- to_wgs84(fish_sites)
forest_sites         <- to_wgs84(forest_sites)

# Combine + sf
sites_all <- bind_rows(
  soils_parkwide_sites,
  noland_sites,
  invert_sites,
  fish_sites,
  forest_sites
) %>% distinct(dataset, site_id, .keep_all = TRUE)

sites_sf <- st_as_sf(sites_all, coords = c("lon","lat"), crs = 4326)
sites_sf$dataset <- factor(
  sites_sf$dataset,
  levels = c("fish_threepass", "forests", "soils_parkwide", "soils_noland", "inverts")
)

######################################################
# 5) PLOT: Watersheds + GRSM datasets (points only)
######################################################
ggplot() +
  geom_sf(data = watershed_outline, fill = NA, color = "black", linewidth = 0.8) +
  geom_sf(data = watershed, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(
    data = sites_sf %>%
      dplyr::mutate(dataset = factor(dataset,
                                     levels = c("fish_threepass","forests","soils_parkwide","soils_noland","inverts"))) %>%
      dplyr::arrange(dataset),
    aes(color = dataset, shape = dataset),
    size = 3, stroke = .75
  ) +
  scale_color_manual(values = c(
    "fish_threepass" = "blue1",
    "inverts"        = "#984ea3",
    "forests"        = "#4daf4a",
    "soils_parkwide" = "orange2",
    "soils_noland"   = "red"
  )) +
  scale_shape_manual(values = c(
    "fish_threepass" = 1,
    "inverts"        = 5,
    "forests"        = 2,
    "soils_parkwide" = 3,
    "soils_noland"   = 3
  )) +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "GRSM Sampling Locations by Dataset (with Watersheds)",
    x = "Longitude", y = "Latitude",
    color = "Dataset", shape = "Dataset"
  ) +
  theme(
    legend.position.inside = c(0.02, 0.98),
    legend.justification   = c("left", "top"),
    legend.background      = element_rect(fill = scales::alpha("white", 0.7), color = NA),
    panel.grid.major       = element_line(color = NA, linewidth = 0.2)
  ) +
  guides(
    color = guide_legend(override.aes = list(fill = NA, alpha = 1, size = 3, stroke = 1),
                         title = "Dataset"),
    shape = guide_legend(override.aes = list(fill = NA, alpha = 1, size = 3, stroke = 1),
                         title = "Dataset")
  )

#################################################
# 6) BUSY PLOT: add Streams over those datasets
#################################################
# Ensure common CRS before overlaying
crs_map <- st_crs(watershed)
streams  <- st_transform(streams,  crs_map)
sites_sf <- st_transform(sites_sf, crs_map)

# Palette per stream name
n_streams   <- n_distinct(na.omit(streams$GNIS_NAME))
pal_streams <- viridis(n_streams, option = "turbo")
names(pal_streams) <- sort(unique(na.omit(streams$GNIS_NAME)))

ggplot() +
  # park
  geom_sf(data = watershed_outline, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = watershed, fill = NA, color = "black", linewidth = 0.3) +
  
  # streams (first color scale)
  geom_sf(
    data = streams %>% filter(!is.na(GNIS_NAME), GNIS_NAME != ""),
    aes(color = GNIS_NAME),
    linewidth = 0.35, alpha = 0.9, show.legend = FALSE
  ) +
  geom_sf_text(
    data = streams %>% filter(!is.na(GNIS_NAME), GNIS_NAME != ""),
    aes(label = GNIS_NAME, color = GNIS_NAME),
    size = 1.5, check_overlap = TRUE, show.legend = FALSE
  ) +
  scale_color_manual(values = pal_streams, guide = "none") +
  
  ggnewscale::new_scale_color() +  # start a new color scale for point datasets
  
  # GRSM datasets (points)
  geom_sf(
    data = sites_sf %>%
      mutate(dataset = factor(dataset,
                              levels = c("fish_threepass","forests","soils_parkwide","soils_noland","inverts"))) %>%
      arrange(dataset),
    aes(color = dataset, shape = dataset),
    size = 3, stroke = .75
  ) +
  scale_color_manual(values = c(
    "fish_threepass" = "blue1",
    "inverts"        = "#984ea3",
    "forests"        = "#4daf4a",
    "soils_parkwide" = "orange2",
    "soils_noland"   = "red"
  ), guide = "legend") +
  scale_shape_manual(values = c(
    "fish_threepass" = 1,
    "inverts"        = 5,
    "forests"        = 2,
    "soils_parkwide" = 3,
    "soils_noland"   = 3
  )) +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "GRSM Sampling Locations + Streams (busy overlay)",
    x = "Longitude", y = "Latitude",
    color = "Dataset", shape = "Dataset"
  ) +
  theme(
    legend.position.inside = c(0.02, 0.98),
    legend.justification   = c("left","top"),
    legend.background      = element_rect(fill = scales::alpha("white", 0.7), color = NA),
    panel.grid.major       = element_line(color = NA, linewidth = 0.2)
  ) +
  guides(
    color = guide_legend(override.aes = list(fill = NA, alpha = 1, size = 3, stroke = 1),
                         title = "Dataset"),
    shape = guide_legend(override.aes = list(fill = NA, alpha = 1, size = 3, stroke = 1),
                         title = "Dataset")
  )

#############################
# 7) NEON DATA + SIMPLE PLOT
#############################
neon_map <- read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Locations/NEON-SiteMap-Table.csv') %>%
  clean_names()
neon_sf <- st_as_sf(neon_map, coords = c("longitude","latitude"), crs = 4326)

# Classify NEON features (aquatic vs terrestrial) from feature_key
neon_sf <- neon_sf %>%
  mutate(feature_key = str_trim(feature_key),
         type = case_when(
           str_detect(feature_key, regex("Sediment|Fish|Riparian|Discharge|Gauge", ignore_case = TRUE)) ~ "neon_aquatic",
           TRUE ~ "neon_terrestrial"
         ))

# Optional: NEON veg plots from VST product
dp <- loadByProduct(site = "GRSM", dpID = "DP1.10098.001", check.size = FALSE)
veg_plots_sf <- dp$vst_perplotperyear %>%
  distinct(plotID, decimalLongitude, decimalLatitude) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# NEON plot (with watershed context)
ggplot() +
  geom_sf(data = watershed_outline, fill = "grey95", color = "black", linewidth = 0.6) +
  geom_sf(data = watershed, fill = "gray95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = dplyr::filter(neon_sf, type == "neon_terrestrial"),
          aes(color = type), shape = 1, size = 2, stroke = .5) +
  geom_sf(data = dplyr::filter(neon_sf, type == "neon_aquatic"),
          aes(color = type), shape = 0, size = 2, stroke = .5) +
  geom_sf(data = veg_plots_sf, color = "forestgreen", shape = 4, size = 2, stroke = 0.5) +
  scale_color_manual(values = c(neon_aquatic = "#00BFFF", neon_terrestrial = "brown"),
                     name = "NEON type") +
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position   = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = scales::alpha("white", 0.8), color = "grey80"),
    panel.grid.major  = element_line(color = "grey90", linewidth = 0.2)
  ) +
  labs(
    title = "NEON Field Installations within GRSM",
    subtitle = "Blue = aquatic; Brown = terrestrial; Green X = VST plot centroids",
    x = "Longitude", y = "Latitude"
  )
