library(nhdplusTools)
library(sf)
library(tidyverse)
library(prism)
library(terra)
library(httr)
library(utils)

#update
drive_path  <- "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU"

# GRSM border
grsm_border0 <- sf::st_read(file.path(drive_path, "Maine/Locations/BOUNDARY_LN/BOUNDARY_LN.shp"))
grsm_border = st_union(grsm_border0)
plot(grsm_border)

# Watershed Data (indlucing HUC 10)
gdb03 <- file.path(drive_path, "/Users/jgradym/Downloads/WBD_03_HU2_GDB/WBD_03_HU2_GDB.gdb") #Region 03 – South Atlantic–Gulf (drains to the Atlantic)
gdb06 <- file.path(drive_path,"WBD_06_HU2_GDB/WBD_06_HU2_GDB.gdb") # Region 06 – Tennessee (drains to the Mississippi)

                   
stream_dir <- file.path(drive_path, "/Maine/Locations/Streams")
st_layers(gdb03)

# Load HUC10 (Watershed) polygons for both regions (watershed boundary dataset = wbd)
wbd10_03 <- st_read(gdb03, "WBDHU10", quiet = TRUE)
wbd10_06 <- st_read(gdb06, "WBDHU10", quiet = TRUE)
wbd10    <- dplyr::bind_rows(wbd10_03, wbd10_06)

# Clip to your bbox
bbox_grsm <- st_as_sfc(st_bbox(c(xmin=-84.01389, ymin=35.42896, xmax=-82.99804, ymax=35.84139), crs=4326))
bbox_wbd  <- st_transform(bbox_grsm, st_crs(grsm_border))

# Clip to GRSM bounding box
bbox_wbd <- st_transform(bbox_grsm, st_crs(wbd10_06 ))

grsm_watershed <- st_intersection(wbd10, bbox_wbd) |>
  dplyr::select(huc10, name, states, areasqkm) 
grsm_watershed$name


# filter out non park watershed names
grsm_watershed = grsm_watershed %>%
  filter(!name %in% c("Sinking Creek-Tennessee River", "French Broad River", "Gulf Fork Big Creek", "Richland Creek-Pigeon River",
                      "Cheoah River", "Alarka Creek-Little Tennessee River", "Middle Tuckasegee River","Upper Tellico Lake",
                      "Lower Tellico Lake")) 
st_write(grsm_watershed, file.path(drive_path, "Maine/Spatial/Watershed/GRSM_watershed.gpkg"), delete_dsn = TRUE)
# Plot with labels
watershed_labels <- st_point_on_surface(watershed_grsm )

ggplot() +
  geom_sf(data = grsm_border,  fill = NA, color = "gray50", linewidth = 1) +
  geom_sf(data = grsm_watershed, fill = NA, color = "blue", linewidth = 0.4) +
  geom_sf(data = grsm_border, color = "gray40", linewidth = 1.1) +
  geom_sf_label(
    data = watershed_labels, aes(label = name),
    fill = "white", color = "blue",
    label.size = 0,         # no box outline
    label.r = unit(0.15, "lines"), # subtle rounding
    size = 4) +
  coord_sf() +
  labs(title = "GRSM HUC10 Watersheds with Park Boundary") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))


#----------------------------------
#------ Add streams --------------
#---------------------------------


# --- NHDPlus HR download + flowlines ---

# USGS source for watersheds and stream hyrdrology of GRSM:
#https://www.usgs.gov/national-hydrography/access-national-hydrography-products
#https://apps.nationalmap.gov/downloader/#/

dir.create(stream_dir, recursive = TRUE, showWarnings = FALSE)
hu4s <- c("0307","0601")
download_nhdplushr(nhd_dir = stream_dir, hu_list = hu4s, download_files = TRUE)

gdb_paths <- list.dirs(stream_dir, recursive = TRUE, full.names = TRUE) |>
  grep("\\.gdb$", x = _, value = TRUE)

flow_list <- lapply(gdb_paths, function(g) {
  st_read(g, "NHDFlowline", quiet = TRUE) |>
    select(any_of(c("GNIS_Name","gnis_name","FType","FCode","ReachCode","reachcode","NHDPlusID","nhdplusid","geometry")))
})

# check crdds
imap(flow_list, ~ { message(.y, ": ", st_crs(.x)$input); NULL })
# drop the vertical crs

flows_norm <- map(flow_list, ~ .x %>%
                     sf::st_zm(drop = TRUE) %>%
                     sf::st_transform(target_crs) %>%
                     dplyr::rename_with(tolower)
)

target_crs <- 26917
unique(flows_norm$ftype)
flows <- dplyr::bind_rows(flows_norm) %>%
  sf::st_make_valid() %>%
  sf::st_simplify(dTolerance = 50)  # meters in UTM

#
# --- Clip streams to GRSM HU10 mask ---
grsm_stream <- st_intersection(flows, st_transform(grsm_mask, st_crs(flows))) |>
  filter(!is.na(gnis_name), gnis_name != "")

grsm_mask <- st_union(st_make_valid(wbd10_grsm))
grsm_stream <- st_intersection(flows, grsm_mask) |>
  filter(!is.na(gnis_name), gnis_name != "")

grsm_stream <- st_intersection(
  st_transform(flows, 26917),
  st_transform(grsm_mask, 26917)
) |> dplyr::filter(!is.na(gnis_name), gnis_name != "")
grsm_stream = grsm_stream[1]

#st_write(grsm_stream, file.path(drive_path, "Maine/Spatial/Streams/GRSM_streams.gpkg"), delete_dsn = TRUE)


watershed_label <- st_point_on_surface(grsm_watershed )

# --- Label points (midpoint per named stream) ---
lines_longest <- grsm_stream |>
  st_transform(26917) |>
  st_collection_extract("LINESTRING") |>
  dplyr::filter(!is.na(gnis_name), gnis_name != "") |>
  dplyr::group_by(gnis_name) |>
  dplyr::summarise(shape = st_line_merge(st_union(shape)), .groups = "drop") |>
  st_cast("LINESTRING", warn = FALSE) |>
  dplyr::group_by(gnis_name) |>                                   # <- re-group here
  dplyr::ungroup()

pts_utm <- st_cast(st_line_sample(lines_longest, sample = 0.5), "POINT")

label_pts <- st_sf(
  gnis_name = lines_longest$gnis_name,
  shape     = pts_utm,
  crs       = st_crs(lines_longest)
) |>
  st_transform(st_crs(grsm_stream))

write

# --- Your plot (unchanged aesthetics) ---
# Reduce stream names  to those in data files
inverts <- readr::read_csv(file.path(drive_path, "Maine/Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export_locations.csv")) # locations added
inverts$StreamName


# Watersheds + Streams that inverts were sampled in
ggplot() +
  geom_sf(data = grsm_border,  fill = NA, color = "gray50", linewidth = 1) +
  geom_sf(data = grsm_watershed, fill = NA, color = "blue", linewidth = 0.4) +
  
  geom_sf(data = grsm_stream %>% filter(!is.na(gnis_name), gnis_name %in% inverts$StreamName),
          aes(color = gnis_name), linewidth = 0.4, show.legend = F) +
  geom_sf_text(data = label_pts %>% filter(gnis_name %in% inverts$StreamName),
               aes(label = gnis_name, color = gnis_name),
               size = 2.6, check_overlap = TRUE, show.legend = F) +
  geom_sf_label(data = watershed_label, aes(label = name),fill = "white", fontface = "bold", sieze = 3, color = "blue",label.size = 0,
                label.r = unit(0.15, "lines"), size = 3) +
  coord_sf() +
  ggtitle("GRSM Watersheds (blue) and Streams with Invert Sampling (color)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank())


#------ Streams with fish -----
three_pass <- read_csv(file.path(drive_path, 'Maine/Data/Aquatics_Fish/Three_Pass/Summary_data/GRSM_Fish_3-Pass_Summary_with_loc.csv'))

ggplot() +
  geom_sf(data = grsm_border,  fill = NA, color = "gray50", linewidth = 1) +
  geom_sf(data = grsm_watershed, fill = NA, color = "blue", linewidth = 0.4) +
  
  geom_sf(data = grsm_stream %>% filter(!is.na(gnis_name), gnis_name %in% three_pass$Stream),
          aes(color = gnis_name), linewidth = 0.4, show.legend = F) +
  geom_sf_text(data = label_pts %>% filter(gnis_name %in% three_pass$Stream),
               aes(label = gnis_name, color = gnis_name),
               size = 2.6, check_overlap = TRUE, show.legend = F) +
  geom_sf_label(data = cent10, aes(label = name),fill = "white", fontface = "bold", sieze = 3, color = "blue",label.size = 0,
                label.r = unit(0.15, "lines"), size = 3) +
  coord_sf() +
  ggtitle("GRSM Watersheds (blue) and Streams with Fish Sampling (color)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank())

#------------------------ Prism Data -------------



# Quick outline check (base graphics)
plot(grsm_border)


# Download PRISM data to a folder
DL  <- file.path(drive_path, 'Maine/Data/Climate/prism')
OUT <- file.path(DL, "annual")
dir.create(DL,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

yrs <- 1980:2023
mos <- 1:12
elements <- c("ppt","tmean")   # add "tmin","tmax" later if you want


grsm <- st_as_sf(st_union(st_make_valid(watershed_outline)))

# --- 1) download monthly 800 m PRISM (tiny, explicit loop; no functions) ----
for (el in elements) {
  for (yy in yrs) for (mm in mos) {
    ym <- sprintf("%04d%02d", yy, mm)
    base <- sprintf("https://services.nacse.org/prism/data/get/us/800m/%s/%s", el, ym)
    
    # clean any half-made target dir to avoid "already exists" bug
    old_dirs <- list.dirs(DL, recursive = FALSE, full.names = TRUE)
    old_dirs <- old_dirs[grepl(paste0("^PRISM_", el, ".*_", ym, "$"), basename(old_dirs))]
    if (length(old_dirs)) unlink(old_dirs, recursive = TRUE, force = TRUE)
    
    # download zip to temp, then unzip to a stable folder name
    tf <- tempfile(fileext = ".zip")
    resp <- GET(base, write_disk(tf, overwrite = TRUE), timeout(300))
    stop_for_status(resp)
    
    cd  <- headers(resp)[["content-disposition"]]
    zn  <- if (!is.null(cd) && grepl("filename=", cd))
      sub('.*filename="?([^";]+).*', "\\1", cd) else paste0("PRISM_", el, "_", ym, ".zip")
    if (!grepl("\\.zip$", zn)) zn <- paste0(zn, ".zip")
    out_zip  <- file.path(DL, zn)
    file.rename(tf, out_zip)
    
    out_dir <- file.path(DL, sub("\\.zip$", "", zn))
    if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
    unzip(out_zip, exdir = out_dir)
    unlink(out_zip)
    message("OK: ", el, " ", ym)
  }
}

# 


# Load data and get a raster stack
ppt_files <- list.files(
  DL,
  pattern = "prism_ppt_us_30s_[0-9]{6}\\.tif$",
  recursive = TRUE, full.names = TRUE, ignore.case = TRUE
)

tmean_files <- list.files(
  DL,
  pattern = "prism_tmean_us_30s_[0-9]{6}\\.tif$",
  recursive = TRUE, full.names = TRUE, ignore.case = TRUE
)

order_by_ym <- function(x){
  ym <- str_extract(basename(x), "[0-9]{6}")
  x[order(ym)]
}
ppt_files   <- order_by_ym(ppt_files)
tmean_files <- order_by_ym(tmean_files)

length(ppt_files); length(tmean_files)   # should both be > 0

# Precipitation
ppt_stack   <- rast(ppt_files)      # mm / month
# Temperature
tmean_stack <- rast(tmean_files)    # °C (GeoTIFFs from NACSE are in °C)

stopifnot(exists("watershed_outline"))
grsm_ll <- st_transform(st_as_sf(st_union(st_make_valid(watershed_outline))),
                        crs(ppt_stack))

ppt_grsm   <- mask(crop(ppt_stack,   vect(grsm_ll)), vect(grsm_ll))
tmean_grsm <- mask(crop(tmean_stack, vect(grsm_ll)), vect(grsm_ll))

# figure out which monthly layers belong to each year from file names
ym  <- str_extract(basename(ppt_files), "[0-9]{6}")
yrs <- substr(ym, 1, 4)
idx_by_year <- split(seq_along(yrs), yrs)
years_vec <- names(idx_by_year)

# annual MEAN (ppt = mean of 12 months; tmean = mean of 12 months)
ppt_annual_mean <- rast(lapply(years_vec, function(yy) app(ppt_grsm[[ idx_by_year[[yy]] ]],
                                                           mean, na.rm = TRUE)))
tmean_annual_mean <- rast(lapply(years_vec, function(yy) app(tmean_grsm[[ idx_by_year[[yy]] ]],
                                                             mean, na.rm = TRUE)))

names(ppt_annual_mean)   <- paste0("ppt_mean_",   years_vec)
names(tmean_annual_mean) <- paste0("tmean_mean_", years_vec)

# Write rasters - mean temp and precip
writeRaster(ppt_annual_mean, file.path(drive_path, 'Maine/Data/Climate/PRISM_GRSM_precip_annual_mean_stack.tif'))

writeRaster(tmean_annual_mean, file.path(drive_path, 'Maine/Data/Climate/PRISM_GRSM_temp_annual_mean_stack.tif'))


plot(ppt_annual_mean$ppt_mean_1980)
plot(tmean_annual_mean$tmean_mean_1980)



# brown → green (for precipitation)
plot(ppt_annual_mean$ppt_mean_1980,
     col = hcl.colors(20, "Terrain", rev = TRUE),   # earthy brown–green
     main = "Annual Precipitation 1980 (mm)")

# red → blue (for temperature)
plot(tmean_annual_mean$tmean_mean_1980,
     col = rev(hcl.colors(20, "RdBu")),             # red–blue diverging
     main = "Annual Mean Temperature 1980 (°C)")

