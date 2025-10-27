
library(prism)
library(terra)
library(httr)
library(utils)
library(dplyr)
library(stringr)
library(sf)

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



DL  <- "/Users/jgradym/Desktop/prism"
OUT <- file.path(DL, "annual")
dir.create(DL,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

yrs <- 1980:2023
mos <- 1:12
elements <- c("ppt","tmean")   # add "tmin","tmax" later if you want

# GRSM polygon you already built earlier
stopifnot(exists("watershed_outline"))
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

library(stringr)
library(terra)
library(sf)

DL  <- "/Users/jgradym/Desktop/prism"
OUT <- file.path(DL, "annual"); dir.create(OUT, showWarnings = FALSE)

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


writeRaster(ppt_annual_mean, 
            '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Climate/PRISM_GRSM_precip_annual_mean_stack.tif')

writeRaster(tmean_annual_mean, 
            '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Climate/PRISM_GRSM_temp_annual_mean_stack.tif')


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
