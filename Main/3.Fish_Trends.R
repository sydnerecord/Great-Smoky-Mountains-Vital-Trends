# ======================================================================
# Libraries
# ======================================================================
library(RColorBrewer)  # palettes
library(readxl)        # read_xlsx() for NCBI file
library(viridisLite)   # color options
library(sf)            # spatial analysis
library(tidyverse)     # data import (readr), manipulation (dplyr), plotting (ggplot2), strings (stringr)
library(scales)



theme_plot <- function(legend = TRUE) {
  ggplot2::theme(
    panel.grid        = element_blank(),                   # remove major/minor grids
    aspect.ratio      = .75,                               # tall aspect for time-series
    axis.text         = element_text(size = 18, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title        = element_text(size = 18),
    axis.title.y      = element_text(margin = margin(r = 10)),
    axis.title.x      = element_text(margin = margin(t = 10)),
    axis.title.x.top  = element_text(margin = margin(b = 5)),
    plot.title        = element_text(size = 18, face = "plain", hjust = 0),  
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background  = element_blank(),
    strip.background  = element_blank(),
    legend.text       = element_text(size = 15),
    legend.title      = element_text(size = 18),
    legend.position   = if (legend) "right" else "none",   # toggle legend on/off
    text              = element_text(family = "Helvetica") # global font
  )
}

# 1) Read spatial data layers
# ----------------------------

# Update path
drive_path  <- "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU"

# fish data
three_pass <- read_csv(file.path(drive_path, 'Maine/Data/Aquatics_Fish/Three_Pass/Summary_data/GRSM_Fish_3-Pass_Summary_with_coordinates.csv'))
str(three_pass)
# Location Categories - Stream and LOC_NAME or Code, Watershed and Streamname
# Response variables :
#YOYpop – Count of young-of-year individuals (juveniles) captured per sample.
#ADTpop – Count of adult individuals captured per sample.
#TOTpop – Total number of individuals (YOY + adult) captured.
#YOYDens – Density of young-of-year fish (individuals per unit area sampled).
#ADTDens – Density of adult fish (individuals per unit area).
#TOTDens – Total fish density (sum of YOY and adult densities).
#YOYmnwt – Mean individual weight of young-of-year fish (g).
#ADTmnwt – Mean individual weight of adult fish (g).
#TOTmnwt – Mean individual weight across all size classes (g).
#YOYBiom – Biomass of young-of-year fish (g m⁻² or g, depending on standardization).
#ADTBiom – Biomass of adult fish.
#TOTBiom – Total fish biomass (YOY + adult)


# --- Species lookup table -------------------------------------
species_lookup <- tibble::tribble(
  ~code, ~common, ~scientific,
  "BEC","Bigeye Chub","Hybopsis amblops",
  "BRH","Black Redhorse","Moxostoma duquesnei",
  "GFD","Greenfin Darter","Etheostoma chlorobranchium",
  "GID","Gilt Darter","Percina evides",
  "NHS","Northern Hog Sucker","Hypentelium nigricans",
  "RIC","River Chub","Nocomis micropogon",
  "RLD","Redline Darter","Etheostoma rufilineatum",
  "ROB","Rock Bass","Ambloplites rupestris",
  "SMB","Smallmouth Bass","Micropterus dolomieu",
  "STR","Striped Shiner","Luxilus chrysocephalus",
  "TES","Telescope Shiner","Notropis telescopus",
  "TNS","Tennessee Shiner","Notropis leuciodus",
  "WPS","Warpaint Shiner","Luxilus coccogenis",
  "WTS","Whitetail Shiner","Cyprinella galactura",
  "BAD","Banded Darter","Etheostoma zonale",
  "RBT","Rainbow Trout","Oncorhynchus mykiss",
  "WHS","White Sucker","Catostomus commersonii",
  "BND","Blacknose Dace","Rhinichthys atratulus",
  "WND","Wounded Darter","Etheostoma vulneratum",
  "CKC","Creek Chub","Semotilus atromaculatus",
  "SMT","Smoky Madtom","Noturus baileyi",
  "CIT","Citico Darter","Etheostoma sitikuense",
  "BNT","Brown Trout","Salmo trutta",
  "RSD","Rosyside Dace","Clinostomus funduloides",
  "SAS","Saffron Shiner","Notropis rubricroceus",
  "FHM","Fathead Minnow","Pimephales promelas",
  "LND","Longnose Dace","Rhinichthys cataractae",
  "FTD","Fantail Darter","Etheostoma flabellare",
  "TND","Tennessee Dace","Chrosomus tennesseensis",
  "BKT","Brook Trout","Salvelinus fontinalis"
)


species_lookup %>%
  dplyr::filter(code %in% species_min) %>%
  dplyr::select(code, common, scientific)

#----------- Add spatial layers -------------------

# Spatial layers (consistent format)
grsm_watershed = st_read(file.path(drive_path, "Maine/Spatial/Watershed/GRSM_watershed.gpkg"))[2]

grsm_streams = st_read(file.path(drive_path, "Maine/Spatial/Streams/GRSM_streams.gpkg"))

grsm_border <- sf::st_read(file.path(drive_path, "Maine/Spatial/BOUNDARY_LN/BOUNDARY_LN.shp"))[6] %>%
  sf::st_transform(sf::st_crs(grsm_watershed))


#-------- Analyze Sampling -------
# --- Optional: restrict to species with ≥ min_sample adult-positive years ---
min_sample <- 5
species_min <- three_pass %>%
  dplyr::mutate(Year = lubridate::year(Date)) %>%
  dplyr::group_by(Species, Year) %>%
  dplyr::summarise(has_adults = any(is.finite(ADTDens) & ADTDens > 0), .groups = "drop") %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(n_years_pos = sum(has_adults), .groups = "drop") %>%
  dplyr::filter(n_years_pos >= min_sample) %>%
  dplyr::pull(Species)

# --- Define “sampling events” as unique LOC_NAME × Date per Stream ---
base_samples <- three_pass %>%
  dplyr::filter(Species %in% species_min) %>%                 # drop if you want all species
  dplyr::mutate(
    Stream = Stream,
    Year   = lubridate::year(Date),
    Month  = factor(lubridate::month(Date, label = TRUE, abbr = TRUE), levels = month.abb),
    doy    = lubridate::yday(Date),
    Day    = as.Date(Date)
  ) %>%
  dplyr::distinct(Stream, LOC_NAME, Day, Year, Month, doy)

# keep streams with ≥5 sampled years
streams_eligible <- base_samples %>%
  dplyr::group_by(Stream) %>%
  dplyr::summarise(n_years = dplyr::n_distinct(Year), .groups = "drop") %>%
  dplyr::filter(n_years >= 5) %>%
  dplyr::pull(Stream)

base_5yr <- base_samples %>% dplyr::filter(Stream %in% streams_eligible)

# --- (1) Stacked bar: count of sampling events by month, colored by stream ---
ggplot(base_5yr, aes(x = Month, fill = Stream)) +
  geom_bar(position = "stack") +
  scale_fill_viridis_d(option = "turbo") +
  labs(x = "Month", y = "Number of Samples", fill = "Stream",
       title = "Fish Three-Pass: Overall Sampling Seasonality (≥5-Year Streams)") +
  theme_plot() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

# --- (2) Scatter: day-of-year vs year, colored by stream ---
ggplot(base_5yr, aes(y = Year, x = doy, color = Stream)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_continuous(
    breaks = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349),
    labels = month.abb
  ) +
  scale_y_continuous(breaks = seq(1980, 2025, by = 5)) +
  scale_color_viridis_d(option = "turbo") +
  labs(x = "Sampling Month", y = "Year", color = "Stream",
       title = "Fish Three-Pass: Sampling Dates by Year (≥5-Year Streams)") +
  theme_plot() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

# ----------------------------
# Identify eligible species with minimum years sampling 
# ----------------------------


# ----------------------------
# Spatial analysis
# ----------------------------
# Use one projection everywhere (NAD83 / UTM 17N)
target_crs <- 26917

# --- Ensure common CRS -------------------------------------------------
grsm_watershed <- st_transform(grsm_watershed, target_crs)   # must have a 'name' col
grsm_border    <- st_transform(grsm_border,    target_crs)
grsm_streams   <- st_transform(grsm_streams,   target_crs)

# --- One coordinate per stream from fish data --------------------------
fish_locs_sf <- three_pass %>%
  filter(is.finite(LAT), is.finite(LON)) %>%
  arrange(Stream, Date) %>%
  group_by(Stream) %>%
  slice_head(n = 1) %>%
  st_as_sf(coords = c("LON","LAT"), crs = 4326, remove = FALSE) %>%
  st_transform(target_crs) %>%
  transmute(Location = Stream, name_norm = str_squish(str_to_lower(Stream)), geometry)

# --- Match hydro lines to those streams (normalize both sides) ----------
streams_matched_sf <- grsm_streams %>%
  filter(!is.na(gnis_name), gnis_name != "") %>%
  mutate(name_norm = str_squish(str_to_lower(gnis_name))) %>%
  semi_join(fish_locs_sf %>% st_drop_geometry() %>% select(name_norm), by = "name_norm")

# --- Labels: one point per stream & per watershed ----------------------
stream_labels_sf <- streams_matched_sf %>%
  group_by(gnis_name) %>%
  summarize(do_union = TRUE, .groups = "drop") %>%
  mutate(geometry = st_point_on_surface(geom))

# If your watershed name column isn’t 'name', change it below
watershed_labels_sf <- grsm_watershed %>%
  mutate(geometry = st_point_on_surface(geom)) %>%
  select(name, geometry)

# --- Plot --------------------------------------------------------------
# With watershed
ggplot() +
  geom_sf(data = grsm_border,     fill = NA, color = "gray50", linewidth = 1) +
  geom_sf(data = grsm_watershed,  fill = NA, color = "blue", linewidth = 0.6) +
  geom_sf_label(data = watershed_labels_sf, aes(label = name),
                fill = "white", color = "blue", fontface = "bold",
                label.size = 0, label.r = grid::unit(0.15, "lines"), size = 3) +
  geom_sf(data = streams_matched_sf, aes(color = gnis_name), linewidth = 0.75, alpha = 0.9) +
  geom_sf(data = fish_locs_sf, aes(color = Location), size = 4) +
  geom_sf_text(data = stream_labels_sf,   aes(label = gnis_name, color = gnis_name),
               size = 3, check_overlap = TRUE) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  coord_sf() +
  labs(title = "GRSM Fish Sampling Locations") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

# Without watershed
ggplot() +
  geom_sf(data = grsm_border,     fill = NA, color = "gray50", linewidth = 1) +
  geom_sf(data = streams_matched_sf, aes(color = gnis_name), linewidth = 0.75, alpha = 0.9) +
  geom_sf(data = fish_locs_sf, aes(color = Location), size = 4) +
  geom_sf_text(data = stream_labels_sf,   aes(label = gnis_name, color = gnis_name),
               size = 3, check_overlap = TRUE) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  coord_sf() +
  labs(title = "GRSM Fish Sampling Locations") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())


# Fit per-(Species, Stream) linear models
# ----------------------------


min_sample <- 5 # same nubmer of species if min_sample = 2

species_min <- three_pass %>%
  dplyr::mutate(Year = lubridate::year(Date)) %>%
  dplyr::group_by(Species, Stream, Year) %>%
  dplyr::summarise(has_adults = any(is.finite(ADTDens) & ADTDens > 0),
                   .groups = "drop_last") %>%
  dplyr::summarise(n_years_pos = sum(has_adults), .groups = "drop") %>%
  dplyr::summarise(max_years_pos = max(n_years_pos), .by = Species) %>%
  dplyr::filter(max_years_pos >= min_sample) %>%
  dplyr::pull(Species)

# Which species with multiyear sampling?
species_lookup %>%
  dplyr::filter(code %in% species_min) %>%
  dplyr::select(code, common, scientific)

lm_species_stream <- three_pass %>%
  mutate(Year = as.integer(format(Date, "%Y"))) %>%
  filter(Species %in% species_min, is.finite(ADTDens), ADTDens > 0) %>%
  group_by(Species, Stream, Code, Year) %>%
  summarise(site_adult_density = mean(ADTDens, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, Stream, Year) %>%
  summarise(mean_adult_density = mean(site_adult_density, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, Stream) %>%
  filter(n_distinct(Year) >= min_sample) %>%
  do({
    m  <- lm(mean_adult_density ~ Year, data = .)
    cs <- summary(m)$coefficients
    tibble(
      slope   = cs["Year", "Estimate"],
      p_value = cs["Year", "Pr(>|t|)"],
      r_sq    = summary(m)$r.squared
    )
  }) %>%
  ungroup()

unique(lm_species_stream$Species)
# "BKT" "BNT" "RBT" "ROB" "SMB"


# ----------------------------
# 7) Plot one species (example: RBT)
# ----------------------------
species_names <- tibble::tibble(
  Species = c("BKT", "BNT", "RBT", "ROB", "SMB", "BND", "LND"),
  Common  = c("Brook Trout", "Brown Trout", "Rainbow Trout",
              "Rock Bass", "Smallmouth Bass",
              "Blacknose Dace", "Longnose Dace")
)

sp <- "BKT"
# Toggle: TRUE = gray out nonsignificant slopes, FALSE = color all
gray_out <- F
p_thresh <- 0.05


sp_name <- "Brook Trout"


# Create a plotting dataframe for one species, with slope filtering logic
dat <- lm_species_stream_sf %>%
  dplyr::left_join(species_names, by = "Species") %>%
  dplyr::filter(Species == sp) %>%
  dplyr::mutate(
    fill_val = dplyr::if_else(
      gray_out & p_value >= p_thresh,
      NA_real_,
      slope
    )
  )

# Plot
ggplot() +
  geom_sf(data = grsm_border,    fill = NA, color = "gray40", linewidth = 1) +
  geom_sf(data = grsm_watershed, fill = NA, color = "grey70", linewidth = 0.3) +
  geom_sf(data = streams_matched_sf, color = "grey70", linewidth = 0.6, alpha = 0.9) +
  geom_sf_text(data = stream_labels_sf, aes(label = gnis_name),
               color = "grey40", size = 3, check_overlap = TRUE) +
  geom_sf(data = dat, aes(fill = fill_val),
          shape = 21, color = "black", size = 5, stroke = 0.3) +
  scale_fill_gradientn(colors = c("#1E90FF", "white", "red"),
                       limits = c(-max_abs, max_abs),
                       na.value = "grey80", name = "Slope") +
  coord_sf() +
  labs(title = paste0("Adult Density Trend — ", sp_name)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())



#-------- loop ------


# ============================================================
# Save adult density trend maps (slope) for each species as PDF
#   - First species: all points colored
#   - Remaining species: nonsignificant points grayed out
# ============================================================

pdf_out <- "~/Downloads/fish_trend_maps_by_species.pdf"
p_thresh <- 0.05
species_list <- sort(unique(lm_species_stream_sf$Species))

# ============================================================
# Save adult density trend maps (slope) for each species as PDF
#   - All species fully colored (no gray-out)
# ============================================================

pdf_out <- "~/Downloads/fish_trend_maps_by_species_check.pdf"
p_thresh <- 0.05
species_list <- sort(unique(lm_species_stream_sf$Species))
gray_out <- FALSE   # toggle: FALSE = color all, TRUE = gray non-significant

pdf(file = pdf_out, width = 10, height = 7)
for (sp in species_list) {
  dat <- lm_species_stream_sf %>%
    dplyr::left_join(species_names, by = "Species") %>%
    dplyr::filter(Species == sp) %>%
    dplyr::mutate(
      fill_val = if (gray_out)
        ifelse(p_value < p_thresh, slope, NA_real_)
      else slope
    )
  
  sp_name <- unique(dat$Common)
  max_abs <- max(abs(dat$slope), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  p <- ggplot() +
    geom_sf(data = grsm_border, fill = NA, color = "gray40", linewidth = 1) +
    geom_sf(data = grsm_watershed, fill = NA, color = "grey70", linewidth = 0.2) +
    geom_sf(data = streams_matched_sf, color = "grey70", linewidth = 0.6, alpha = 0.9) +
    geom_sf_text(data = stream_labels_sf, aes(label = gnis_name),
                 color = "grey40", size = 3, check_overlap = TRUE) +
    geom_sf(data = dat, aes(fill = fill_val),
            shape = 21, color = "black", size = 4, stroke = 0.3) +
    scale_fill_gradientn(colors = c("#1E90FF", "white", "red"),
                         limits = c(-max_abs, max_abs),
                         na.value = "grey80", name = "Slope") +
    coord_sf() +
    labs(title = paste0("Adult Density Trend — ", sp_name, " (", sp, ")")) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank())
  
  print(p)
}
dev.off()

cat("PDF saved to:", normalizePath(pdf_out), "\n")


# grayed out if non signicant

# ============================================================
# Save adult density trend maps (slope) for each species as PDF
#   - All species fully colored (no gray-out)
# ============================================================
# update path
pdf_out <- "~/Downloads/fish_trend_maps_by_species_gray.pdf"
p_thresh <- 0.05
species_list <- sort(unique(lm_species_stream_sf$Species))
gray_out <- T   # toggle: FALSE = color all, TRUE = gray non-significant

pdf(file = pdf_out, width = 10, height = 7)
for (sp in species_list) {
  dat <- lm_species_stream_sf %>%
    dplyr::left_join(species_names, by = "Species") %>%
    dplyr::filter(Species == sp) %>%
    dplyr::mutate(
      fill_val = if (gray_out)
        ifelse(p_value < p_thresh, slope, NA_real_)
      else slope
    )
  
  sp_name <- unique(dat$Common)
  max_abs <- max(abs(dat$slope), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  p <- ggplot() +
    geom_sf(data = grsm_border, fill = NA, color = "gray40", linewidth = 1) +
    geom_sf(data = grsm_watershed, fill = NA, color = "grey70", linewidth = 0.2) +
    geom_sf(data = streams_matched_sf, color = "grey70", linewidth = 0.6, alpha = 0.9) +
    geom_sf_text(data = stream_labels_sf, aes(label = gnis_name),
                 color = "grey40", size = 3, check_overlap = TRUE) +
    geom_sf(data = dat, aes(fill = fill_val),
            shape = 21, color = "black", size = 4, stroke = 0.3) +
    scale_fill_gradientn(colors = c("#1E90FF", "white", "red"),
                         limits = c(-max_abs, max_abs),
                         na.value = "grey80", name = "Slope") +
    coord_sf() +
    labs(title = paste0("Adult Density Trend — ", sp_name, " (", sp, ")")) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank())
  
  print(p)
}
dev.off()

cat("PDF saved to:", normalizePath(pdf_out), "\n")


# ===========================================
# Fish: trends per species (each stream plotted)
# - Exclude zeros
# - Average per Code–Year, then Stream–Year
# - Require ≥5 years per stream for plotting 
# loop through species
# ===========================================


# Species to plot
species_list <- sort(unique(lm_species_stream_sf$Species))

# Year range for axes
yrs <- range(fish_species_stream_year$Year, na.rm = TRUE)

# Output path
pdf_path <- "~/Downloads/adult_density_all_species.pdf" # update to your computer

# Open a multi-page PDF device
pdf(pdf_path, width = 8, height = 6, useDingbats = FALSE)  # useDingbats=FALSE for clean text rendering

for (sp in species_list) {
  
  # Lookup common name
  sp_common <- species_names %>%
    filter(Species == sp) %>%
    pull(Common)
  
  title_text <- if (length(sp_common) == 0) sp else paste0(sp, " — ", sp_common)
  
  # Build plot
  p <- ggplot(
    fish_species_stream_year %>%
      filter(Species == sp, n_years_stream >= 5),
    aes(x = Year, y = mean_adult_density, color = Stream, group = Stream)
  ) +
    geom_line() +
    geom_point(size = 2) +
    geom_smooth(method = "lm", linetype = "dashed", se = FALSE, linewidth = 0.75) +
    scale_y_log10(
      breaks = log_breaks(base = 10),
      labels = label_number(accuracy = 1, trim = TRUE)
    ) +
    scale_x_continuous(
      breaks = seq(yrs[1], yrs[2], by = 5),
      limits = c(yrs[1], yrs[2])
    ) +
    scale_color_viridis_d(option = "turbo") +
    labs(x = "Year", y = "Adult Density", color = "Stream") +
    theme_plot() +
    ggtitle(paste0("Adult Density — ", title_text))
  
  print(p)  # prints to the open PDF device
}

dev.off()  # close the PDF device


# 1) Build per-species, per-stream, per-year means (exclude 0s first)
fish_species_stream_year <- three_pass %>%
  mutate(Year = as.integer(format(Date, "%Y"))) %>%
  filter(is.finite(ADTDens), ADTDens > 0) %>%
  group_by(Species, Stream, Code, Year) %>%
  summarise(site_adult_density = mean(ADTDens, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, Stream, Year) %>%
  summarise(mean_adult_density = mean(site_adult_density, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, Stream) %>%
  mutate(n_years_stream = n_distinct(Year)) %>%
  ungroup()

# 2) Pick a species (change "RBT" as needed)
sp <- "RBT"

yrs <- fish_species_stream_year %>%
  filter(Species == sp) %>%
  summarise(rng = range(Year, na.rm = TRUE)) %>%
  pull(rng)

# 3) Plot (streams = colors; dashed LM fit; log10 y with clean ticks)
ggplot(
  fish_species_stream_year %>%
    filter(Species == sp, n_years_stream >= 5),
  aes(x = Year, y = mean_adult_density, color = Stream, group = Stream)) +
  geom_line() +
  geom_point(size = 2) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.75) +
  scale_y_log10(breaks = scales::log_breaks(base = 10),
                labels = scales::label_number(accuracy = 1, trim = TRUE)) +
  scale_x_continuous(breaks = seq(yrs[1], yrs[2], by = 5),
                     limits = c(yrs[1], yrs[2])) +
  scale_color_viridis_d(option = "turbo") +
  labs(x = "Year", y = "Adult Density", color = "Stream") +
  theme_plot() +
  ggtitle(paste0("Adult Density — ", sp))





# ================================================================
# Fish: per-species adult density trends (each stream = one line)
# - Exclude zeros
# - Average per Code–Year, then per Stream–Year
# - Require ≥5 sampled years per stream
# - No linear model overlay (lines + points only)
# - Add full species name in title
# - Output: one PDF page per species
# ================================================================



# 1) Prepare per-species, per-stream, per-year data
fish_species_stream_year <- three_pass %>%
  mutate(Year = as.integer(format(Date, "%Y"))) %>%
  filter(is.finite(ADTDens), ADTDens > 0) %>%
  group_by(Species, Stream, Code, Year) %>%
  summarise(site_adult_density = mean(ADTDens, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, Stream, Year) %>%
  summarise(mean_adult_density = mean(site_adult_density, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, Stream) %>%
  mutate(n_years_stream = n_distinct(Year)) %>%
  ungroup()

# 0) Save loop
pdf_out <- path.expand("~/Downloads/fish_trends_per_species_streams.pdf") #update output path
dir.create(dirname(pdf_out), showWarnings = FALSE, recursive = TRUE)

# 1) Open PDF device
pdf(file = pdf_out, width = 11, height = 7, useDingbats = FALSE)
on.exit(dev.off(), add = TRUE)

# 2) Loop and PLOT at least something
species_list <- sort(unique(fish_species_stream_year$Species))
n_pages <- 0L

for (sp in species_list) {
  dat <- fish_species_stream_year %>%
    filter(Species == sp, n_years_stream >= 5) %>%
    left_join(species_names, by = "Species")
  
  if (nrow(dat) == 0L || n_distinct(dat$Year) < 2L) next
  
  yrs <- range(dat$Year, na.rm = TRUE)
  sp_name  <- unique(dat$Common)
  title_txt <- if (length(sp_name)==0 || is.na(sp_name)) sp else paste0(sp_name, " (", sp, ")")
  
  p <- ggplot(dat, aes(Year, mean_adult_density, color = Stream, group = Stream)) +
    geom_line() +
    geom_point(size = 2) +
    geom_smooth(method = "lm", linetype = "dashed", se = FALSE, linewidth = 0.75) +
    scale_y_log10(breaks = log_breaks(10), labels = label_number(accuracy = 1, trim = TRUE)) +
    scale_x_continuous(breaks = seq(yrs[1], yrs[2], by = 10)) +  # no hard limits
    scale_color_viridis_d(option = "turbo") +
    labs(x = "Year", y = "Adult Density (mean per stream-year)", color = "Stream",
         title = paste0("Adult Density — ", title_txt)) +
    theme_plot(legend = TRUE)
  
  print(p)
  n_pages <- n_pages + 1L
}

# 3) Close device now (ensured by on.exit too)
dev.off()


