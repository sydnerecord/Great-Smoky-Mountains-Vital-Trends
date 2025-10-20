# ======================================================================
# Libraries
# ======================================================================
library(RColorBrewer)  # color palettes 
library(readxl)        # read_xlsx() 
library(viridisLite)   # color options
library(tidyverse)     # data import (readr), manipulation (dplyr), plotting (ggplot2), strings (stringr)
library(sf)            # spatial data

# ======================================================================
# Plot theme
# - Compact, high-contrast theme for time-series panels
# - 'legend' arg toggles legend visibility without re-specifying theme
# ======================================================================
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

# ======================================================================
# Input data: GRSM macroinvertebrate specimen export
# ===============================================================
inverts = read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export_locations.csv')
identical(inverts$station_code, inverts$station_code2)
ncbi = read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/NCBI_EPT/NCBI_results.csv')
ept = read_csv("~/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/NCBI_EPT/EPT_summary.csv")


# ----------------------------------------------------------------------
# Whitespace cleanup
# - Trim only character columns; leave numeric/integer columns unchanged
# ----------------------------------------------------------------------
inverts <- inverts %>%
  mutate(across(where(is.character), str_squish))

# ----------------------------------------------------------------------
# Helper parsing columns
# - Location/Site parsed from LOC_NAME = "Location, Site, ..."
# - Genus parsed from first token of Lab_Scientific_Name
# - Year parsed from first token of Start_Date (assumes "YYYY-..." or "YYYY")
# ----------------------------------------------------------------------
inverts$Location <- stringr::word(inverts$LOC_NAME, 1, sep = ",")
inverts$Site     <- stringr::word(inverts$LOC_NAME, 2, sep = ",")
inverts$Genus    <- stringr::word(inverts$Lab_Scientific_Name, 1, sep = " ")
inverts$Year     <- as.numeric(stringr::word(inverts$Start_Date, 1, sep = "-"))

# ----------------------------------------------------------------------
# ---------------- Plotted Trends ---------------------------------------
# ----------------------------------------------------------------------

#---- Inverts abundance and richness of general across sites --
inverts_site_year <- inverts %>%
  group_by(Location, Site, Year) %>%
  summarise(
    n_genera = n_distinct(Genus, na.rm = TRUE),
    total_abundance = sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

#---- Averged sites per stream --
inverts_avg_stream <- inverts_site_year %>%
  group_by(Location, Year) %>%
  summarise(
    mean_richness   = mean(n_genera, na.rm = TRUE),          # average genera per site
    mean_abundance  = mean(total_abundance, na.rm = TRUE),   # average abundance per site
    n_sites         = n(),                                   # number of sites contributing
    .groups = "drop"
  )

# Add years sampled
inverts_avg_loc <- inverts_avg_stream  %>%
  dplyr::group_by(Location) %>%
  dplyr::mutate(n_years_sampled = dplyr::n_distinct(Year)) %>%
  dplyr::ungroup()

invert_summary = inverts_avg_loc 

invert_summary <- invert_summary %>%
  left_join(ept %>% select(Location, Year, EPT_rich), by = c("Location", "Year")) %>%
  left_join(ncbi %>% select(Location, Year, NCBI_mean), by = c("Location", "Year"))
inverts
invert_summary$log_richness = log(invert_summary$mean_richness)
invert_summary$log_abundance = log(invert_summary$mean_abundance)
invert_summary$log_richness [!is.finite(invert_summary$log_richness)]  <- NA
invert_summary$log_abundance[!is.finite(invert_summary$log_abundance)] <- NA
# ------ Plots ------

# --- Richness
p_richness <- ggplot(
  invert_summary %>% filter(n_years_sampled >= 5, mean_richness > 0),
  aes(x = Year, y = mean_richness, color = Location, group = Location)) +
  geom_line() +
  geom_point(size = 2) +
  geom_smooth(method = "lm", linetype = "dashed", se = F, size = .75) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = seq(1995, 2020, by = 5),
    limits = c(1995, 2020)) +
  labs(x = "Year", y = "Mean Genus Richness per Stream", color = "Stream") +
  theme_plot() +
  ggtitle("Genera Richness") +
  scale_color_viridis_d(option = "turbo") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right") +
  annotate("segment", x = 2017.2, xend = 2020, y = 112, yend = 112,
           color = "black", linetype = "dashed") +
  annotate("text", x = 2013.5, y = 110, label = "Linear Fit", size = 5,
           hjust = 0, vjust = 0, color = "black")
p_richness


# --- Abundance
p_abundance <- ggplot(
  invert_summary %>% filter(n_years_sampled >= 5, mean_richness > 0),
  aes(x = Year, y = mean_abundance, color = Location, group = Location)) +
  geom_line() +
  geom_point(size = 2) +
  geom_smooth(method = "lm", linetype = "dashed", se = F, size = .75) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = seq(1985, 2020, by = 5),
    limits = c(1984, 2020)) +
  labs(x = "Year", y = "Genera Abundance per Stream", color = "Stream") +
  theme_plot() +
  ggtitle("Mean Genera Abundance") +
  scale_color_viridis_d(option = "turbo") +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right") +
  annotate("segment", x = 2017, xend = 2019.5, y = 3100, yend = 3100,
           color = "black", linetype = "dashed") +
  annotate("text", x = 2011.5, y = 3000, label = "Linear Fit", size = 5,
           hjust = 0, vjust = 0, color = "black")
p_abundance

# --- NCBI
ref_lines <- data.frame(
  y     = c(4.18, 5.09, 5.91, 7.05),
  label = c("Excellent", "Good", "Good–Fair", "Poor")) 
x_left <- min(invert_summary$Year, na.rm = TRUE)

p_ncbi <- ggplot(
  invert_summary %>% filter(n_years_sampled >= 5, NCBI_mean > 0),
  aes(x = Year, y = NCBI_mean, color = Location, group = Location)) +
  geom_line() +
  geom_point(size = 2) +
  geom_smooth(method = "lm", linetype = "dashed", se = F, size = .75) +
  scale_x_continuous(
    breaks = seq(1985, 2020, by = 5),
    limits = c(1984, 2020)) +
  labs(x = "Year", y = "NCBI", color = "Stream") +
  scale_y_continuous(limits = c(min(ncbi$NCBI_mean, na.rm = TRUE), pmax(max(ncbi$NCBI_mean, na.rm = TRUE), 4.4))) +
  geom_text(
    data = ref_lines, aes(x = x_left, y = y + 0.1, label = label),
    color  = "grey30", size = 4, hjust = 0, inherit.aes = FALSE) +
  theme_plot() +
  scale_color_viridis_d(option = "turbo") +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right") +
  ggtitle("Stream NCBI") +
  annotate("segment", x = 2017, xend = 2019.5, y = 3.9, yend = 3.9,
           color = "black", linetype = "dashed") +
  annotate("text", x = 2011.5, y = 3.85, label = "Linear Fit", size = 5,
           hjust = 0, vjust = 0, color = "black") +
  geom_hline(data = ref_lines, aes(yintercept = y), color = "grey40", 
             linetype = "dashed", linewidth = 0.8,inherit.aes = FALSE) 
p_ncbi

# EPT

p_ept <- ggplot(
  invert_summary %>% filter(n_years_sampled >= 5, EPT_rich> 0),
  aes(x = Year, y = EPT_rich, color = Location, group = Location)) +
  geom_line() +
  geom_point(size = 2) +
  geom_smooth(method = "lm", linetype = "dashed", se = F, size = .75) +
  scale_x_continuous(
    breaks = seq(1985, 2020, by = 5),
    limits = c(1984, 2020)) +
  labs(x = "Year", y = "EPT", color = "Stream") +
  theme_plot() +
  ggtitle("EPT") +
  scale_color_viridis_d(option = "turbo") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right") +
  annotate("segment", x = 2017, xend = 2019.5, y = 60, yend = 60,
           color = "black", linetype = "dashed") +
  annotate("text", x = 2011, y = 59, label = "Linear Fit", size = 5,
           hjust = 0, vjust = 0, color = "black")
p_ept


#------ Attach geometry for plotting

coords_sf <- inverts %>%
  st_as_sf(coords = c("LON","LAT"), crs = 4326, remove = FALSE) %>%
  mutate(stream = str_squish(str_remove(LOC_NAME, ",?\\s*Site\\s*\\w+$"))) %>%
  group_by(stream) %>%
  summarise(geometry = sf::st_centroid(sf::st_union(geometry)), .groups = "drop")

invert_summary_sf <- invert_summary %>%
  left_join(coords_sf %>% dplyr::rename(Location = stream), by = "Location") %>%
  sf::st_as_sf()

# https://grsm-nps.opendata.arcgis.com/datasets/60eb34ba0a354554ada335a11b83180f_0/explore?location=35.641900%2C-83.544595%2C9.86
#https://grsm-nps.opendata.arcgis.com/datasets/60eb34ba0a354554ada335a11b83180f_0/explore?location=35.641900%2C-83.544595%2C9.86

# Read watershed polygons (index [2] keeps the desired layer per your workflow)
watershed <- st_read('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Locations/GRSM_WATERSHEDS/GRSM_WATERSHEDS.shp')[2]
plot(watershed)  # quick base plot check

grsm_border <- st_read('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Locations/BOUNDARY_LN/BOUNDARY_LN.shp')[6] %>%
  st_transform(st_crs(watershed))
str(grsm_border)

grsm_border$GlobalID


# 1) Make everything the same CRS as watershed

invert_summary_sf <- st_transform(invert_summary_sf, st_crs(watershed))

# 2) Keep only streams that appear in invert_summary_sf$Location (name match)
streams <- st_read('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Locations/Streams/GRSM_HYDROLOGY.geojson')[3]

# Minimal prep for colorizing by name
streams_named <- streams %>%
  filter(!is.na(GNIS_NAME), GNIS_NAME != "") %>%
  mutate(color_id = factor(GNIS_NAME))
streams_used <- streams_named %>%
  mutate(name_norm = stringr::str_squish(stringr::str_to_lower(GNIS_NAME)))

loc_used <- invert_summary_sf %>%
  st_drop_geometry() %>%
  dplyr::distinct(Location) %>%
  dplyr::mutate(name_norm = stringr::str_squish(stringr::str_to_lower(Location)))

streams_used <- streams_used %>%
  dplyr::semi_join(loc_used, by = "name_norm")


# join locations to coords; keep one point per Location
sites_pts <- coords_sf %>%
  inner_join(loc_used, by = c("stream" = "Location")) %>%  # keep only sampled sites
  mutate(Location = stream) %>%
  dplyr::select(Location, geometry) %>%
  sf::st_transform(sf::st_crs(watershed))

# one label point per stream name (keep all stream lines)
stream_labels <- streams_used %>%
  dplyr::group_by(GNIS_NAME, color_id) %>%
  dplyr::summarise(do_union = TRUE, .groups = "drop") %>%
  sf::st_transform(26917) %>%
  dplyr::mutate(geometry = sf::st_point_on_surface(geometry)) %>%  # one point near the middle
  sf::st_transform(sf::st_crs(streams_used))


# 4) Plot: watersheds + border + filtered streams + locations
ggplot() +
  geom_sf(data = grsm_border,  fill = NA, color = "gray50", linewidth = 1) +
  geom_sf(data = watershed,    fill = NA, color = "grey30", linewidth = 0.15) +
  geom_sf(data = streams_used, aes(color = GNIS_NAME), linewidth = .75, alpha = 0.9) +
  geom_sf_text(data = stream_labels, aes(label = GNIS_NAME, color = GNIS_NAME),
               size = 3, check_overlap = TRUE) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  geom_sf(data = sites_pts, aes(color = Location), size = 2) +
  coord_sf() +
  labs(title = "GRSM Invert Locations") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank()) 


ggplot() +
  geom_sf(data = grsm_border,  fill = NA, color = "gray50", linewidth = 1) +
  #geom_sf(data = watershed,    fill = NA, color = "grey30", linewidth = 0.15) +
  geom_sf(data = streams_used, aes(color = GNIS_NAME), linewidth = .75, alpha = 0.9) +
  geom_sf_text(data = stream_labels, aes(label = GNIS_NAME, color = GNIS_NAME),
               size = 4, check_overlap = TRUE) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  geom_sf(data = sites_pts, aes(color = Location), size = 3) +
  coord_sf() +
  labs(title = "GRSM Invert Locations") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())



#------ Get Slopes -----


responses <- c("log_richness", "log_abundance", "EPT_rich", "NCBI_mean")

invert_summary$Year <- as.numeric(invert_summary$Year)


fit_stream_trend <- function(stream_data, response_var){
  m  <- lm(stats::reformulate("Year", response_var), data=stream_data)
  cs <- coef(summary(m))
  i  <- match("Year", rownames(cs))
  tibble(slope = ifelse(is.na(i), NA_real_, cs[i,"Estimate"]),
         p_value = ifelse(is.na(i), NA_real_, cs[i,"Pr(>|t|)"]),
         r_sq = summary(m)$r.squared)
}

responses <- c("log_richness","log_abundance","EPT_rich","NCBI_mean")

lm_per_stream <- purrr::map_dfr(responses, function(var){
  invert_summary %>%
    filter(n_years_sampled>=5, is.finite(.data[[var]])) %>%
    group_by(Location) %>%
    group_modify(~ fit_stream_trend(.x, var)) %>%
    mutate(response=var, significant=ifelse(!is.na(p_value) & p_value<0.05,"yes","no")) %>%
    ungroup()
})

loc_geom <- invert_summary_sf %>% group_by(Location) %>% slice_head(n=1) %>% select(Location, geometry)
lm_per_stream_sf <- lm_per_stream %>% left_join(loc_geom, by="Location") %>% sf::st_as_sf()



# Richness
pts <- lm_per_stream_sf %>% filter(response == "log_richness")
max_abs <- max(abs(pts$slope), na.rm = TRUE)
gray_out <- F  # TRUE = gray nonsignificant, FALSE = color all
pts <- pts %>% mutate(fill_val = if (gray_out) ifelse(significant == "yes", slope, NA_real_) else slope)

# one map per species (prints each to the plotting device)
for (sp in sort(unique(lm_species_stream_sf$Species))) {
  dat <- lm_species_stream_sf %>% dplyr::filter(Species == sp)
  
  if (nrow(dat) == 0) next
  
  max_abs <- max(abs(dat$slope), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  p <- ggplot() +
    geom_sf(data = grsm_border %>% sf::st_transform(sf::st_crs(watershed)),
            fill = NA, color = "gray50", linewidth = 1) +
    geom_sf(data = watershed, fill = NA, color = "grey30", linewidth = 0.15) +
    geom_sf(data = streams_matched_sf, color = "grey70", linewidth = 0.5, alpha = 0.7) +
    geom_sf(data = dat, aes(fill = slope),
            shape = 21, color = "black", size = 4, stroke = 0.3) +
    scale_fill_gradientn(colors = c("#1E90FF", "white", "red"),
                         limits = c(-max_abs, max_abs),
                         na.value = "grey80", name = "Slope") +
    coord_sf() +
    labs(title = paste("Adult Density Trend (slope) —", sp)) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank())
  
  print(p)
}



# Abundance
gray_out <- F # toggle: TRUE = gray nonsignificant, FALSE = color all

pts <- lm_per_stream_sf %>% filter(response == "log_abundance")
max_abs <- max(abs(pts$slope), na.rm = TRUE)
pts <- pts %>% mutate(fill_val = if (gray_out) ifelse(significant == "yes", slope, NA_real_) else slope)

ggplot() +
  geom_sf(data = grsm_border, fill = NA, color = "gray50", linewidth = 1) +
  geom_sf(data = streams_used, color = "grey70", linewidth = 0.5, alpha = 0.7) +
  geom_sf_text(data = stream_labels, aes(label = GNIS_NAME), color = "black",
               size = 3, check_overlap = TRUE, show.legend = F) +
  geom_sf(data = pts, aes(fill = fill_val),
          shape = 21, color = "black", size = 7, stroke = 0.3) +
  scale_fill_gradientn(
    colors = c("#1E90FF", "white", "red"),
    limits = c(-max_abs, max_abs),
    na.value = "grey70") +
  coord_sf() +
  labs(title = "Invertebrate Abundance Trends", fill = "Slope") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())


# EPT richness
pts <- lm_per_stream_sf %>% filter(response == "EPT_rich")
max_abs <- max(abs(pts$slope), na.rm = TRUE)
pts <- pts %>% mutate(fill_val = if (gray_out) ifelse(significant == "yes", slope, NA_real_) else slope)
gray_out <- F

ggplot() +
  geom_sf(data = grsm_border, fill = NA, color = "gray50", linewidth = 1) +
  geom_sf(data = streams_used, color = "grey70", linewidth = 0.5, alpha = 0.7) +
  geom_sf_text(data = stream_labels, aes(label = GNIS_NAME), size = 3, color = "grey10", check_overlap = TRUE) +
  geom_sf(data = pts, aes(fill = fill_val), shape = 21, color = "black", size = 7, stroke = 0.3) +
  scale_fill_gradientn(colors = c("#1E90FF", "white", "red"),
                       limits = c(-max_abs, max_abs),
                       na.value = "grey70") +
  coord_sf() +
  labs(title = "EPT Trends", fill = "Slope") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())


# NCBI

gray_out <- F  # TRUE = gray nonsignificant, FALSE = color all
pts <- lm_per_stream_sf %>% filter(response == "NCBI_mean")
max_abs <- max(abs(pts$slope), na.rm = TRUE)
pts <- pts %>% mutate(fill_val = if (gray_out) ifelse(significant == "yes", slope, NA_real_) else slope)

ggplot() +
  geom_sf(data = grsm_border, fill = NA, color = "gray50", linewidth = 1) +
  geom_sf(data = streams_used, color = "grey70", linewidth = 0.5, alpha = 0.7) +
  geom_sf_text(data = stream_labels, aes(label = GNIS_NAME),
               size = 3, color = "grey10", check_overlap = TRUE) +
  geom_sf(data = pts, aes(fill = fill_val),
          shape = 21, color = "black", size = 7, stroke = 0.3) +
  scale_fill_gradientn(colors = c("#1E90FF", "white", "red"), limits = c(-max_abs, max_abs),na.value = "grey70") +
  coord_sf() +
  labs(title = "NCBI Trends", fill = "Slope") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

