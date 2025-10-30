##############################################################
#  GRSM Integrated Cross-Taxon Dataset Assembly & Diagnostics
#  Author: John Grady
#  Date:   2025-10-16
##############################################################

#==============================================================
# 1. Libraries and global settings
#==============================================================

library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(stringr)

#--- Base path - upudate to your computer ---
general_path <- "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine"  

# Helper to shorten file.path calls
gp <- function(...) file.path(general_path, ...)

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
    plot.title        = element_text(size = 18, face = "plain", hjust = 1),  # NOTE: hjust=10 intentionally pushes title far right
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background  = element_blank(),
    strip.background  = element_blank(),
    legend.text       = element_text(size = 15),
    legend.title      = element_text(size = 18),
    legend.position   = if (legend) "right" else "none",   # toggle legend on/off
    text              = element_text(family = "Helvetica") # global font
  )
}
#==============================================================
# 2. Load & summarize individual datasets
#==============================================================

#--------------------------------------------------------------
# 2.1  FISH (Three-pass data)
#--------------------------------------------------------------
# Add location data
three_pass_data <- read_xlsx(
  gp("Data/Aquatics_Fish/Three_Pass/Summary_data/GRSM_Fish_3-Pass_Summary.xlsx"),
  sheet = "Summary") %>%
  rename(STATION_NAME = Site) %>%
  distinct(.) %>% # exclude duplicated rows
  left_join(
    read_csv(gp("Data/Aquatics_Fish/Three_Pass/Locations/GRSM_THREE_PASS.csv")) %>%
      select(LAT, LON, WATERSHED, STATION_NAME, STREAMNAME) %>%
      distinct(),
    by = "STATION_NAME") %>%
  rename(LOC_NAME = STATION_NAME)

# Write three pass data with location 
write_csv(three_pass_data, gp("Data/Aquatics_Fish/Three_Pass/Summary_data/GRSM_Fish_3-Pass_Summary_with_loc.csv"))

three_pass_summary <- three_pass_data %>%
  mutate(Year = year(Date)) %>%
  group_by(LOC_NAME, Year, WATERSHED, STREAMNAME, LAT, LON) %>%
  summarise(
    Richness = n_distinct(Species),
    ADTDens  = mean(ADTDens, na.rm = TRUE),
    ADTBiom  = mean(ADTBiom, na.rm = TRUE),
    .groups  = "drop"
  )

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


#--------------------------------------------------------------
# 2.2  MACROINVERTEBRATES
#--------------------------------------------------------------

# download and prepare location data for a join
invert_loc_codes <- read_csv(gp("Data/Aquatics_Macroinverts/Documents/Locations.csv")) %>%
  mutate(
    station_code = STATION_NAME %>%
      str_remove("(I&M|I\\&M)$") %>%
      str_replace_all("[^A-Z0-9]", "")) %>%
  select(station_code, STATION_NAME, LOC_NAME, Watershed, StreamName, LAT, LON) %>%
  distinct()

inverts_data <- read_csv(gp("Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export.csv")) %>%
  left_join(invert_loc_codes, by = "LOC_NAME")

# Response varaibles = Count (number of individuals per taxon)
# Categorical variables station code (sampling lcoation), LOC_NAME more descriptive but sometimes has multiple sampling codes

# Genus richness per site-year
genus_rich <- inverts_data %>%
  mutate(
    Year  = year(coalesce(Start_Date, as.Date(EventDate))),
    Genus = sub("\\s.*", "", Lab_Scientific_Name)
  ) %>%
  group_by(LOC_NAME, Year) %>%
  summarise(Genus_richness = n_distinct(Genus, na.rm = TRUE), .groups = "drop")

# Combine EPT & NCBI metrics
ept_clean <- read_csv(gp("Data/Aquatics_Macroinverts/NCBI_EPT/EPT_site_summary.csv")) %>%
  mutate(LOC_NAME = paste0(Location, ", ", Site)) %>%
  select(LOC_NAME, Year, EPT_rich)

ncbi_clean <- read_csv(gp("Data/Aquatics_Macroinverts/NCBI_EPT/NCBI_site_results.csv")) %>%
  mutate(LOC_NAME = paste0(Location, ", ", Site)) %>%
  select(LOC_NAME, Year, NCBI)

inv_metrics <- full_join(ept_clean, ncbi_clean, by = c("LOC_NAME", "Year"))

inverts_summary <- inv_metrics %>%
  left_join(genus_rich, by = c("LOC_NAME", "Year")) %>%
  left_join(
    invert_loc_codes %>% select(LOC_NAME, Watershed, StreamName, LAT, LON),
    by = "LOC_NAME"
  ) %>%
  select(LOC_NAME, Watershed, StreamName, LAT, LON, Year,
         EPT_rich, NCBI, Genus_richness)

#--------------------------------------------------------------
# 2.3  SOILS
#--------------------------------------------------------------
soils_data <- read_csv(gp("Data/Soil_Quality/Soils.csv")) %>%
  left_join(read_csv(gp("Data/Soil_Quality/soil_locations_all.csv")), by = "LOC_NAME") %>%
  mutate(Year = year(Event_Date))
#LOC_NAME is location name, response variables are other columns
#--------------------------------------------------------------
# 2.4  VEGETATION (trees + seedlings)
#--------------------------------------------------------------
forest_loc <- read_csv(gp("Data/Forest_Health/Locations.csv")) %>%
  select(LOC_NAME, VS_WATERSHED, LAT, LON)

trees_summary <- read_csv(gp("Data/Forest_Health/Trees.csv")) %>%
  left_join(forest_loc, by = "LOC_NAME") %>%
  mutate(Year = year(Event_Date)) %>%
  filter(is.na(Vigor_Desc) | Vigor_Desc != "Dead") %>%
  group_by(LOC_NAME, Year, VS_WATERSHED, LAT, LON) %>%
  summarise(Tree_richness = n_distinct(SpeciesCode),
            Tree_abundance = n(), .groups = "drop")

seedlings_summary <- read_csv(gp("Data/Forest_Health/Seedlings.csv")) %>%
  left_join(forest_loc, by = "LOC_NAME") %>%
  mutate(Year = year(Event_Date)) %>%
  group_by(LOC_NAME, Year, VS_WATERSHED, LAT, LON) %>%
  summarise(
    Seedling_richness  = n_distinct(SpeciesCode),
    Seedling_abundance = sum(Stem_Count * CORRECTION_FACTOR, na.rm = TRUE),
    .groups = "drop"
  )

plant_diversity <- full_join(
  trees_summary, seedlings_summary,
  by = c("LOC_NAME", "Year", "VS_WATERSHED", "LAT", "LON")
) %>%
  mutate(
    Total_richness  = rowSums(select(., Tree_richness, Seedling_richness), na.rm = TRUE),
    Total_abundance = rowSums(select(., Tree_abundance, Seedling_abundance), na.rm = TRUE)
  )
# Columns self explanatory

#--------------------------------------------------------------
# 2.5  CLIMATE (PRISM rasters)
#--------------------------------------------------------------
precip <- rast(gp("Data/Climate/PRISM_GRSM_precip_annual_mean_stack.tif"))
temp   <- rast(gp("Data/Climate/PRISM_GRSM_temp_annual_mean_stack.tif"))


# brown → green (for precipitation)
plot(precip$ppt_mean_1980,
     col = hcl.colors(20, "Terrain", rev = TRUE),   # earthy brown–green
     main = "Annual Precipitation 1980 (mm)")

# red → blue (for temperature)
plot(temp$tmean_mean_1980,
     col = rev(hcl.colors(20, "RdBu")),             # red–blue diverging
     main = "Annual Mean Temperature 1980 (°C)")

#==============================================================
# 3. Spatial integration & clustering (3.2 km)
#==============================================================
# -----------------------------------------------------------
# Coverage vs radius (justify 3.2 km) — place BEFORE clustering
# -----------------------------------------------------------
# Convert to sf, drop NA coords
fish_sf   <- three_pass_summary %>%
  filter(!is.na(LON), !is.na(LAT)) %>%
  st_as_sf(coords = c("LON","LAT"), crs = 4269)

invert_sf <- inverts_summary %>%
  filter(!is.na(LON), !is.na(LAT)) %>%
  st_as_sf(coords = c("LON","LAT"), crs = 4269)

plant_sf  <- plant_diversity %>%
  filter(!is.na(LON), !is.na(LAT)) %>%
  st_as_sf(coords = c("LON","LAT"), crs = 4269)

soil_sf   <- soils_data %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4269)

# Project each taxon to meters (UTM 17N)
fish_p  <- st_transform(fish_sf,   26917)
inv_p   <- st_transform(invert_sf, 26917)
plant_p <- st_transform(plant_sf,  26917)
soils_p <- st_transform(soil_sf,   26917)

# ---------- Fish baseline: nearest inverts / plants / soils ----------
idx_inv_for_fish   <- st_nearest_feature(fish_p,  inv_p)
idx_plant_for_fish <- st_nearest_feature(fish_p,  plant_p)
idx_soil_for_fish  <- st_nearest_feature(fish_p,  soils_p)

dist_invert_m_fish <- as.numeric(st_distance(fish_p,  inv_p[idx_inv_for_fish,  ], by_element = TRUE))
dist_plant_m_fish  <- as.numeric(st_distance(fish_p,  plant_p[idx_plant_for_fish,], by_element = TRUE))
dist_soil_m_fish   <- as.numeric(st_distance(fish_p,  soils_p[idx_soil_for_fish, ], by_element = TRUE))

# Minimal fish-baseline table for coverage calc
cross_taxon <- fish_sf %>%
  st_drop_geometry() %>%
  mutate(
    nearest_invert = inverts_summary$LOC_NAME[idx_inv_for_fish],
    nearest_plant  = plant_diversity$LOC_NAME[idx_plant_for_fish],
    nearest_soil   = soils_data$LOC_NAME[idx_soil_for_fish],
    dist_invert_m  = dist_invert_m_fish,
    dist_plant_m   = dist_plant_m_fish,
    dist_soil_m    = dist_soil_m_fish,
    # watershed checks (names differ across tables)
    same_ws_invert = (WATERSHED == inverts_summary$Watershed[idx_inv_for_fish]),
    same_ws_plant  = (WATERSHED == plant_diversity$VS_WATERSHED[idx_plant_for_fish])
  )

# ---------- Invert baseline: nearest fish / plants / soils ----------
idx_fish_for_inv   <- st_nearest_feature(inv_p,   fish_p)
idx_plant_for_inv  <- st_nearest_feature(inv_p,   plant_p)
idx_soil_for_inv   <- st_nearest_feature(inv_p,   soils_p)

dist_fish_m_inv    <- as.numeric(st_distance(inv_p, fish_p[idx_fish_for_inv,   ], by_element = TRUE))
dist_plant_m_inv   <- as.numeric(st_distance(inv_p, plant_p[idx_plant_for_inv, ], by_element = TRUE))
dist_soil_m_inv    <- as.numeric(st_distance(inv_p, soils_p[idx_soil_for_inv,  ], by_element = TRUE))

inv_nn <- invert_sf %>%
  st_drop_geometry() %>%
  mutate(
    nearest_fish   = three_pass_summary$LOC_NAME[idx_fish_for_inv],
    nearest_plant  = plant_diversity$LOC_NAME[idx_plant_for_inv],
    nearest_soil   = soils_data$LOC_NAME[idx_soil_for_inv],
    dist_fish_m    = dist_fish_m_inv,
    dist_plant_m   = dist_plant_m_inv,
    dist_soil_m    = dist_soil_m_inv,
    same_ws_fish   = (Watershed == three_pass_summary$WATERSHED[idx_fish_for_inv]),
    same_ws_plant  = (Watershed == plant_diversity$VS_WATERSHED[idx_plant_for_inv])
  )

# ---------- Coverage tables across candidate radii ----------
rads <- c(800, 1200, 1600, 2000, 3200, 5000)

coverage_fish_tbl <- purrr::map_dfr(rads, function(r){
  df <- cross_taxon %>%
    mutate(
      has_inverts = (dist_invert_m <= r) | (same_ws_invert %in% c(TRUE, "Y")),
      has_plants  = (dist_plant_m  <= r) | (same_ws_plant  %in% c(TRUE, "Y")),
      has_soils   =  dist_soil_m   <= r
    )
  tibble(
    baseline       = "fish",
    radius_m       = r,
    n_rows         = nrow(df),
    pct_has_inverts= mean(df$has_inverts, na.rm = TRUE),
    pct_has_plants = mean(df$has_plants,  na.rm = TRUE),
    pct_has_both   = mean(df$has_inverts & df$has_plants, na.rm = TRUE),
    pct_has_all    = mean(df$has_inverts & df$has_plants & df$has_soils, na.rm = TRUE),
    med_d_invert_m = median(df$dist_invert_m, na.rm = TRUE),
    med_d_plant_m  = median(df$dist_plant_m,  na.rm = TRUE),
    med_d_soil_m   = median(df$dist_soil_m,   na.rm = TRUE)
  )
})
coverage_fish_tbl

coverage_invert_tbl <- purrr::map_dfr(rads, function(r){
  df <- inv_nn %>%
    mutate(
      has_fish   = (dist_fish_m  <= r) | (same_ws_fish  %in% c(TRUE, "Y")),
      has_plants = (dist_plant_m <= r) | (same_ws_plant %in% c(TRUE, "Y")),
      has_soils  =  dist_soil_m  <= r
    )
  tibble(
    baseline     = "invert",
    radius_m     = r,
    n_rows       = nrow(df),
    pct_has_fish = mean(df$has_fish,   na.rm = TRUE),
    pct_has_plants = mean(df$has_plants, na.rm = TRUE),
    pct_has_both = mean(df$has_fish & df$has_plants, na.rm = TRUE),
    pct_has_all  = mean(df$has_fish & df$has_plants & df$has_soils, na.rm = TRUE),
    med_d_fish_m = median(df$dist_fish_m,  na.rm = TRUE),
    med_d_plant_m= median(df$dist_plant_m, na.rm = TRUE),
    med_d_soil_m = median(df$dist_soil_m,  na.rm = TRUE)
  )
})
coverage_invert_tbl

# quick comparison plot
coverage_both <- bind_rows(coverage_fish_tbl, coverage_invert_tbl)
ggplot(coverage_both, aes(radius_m, pct_has_both, color = baseline)) +
  geom_line() + geom_point(size = 3) +
  theme_plot() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Radius (m)", y = "% with both biotic partners",
       title = "Coverage vs radius from fish and invert baselines")

#------ Build clusters at 3200 radius

# Convert to sf, drop NA coords
fish_sf   <- three_pass_summary %>% filter(!is.na(LON), !is.na(LAT)) %>%
  st_as_sf(coords = c("LON","LAT"), crs = 4269)
invert_sf <- inverts_summary    %>% filter(!is.na(LON), !is.na(LAT)) %>%
  st_as_sf(coords = c("LON","LAT"), crs = 4269)
plant_sf  <- plant_diversity    %>% filter(!is.na(LON), !is.na(LAT)) %>%
  st_as_sf(coords = c("LON","LAT"), crs = 4269)
soil_sf   <- soils_data         %>% filter(!is.na(lon), !is.na(lat)) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4269)

# Combine taxa & project to meters (UTM 17N)
all_pts <- bind_rows(
  fish_sf %>% mutate(taxon="fish"),
  invert_sf %>% mutate(taxon="invert"),
  plant_sf %>% mutate(taxon="plant"),
  soil_sf %>% mutate(taxon="soil")
)
all_m <- st_transform(all_pts, 26917)

# Build clusters within 3.2 km
r_m <- 3200
within_ix <- st_is_within_distance(all_m, all_m, dist = r_m)
cluster_id <- rep(NA_integer_, length(within_ix))
current <- 0L
unassigned <- seq_along(within_ix)

while(length(unassigned)){
  current <- current + 1L
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

all_with_cluster <- all_pts %>% mutate(cluster_id = cluster_id)

# Cluster centroids
cluster_centroids <- all_with_cluster %>%
  group_by(cluster_id) %>%
  summarise(geometry = st_centroid(st_union(geometry))) %>%
  st_as_sf()

#==============================================================
# 4. Aggregate to cluster-year
#==============================================================

# Lookups (LOC_NAME → cluster_id)
lookup_tbl <- function(df, tax) df %>% filter(taxon==tax) %>%
  st_drop_geometry() %>% distinct(LOC_NAME, cluster_id)
fish_lookup  <- lookup_tbl(all_with_cluster,"fish")
inv_lookup   <- lookup_tbl(all_with_cluster,"invert")
plant_lookup <- lookup_tbl(all_with_cluster,"plant")
soil_lookup  <- lookup_tbl(all_with_cluster,"soil")

# Prep per-taxon summaries joined to clusters
fish_dat <- three_pass_summary %>%
  left_join(fish_lookup, by="LOC_NAME") %>%
  rename(Fish_Richness=Richness, Fish_ADTDens=ADTDens, Fish_ADTBiom=ADTBiom) %>%
  select(cluster_id, Year, Fish_Richness, Fish_ADTDens, Fish_ADTBiom)

inv_dat <- inverts_summary %>%
  left_join(inv_lookup, by="LOC_NAME") %>%
  rename(Invert_EPT_rich=EPT_rich, Invert_NCBI=NCBI, Invert_Genus_richness=Genus_richness) %>%
  select(cluster_id, Year, Invert_EPT_rich, Invert_NCBI, Invert_Genus_richness)

plant_dat <- plant_diversity %>%
  left_join(plant_lookup, by="LOC_NAME") %>%
  rename(Plant_Tree_richness=Tree_richness, Plant_Tree_abundance=Tree_abundance,
         Plant_Seedling_richness=Seedling_richness, Plant_Seedling_abundance=Seedling_abundance) %>%
  select(cluster_id, Year, Plant_Tree_richness, Plant_Tree_abundance,
         Plant_Seedling_richness, Plant_Seedling_abundance)

soils_site_year <- soils_data %>%
  mutate(Year=year(Event_Date)) %>%
  group_by(LOC_NAME,Year) %>%
  summarise(across(c(SoilpH,CEC,pOM,Ca_ppm,Mg_ppm,K_ppm), mean, na.rm=TRUE),
            .groups="drop") %>%
  rename(Soil_pH=SoilpH, Soil_CEC=CEC, Soil_pOM=pOM) %>%
  left_join(soil_lookup, by="LOC_NAME") %>%
  select(cluster_id,Year,Soil_pH,Soil_CEC,Soil_pOM,Soil_Ca_ppm=Ca_ppm,
         Soil_Mg_ppm=Mg_ppm,Soil_K_ppm=K_ppm)

# Aggregate by cluster-year
agg_mean <- function(df) summarise(df, across(where(is.numeric), mean, na.rm=TRUE), .groups="drop")
agg_sum  <- function(df) summarise(df, across(where(is.numeric), sum, na.rm=TRUE), .groups="drop")

cluster_fish  <- fish_dat  %>% group_by(cluster_id,Year) %>% agg_mean()
cluster_inv   <- inv_dat   %>% group_by(cluster_id,Year) %>% agg_mean()
cluster_plant <- plant_dat %>% group_by(cluster_id,Year) %>%
  summarise(across(ends_with("richness"), mean, na.rm=TRUE),
            across(ends_with("abundance"), sum, na.rm=TRUE), .groups="drop")
cluster_soil  <- soils_site_year %>% group_by(cluster_id,Year) %>% agg_mean()

cluster_year <- reduce(list(cluster_fish, cluster_inv, cluster_plant, cluster_soil),
                       full_join, by=c("cluster_id","Year"))

#==============================================================
# 5. Attach climate to cluster centroids
#==============================================================
clusters_ll <- st_transform(cluster_centroids, 4269)
clusters_v  <- vect(clusters_ll)
ppt_years   <- as.integer(str_extract(names(precip),"\\d{4}"))
tmean_years <- as.integer(str_extract(names(temp),"\\d{4}"))
yrs <- sort(unique(cluster_year$Year))
out_list <- list()

for(y in yrs){
  i_ppt <- which(ppt_years==y)
  i_tmn <- which(tmean_years==y)
  if(!length(i_ppt) || !length(i_tmn)) next
  vals_ppt <- terra::extract(precip[[i_ppt[1]]], clusters_v)
  vals_tmn <- terra::extract(temp[[i_tmn[1]]], clusters_v)
  out_list[[length(out_list)+1]] <-
    tibble(cluster_id=clusters_ll$cluster_id, Year=y,
           Precip_mm=as.numeric(vals_ppt[[2]]),
           Temp_C=as.numeric(vals_tmn[[2]]))
}
clusters_climate <- bind_rows(out_list)
cluster_year <- left_join(cluster_year, clusters_climate, by=c("cluster_id","Year"))

#==============================================================
# 6. Diagnostics
#==============================================================
cluster_basic <- tibble(
  n_points=nrow(all_with_cluster),
  n_clusters=n_distinct(all_with_cluster$cluster_id),
  n_missing=sum(is.na(all_with_cluster$cluster_id))
)
print(cluster_basic)

#==============================================================
# 7. Visualization functions
#==============================================================

vars_present <- function(df, vars){
  keep <- vars[vars %in% names(df)]
  if(!length(keep)) stop("None of the requested variables found.")
  keep[vapply(keep,function(v) any(!is.na(df[[v]])),logical(1))]
}

plot_yearly_trend <- function(df, vars, title=NULL){
  keep <- vars_present(df,vars)
  df_long <- df %>%
    select(cluster_id,Year,all_of(keep)) %>%
    pivot_longer(all_of(keep),names_to="variable",values_to="value") %>%
    group_by(Year,variable) %>%
    summarise(med=median(value,na.rm=TRUE),
              p25=quantile(value,0.25,na.rm=TRUE),
              p75=quantile(value,0.75,na.rm=TRUE),.groups="drop")
  ggplot(df_long,aes(Year,med))+
    geom_ribbon(aes(ymin=p25,ymax=p75),alpha=.25)+
    geom_line(linewidth=1)+
    facet_wrap(~variable,scales="free_y",ncol=3)+
    labs(title=title,subtitle="Median ± IQR")
}

plot_spaghetti <- function(df, vars, title=NULL,sample_n=NULL,seed=1){
  keep <- vars_present(df,vars)
  df_use <- df %>% select(cluster_id,Year,all_of(keep))
  if(!is.null(sample_n)){
    set.seed(seed)
    ids <- df_use %>% distinct(cluster_id) %>% slice_sample(n=sample_n) %>% pull()
    df_use <- df_use %>% filter(cluster_id %in% ids)
  }
  df_long <- df_use %>%
    pivot_longer(all_of(keep),names_to="variable",values_to="value")
  med_year <- df_long %>%
    group_by(Year,variable) %>%
    summarise(med=median(value,na.rm=TRUE),.groups="drop")
  ggplot(df_long,aes(Year,value,group=cluster_id))+
    geom_line(alpha=.25)+
    theme_plot() +
    geom_line(data=med_year,aes(Year,med,group=variable),
              linewidth=1.2,inherit.aes=FALSE)+
    scale_y_log10()+
    facet_wrap(~variable,scales="free_y",ncol=3)+
    labs(title=title,subtitle="Thin = clusters; Thick = median")
}

plot_z_heatmap <- function(df, vars, title=NULL){
  keep <- vars_present(df,vars)
  df_long <- df %>% select(cluster_id,Year,all_of(keep)) %>%
    pivot_longer(all_of(keep),names_to="variable",values_to="value") %>%
    group_by(variable)%>%
    mutate(z=(value-mean(value,na.rm=TRUE))/sd(value,na.rm=TRUE))%>%
    ungroup()
  ref_year <- max(df_long$Year,na.rm=TRUE)
  clust_order <- df_long%>%filter(Year==ref_year)%>%
    group_by(cluster_id)%>%summarise(mz=mean(z,na.rm=TRUE),.groups="drop")%>%
    arrange(desc(mz))%>%pull(cluster_id)
  df_long$cluster_id <- factor(df_long$cluster_id,levels=clust_order)
  ggplot(df_long,aes(Year,cluster_id,fill=z))+
    geom_tile()+facet_wrap(~variable,scales="free",ncol=3)+
    labs(title=title,subtitle=paste("Z-scores; ref year =",ref_year))
}


#==============================================================
# 8. Correlations among variables
#==============================================================

var_groups <- list(
  fish=c("Fish_Richness","Fish_ADTDens","Fish_ADTBiom"),
  invert=c("Invert_EPT_rich","Invert_NCBI","Invert_Genus_richness"),
  plant=c("Plant_Tree_richness","Plant_Seedling_richness",
          "Plant_Tree_abundance","Plant_Seedling_abundance"),
  soil=c("Soil_pH","Soil_CEC","Soil_pOM","Soil_Ca_ppm","Soil_Mg_ppm","Soil_K_ppm"),
  clim=c("Precip_mm","Temp_C")
)

vars_for_corr <- unlist(var_groups)
df_corr <- cluster_year %>% select(any_of(vars_for_corr)) %>% drop_na()
corr_mat <- cor(df_corr,use="pairwise.complete.obs",method="spearman")

as.data.frame(as.table(corr_mat)) %>%
  rename(var1=Var1,var2=Var2,rho=Freq) %>%
  ggplot(aes(var1,var2,fill=rho))+
  geom_tile()+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,limits=c(-1,1))+
  coord_fixed()+
  theme_minimal(base_size=10)+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  labs(title="Spearman correlations among cluster-year variables",fill="ρ")



##############################################################
# plots - generic plot to show integration of data worked 
##############################################################
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
