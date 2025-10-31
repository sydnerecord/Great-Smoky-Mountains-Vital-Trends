

library(tidyverse)
library(ggforce)
library(forcats)


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
#--- Base path - update to your computer ---
general_path <- "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine"  

# Helper to shorten file.path calls
gp <- function(...) file.path(general_path, ...)


# Plants
#--------------------------------------------------------------

#LOC_NAME is location name, response variables are other columns
#--------------------------------------------------------------
# 2.4  VEGETATION (trees + seedlings)
#--------------------------------------------------------------

trees <- read_csv(gp("Data/Forest_Health/Trees_with_coordinates.csv")) %>%
  mutate(Year = as.integer(lubridate::year(Event_Date)))

seedlings <- read_csv(gp("Data/Forest_Health/Seedlings_with_coordinates.csv")) %>%
  mutate(Year = as.integer(lubridate::year(Event_Date)))

woody <- read_csv(gp("Data/Forest_Health/Woody_Stems_with_coordinates.csv")) %>%
  mutate(Year = as.integer(lubridate::year(Event_Date)))

canopy = st_read(gp('Data/Forest_Health/veg_comm_struct_with_coordinates.csv'))

# Years sampled
unique(seedlings$Year)
unique(trees$Year)
unique(woody$Year)
#[1] 2015 2023 2017 2019 2020 2024 2018 2016 2021


#--------- Calculate richness and abundance -----
trees_rich_abun <- trees %>%
  group_by(LOC_NAME, Year, VS_WATERSHED, LAT, LON) %>%
  summarise(
    Tree_richness  = n_distinct(SpeciesCode),
    Tree_abundance = n(),
    .groups = "drop"
  )

seedlings_rich_abun <- seedlings %>%
  group_by(LOC_NAME, Year, VS_WATERSHED, LAT, LON) %>%
  summarise(
    Seedling_richness  = n_distinct(SpeciesCode),
    Seedling_abundance = sum(Stem_Count * CORRECTION_FACTOR, na.rm = TRUE),
    .groups = "drop"
  )


# add a per-record Stem_Count from the size-class columns, then aggregate
woody <- woody %>%
  mutate(
    CORRECTION_FACTOR = tidyr::replace_na(CORRECTION_FACTOR, 1), #if NA replace with 1 so multiplying yields itself
    Stem_Count = rowSums(dplyr::across(starts_with("Count"), ~ tidyr::replace_na(., 0))) # For each row (species × site × year), this sums across all the size-class count columns:
  )

woody_rich_abun <- woody %>%
  dplyr::group_by(LOC_NAME, Year, VS_WATERSHED, LAT, LON) %>%
  dplyr::summarise(
    woody_richness  = dplyr::n_distinct(SpeciesCode),
    woody_abundance = sum(Stem_Count * CORRECTION_FACTOR, na.rm = TRUE),
    .groups = "drop"
  )

plant_diversity <- full_join(
  trees_rich_abun, seedlings_rich_abun, woody_rich_abun,
  by = c("LOC_NAME", "Year", "VS_WATERSHED", "LAT", "LON")
) %>%
  mutate(
    Total_richness  = rowSums(dplyr::select(., Tree_richness, Seedling_richness), na.rm = TRUE),
    Total_abundance = rowSums(dplyr::select(., Tree_abundance, Seedling_abundance), na.rm = TRUE)
  )



#-----------------------------------------------
# Presence/absence for seedlings across sites
#-----------------------------------------------
# Which sites have more than one year sampled? 


year_summary <- seedlings %>%
  group_by(LOC_NAME) %>%
  summarise(
    n_years = n_distinct(Year),
    years_sampled = paste(sort(unique(Year)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_years))

# View the top sites
print(year_summary, n = 20)

# Sites with more than one year of sampling
multi_year_sites <- year_summary %>%
  filter(n_years > 1)

multi_year_sites

#------------------------------------------
#------- Plot Seedling presence absence matrix
#-------------------------------------------
# Presence–absence heatmaps for ALL sites, 9 sites per page (3x3 facets).
# Colors: green = present, red = absent (in sampled years), gray = unsampled year.


plot_all_sites_pa <- function(df, growth_form = "Tree", pdf_out = NULL) {
  df <- df %>%
    mutate(Year = as.integer(Year)) %>%
    filter(GrowthForm == growth_form)
  
  # global year axis (so unsampled years can be grayed)
  all_years <- sort(unique(df$Year))
  
  # species seen at each site (limits rows per facet to locally observed spp)
  species_by_site <- df %>%
    distinct(LOC_NAME, SciName)
  
  # site-years that were sampled (for this growth form)
  site_years <- df %>%
    distinct(LOC_NAME, Year) %>%
    mutate(sampleed = TRUE) %>%
    rename(sampled = sampleed)
  
  # species x year presence within site
  pres <- df %>%
    group_by(LOC_NAME, SciName, Year) %>%
    summarise(present = as.integer(sum(Stem_Count, na.rm = TRUE) > 0), .groups = "drop")
  
  # full grid: (LOC_NAME, SciName seen at that site) x all years
  pa_full <- species_by_site %>%
    crossing(Year = all_years) %>%
    left_join(pres, by = c("LOC_NAME", "SciName", "Year")) %>%
    left_join(site_years, by = c("LOC_NAME", "Year")) %>%
    mutate(
      fill_color = case_when(
        is.na(sampled)                ~ "unsampled",         # no sampling that year at this site
        !is.na(sampled) & present == 1 ~ "present",          # sampled & seen
        !is.na(sampled) & (is.na(present) | present == 0) ~ "absent" # sampled & not seen
      )
    )
  
  base_plot <- ggplot(pa_full, aes(x = factor(Year),
                                   y = fct_rev(fct_inorder(SciName)),
                                   fill = fill_color)) +
    geom_tile(color = "grey60", linewidth = 0.2) +
    scale_fill_manual(
      values = c(present = "darkgreen", absent = "red3", unsampled = "grey80"),
      breaks = c("present", "absent", "unsampled"),
      name = "Status",
      labels = c("Present", "Absent", "Not sampled")
    ) +
    labs(x = "Year", y = "Species",
         title = "Tree presence–absence by site (9 sites per page)") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "top",
      axis.text.y = element_text(size = 10)
    )
  
  # how many pages needed at 3x3 facets
  n_pages <- ggforce::n_pages(
    base_plot + ggforce::facet_wrap_paginate(~ LOC_NAME, ncol = 3, nrow = 3, page = 1, scales = "free_y")
  )
  
  if (!is.null(pdf_out)) {
    pdf(file = pdf_out, width = 11, height = 8.5)
    for (i in seq_len(n_pages)) {
      print(base_plot + ggforce::facet_wrap_paginate(~ LOC_NAME, ncol = 3, nrow = 3, page = i, scales = "free_y"))
    }
    dev.off()
  } else {
    # print all pages to the display device
    for (i in seq_len(n_pages)) {
      print(base_plot + ggforce::facet_wrap_paginate(~ LOC_NAME, ncol = 3, nrow = 3, page = i, scales = "free_y"))
    }
  }
}

# Example:
 plot_all_sites_pa(seedlings, growth_form = "Tree")

 # or save to PDF (set your own path explicitly):
plot_all_sites_pa(seedlings, pdf_out = "~/Downloads/tree_presence_absence_by_site.pdf")

#------------------------------------------
#------- Plot TREE presence–absence matrix
#------------------------------------------
# Presence = any record for (SciName, LOC_NAME) in a given Year.
# "Sampled" = any tree records exist for (LOC_NAME, Year).

plot_all_sites_pa_trees <- function(df, pdf_out = NULL) {
  df <- df %>%
    mutate(Year = as.integer(Year))
  
  # global year axis (so unsampled years can be grayed)
  all_years <- sort(unique(df$Year))
  
  # species seen at each site (limits rows per facet to locally observed spp)
  species_by_site <- df %>%
    distinct(LOC_NAME, SciName)
  
  # site-years that were sampled (any tree record)
  site_years <- df %>%
    distinct(LOC_NAME, Year) %>%
    mutate(sampled = TRUE)
  
  # species x year presence within site (present if any row exists)
  pres <- df %>%
    group_by(LOC_NAME, SciName, Year) %>%
    summarise(present = 1L, .groups = "drop")
  
  # full grid: (LOC_NAME, SciName seen at that site) x all years
  pa_full <- species_by_site %>%
    tidyr::crossing(Year = all_years) %>%
    left_join(pres,      by = c("LOC_NAME", "SciName", "Year")) %>%
    left_join(site_years,by = c("LOC_NAME", "Year")) %>%
    mutate(
      fill_color = dplyr::case_when(
        is.na(sampled)                  ~ "unsampled",             # no sampling that year at this site
        !is.na(sampled) & present == 1L ~ "present",               # sampled & seen
        !is.na(sampled) & (is.na(present) | present == 0L) ~ "absent"  # sampled & not seen
      )
    )
  
  base_plot <- ggplot(pa_full, aes(x = factor(Year),
                                   y = forcats::fct_rev(forcats::fct_inorder(SciName)),
                                   fill = fill_color)) +
    geom_tile(color = "grey60", linewidth = 0.2) +
    scale_fill_manual(
      values = c(present = "darkgreen", absent = "red3", unsampled = "grey80"),
      breaks = c("present", "absent", "unsampled"),
      name = "Status",
      labels = c("Present", "Absent", "Not sampled")
    ) +
    labs(x = "Year", y = "Species",
         title = "Tree presence–absence by site (9 sites per page)") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "top",
      axis.text.y = element_text(size = 10)
    )
  
  # pages needed at 3x3 facets
  n_pages <- ggforce::n_pages(
    base_plot + ggforce::facet_wrap_paginate(~ LOC_NAME, ncol = 3, nrow = 3, page = 1, scales = "free_y")
  )
  
  if (!is.null(pdf_out)) {
    pdf(file = pdf_out, width = 11, height = 8.5)
    for (i in seq_len(n_pages)) {
      print(base_plot + ggforce::facet_wrap_paginate(~ LOC_NAME, ncol = 3, nrow = 3, page = i, scales = "free_y"))
    }
    dev.off()
  } else {
    for (i in seq_len(n_pages)) {
      print(base_plot + ggforce::facet_wrap_paginate(~ LOC_NAME, ncol = 3, nrow = 3, page = i, scales = "free_y"))
    }
  }
}

# Examples:
# Display:
plot_all_sites_pa_trees(trees)

# Save to PDF (set explicit path):
plot_all_sites_pa_trees(trees, pdf_out = "~/Downloads/tree_presence_absence_by_site.pdf")
#------------------------
# ----Trends for seedlings 
#--------------------------

abund <- seedlings %>%
  mutate(
    Year = as.integer(Year),
    Abundance = round(Stem_Count * CORRECTION_FACTOR, 1)   # or use Stem_Count if you prefer raw
  ) %>%
  filter(GrowthForm == "Tree") %>%
  group_by(SciName, LOC_NAME, Year) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(Abundance)) %>%
  # require >= 2 years per species–site (removes single-point sites)
  group_by(SciName, LOC_NAME) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  # also require >= 2 points per species overall (so geom_smooth(method = "lm") is valid)
  group_by(SciName) %>%
  filter(n() >= 2) %>%
  ungroup()


ranked_species <- abund %>%
  group_by(SciName) %>%
  summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(rank = row_number())

plot_abundance_range <- function(df_ranked, ranks_vec) {
  sp <- df_ranked %>% filter(rank %in% ranks_vec) %>% pull(SciName)
  
  abund %>%
    filter(SciName %in% sp) %>%
    ggplot(aes(x = Year, y = Abundance, color = SciName)) +
    geom_jitter(size = 3) +
    geom_smooth(se = FALSE, method = "lm", linewidth = 1) +
    scale_y_log10() +
    scale_x_continuous(breaks = c(2015, 2018, 2021, 2024)) +
    scale_color_viridis_d(option = "turbo", end = 0.9) +
    labs(
      title = paste0("Tree abundance trends (ranks ", min(ranks_vec), "–", max(ranks_vec), ")"),
      x = "Year", y = "Stem count (corrected)"
    ) +
    theme_plot() +
    theme(panel.grid.minor = element_blank())
}

# plot the abundance trends (after filtering and ranking)
plot_abundance_range(ranked_species, 1:10)
plot_abundance_range(ranked_species, 11:20)

# Plot all
step <- 10 # number of species per plot
rank_chunks <- split(ranked_species$rank, ceiling(ranked_species$rank / step))

# Plot all
for (r in rank_chunks) {
  print(plot_abundance_range(ranked_species, r))
}

# Save plots #change path
pdf_out <- "~/Downloads/tree_abundance_trends_by_rank.pdf"
pdf(file = pdf_out, width = 10, height = 6)
for (r in rank_chunks) {
  print(plot_abundance_range(ranked_species, r))
}
dev.off()



#---------------------------------
#--------- Trends for Trees -----
#---------------------------------

#------------------------
# ---- Trends for TREES
#------------------------

# Option A (default): corrected abundance = sum of CORRECTION_FACTOR across tagged stems
# Option B (uncomment to use raw stem counts): set Abundance = 1 for each row before summarise

abund_trees <- trees %>%
  mutate(
    Year = as.integer(Year),
    Abundance = CORRECTION_FACTOR        # <-- Option A (default)
    # Abundance = 1                      # <-- Option B (raw stem count per tag)
  ) %>%
  group_by(SciName, LOC_NAME, Year) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(Abundance)) %>%
  group_by(SciName, LOC_NAME) %>%
  filter(n() >= 2) %>%                   # require >= 2 years per species–site
  ungroup() %>%
  group_by(SciName) %>%
  filter(n() >= 2) %>%                   # require >= 2 total points per species
  ungroup()

ranked_species_trees <- abund_trees %>%
  group_by(SciName) %>%
  summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(rank = row_number())

plot_abundance_range_trees <- function(df_ranked, ranks_vec, abund_df = abund_trees) {
  sp <- df_ranked %>% filter(rank %in% ranks_vec) %>% pull(SciName)
  
  abund_df %>%
    filter(SciName %in% sp) %>%
    ggplot(aes(x = Year, y = Abundance, color = SciName)) +
    geom_jitter(size = 3) +
    geom_smooth(se = FALSE, method = "lm", linewidth = 1) +
    scale_y_log10() +
    scale_x_continuous(breaks = c(2015, 2018, 2021, 2024)) +
    scale_color_viridis_d(option = "turbo", end = 0.9) +
    labs(
      title = paste0("Tree abundance trends (ranks ", min(ranks_vec), "–", max(ranks_vec), ")"),
      x = "Year", y = "Stem count (corrected)"
    ) +
    theme_plot() +
    theme(panel.grid.minor = element_blank())
}

# Preview top ranks
plot_abundance_range_trees(ranked_species_trees, 1:10)
plot_abundance_range_trees(ranked_species_trees, 11:20)

# Chunk ranks and iterate all
step <- 10
rank_chunks <- split(ranked_species_trees$rank, ceiling(ranked_species_trees$rank / step))

for (r in rank_chunks) {
  print(plot_abundance_range_trees(ranked_species_trees, r))
}

# Save all to PDF  (set your preferred path)
pdf_out <- "~/Downloads/tree_abundance_trends_by_rank.pdf"
pdf(file = pdf_out, width = 10, height = 6)
for (r in rank_chunks) {
  print(plot_abundance_range_trees(ranked_species_trees, r))
}
dev.off()
#------------------------------------
#------------Canopy Trends ----------
#------------------------------------
canopy <- canopy %>%
  mutate(across(
    c(
      CanopyCoverValue,
      SubcanopyMaxHeight_m,
      SubcanopyCoverValues,
      ShrubMaxHeight_m,
      ShrubCoverValue,
      HerbaceousMaxHeight_m,
      HerbaceousCoverValue,
      NonvascCoverValue,
      Unveg_Bedock_Perc,
      Unveg_Boulder_Perc,
      Unveg_LeafLitterDuff_Perc,
      Unveg_DecayingWood_Perc,
      Unveg_Water_Perc
    ),
    ~ suppressWarnings(as.numeric(.))
  ))

str(canopy)



p_canopy <- ggplot(canopy %>%
    group_by(LOC_NAME) %>%
    filter(dplyr::n_distinct(Year) >= 2) %>%
    ungroup(),
    aes(x = Year, y = HerbaceousCoverValue, color = LOC_NAME, group = LOC_NAME)) +
  geom_line(alpha = 0.6, show.legend = F) +
  geom_point(size = 2, show.legend = F) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, linewidth = 0.75, show.legend = F) +
  scale_y_continuous(name = "SubcanopyMaxHeight_m", limits = c(0, 100)) +
  scale_x_continuous() +
  labs(x = "Year", color = "LOC_NAME") +
  theme_plot() +
  ggtitle("SubcanopyMaxHeight_m") +
  scale_color_viridis_d(option = "turbo") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position  = "right")

p_canopy
