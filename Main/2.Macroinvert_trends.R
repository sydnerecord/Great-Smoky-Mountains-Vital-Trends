# ======================================================================
# Libraries
# ======================================================================
library(RColorBrewer)  # color palettes 
library(readxl)        # read_xlsx() 
library(viridisLite)   # color options
library(tidyverse)     # data import (readr), manipulation (dplyr), plotting (ggplot2), strings (stringr)
library(sf)            # spatial data
library(lubridate)
library(forcats)

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

# ------------------------------------------------------------------
# Define root paths (user-specific)
# ------------------------------------------------------------------
drive_path  <- "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU"

# ------------------------------------------------------------------
# Datasets
# ------------------------------------------------------------------

# Main dataset

inverts <- readr::read_csv(file.path(drive_path, "Maine/Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export_with_coordinates.csv")) # locations added

# ===============================================================
# Add NCBI
# ==============================================================
ncbi <- readxl::read_xlsx(file.path(drive_path, "Maine/Data/Aquatics_Macroinverts/NCBI_EPT/NCBI_taxa.xlsx"))

# add genus to NCBI from species name
ncbi$Genus <- stringr::word(ncbi$Species, 1, sep = " ")
unique(ncbi$Genus)   # quick glance at parsed genera


# Get average ncbi value of species per genus
ncbi_genus <- ncbi %>%
  group_by(Genus) %>%
  summarise(
    Order        = first(str_split_fixed(na.omit(Order), ":", 2)[, 1]), # keep order - first entry per genus and first word -  removing extra words like this Diptera: Chironomidae 
    genus_ncbi   = mean(NCBI, na.rm = TRUE),
    n_species    = n(),
    n_scored     = sum(!is.na(NCBI)),
    .groups = "drop")

# Add NCBI to inverts
inverts <- inverts %>%
  left_join(ncbi_genus %>% dplyr::select(Genus, genus_ncbi), by = "Genus")


# ----------------------------------------------------------------------
# ---------------- Sampling Checks ---------------------------------------
# ----------------------------------------------------------------------
# Sampling over the years - unique sample codes per month

time_table_samplecodes <- inverts %>%
  filter(!is.na(Start_Date), !is.na(Sample_Code)) %>%
  transmute(
    Year  = year(Start_Date),
    Month = factor(month(Start_Date), levels = 1:12, labels = month.abb),
    Sample_Code
  ) %>%
  distinct(Year, Month, Sample_Code) %>%          # one count per unique sample code in a Year–Month
  count(Year, Month, name = "n_samplecodes") %>%  # totals per Year–Month
  pivot_wider(
    names_from  = Month,
    values_from = n_samplecodes,
    values_fill = 0
  ) %>%
  arrange(Year) %>%
  relocate(Year, all_of(month.abb))

time_table_samplecodes %>% print(n = 27)


# by site, month and year
time_table2 <- inverts %>%
  dplyr::filter(!is.na(Start_Date)) %>%
  dplyr::distinct(
    LOC_NAME,
    Year  = lubridate::year(Start_Date),
    Month = lubridate::month(Start_Date, label = TRUE, abbr = TRUE),
    Sample_Code
  ) %>%                                  # each Sample_Code counts once within LOC_NAME × Year × Month
  dplyr::count(LOC_NAME, Year, Month, name = "n_samples") %>%
  tidyr::pivot_wider(
    names_from  = Month,                  # creates Jan..Dec columns
    values_from = n_samples,
    values_fill = 0
  ) %>%
  dplyr::arrange(LOC_NAME, Year) %>%
  dplyr::relocate(LOC_NAME, Year, dplyr::all_of(month.abb))

time_table2%>% print(n = 36)

# Multiple samples?
# Count number of distinct sample codes per site-year
multisite_samples <- inverts %>%
  group_by(LOC_NAME, Year) %>%
  summarise(
    n_samplecodes = n_distinct(Sample_Code),
    n_days        = n_distinct(Start_Date),
    sample_codes  = paste(unique(Sample_Code), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_samplecodes > 1)

# 
multisite_samples %>% print(n = 11)

# ----------------------------------------------------------------------
#Plot by month
# ----------------------------------------------------------------------

inverts %>%
  mutate(
    Stream = stringr::word(LOC_NAME, 1, sep = ","),
    Year   = year(Start_Date)
  ) %>%
  group_by(Stream) %>%
  filter(n_distinct(Year) >= 5) %>%              # keep only ≥5 years sampled
  ungroup() %>%
  mutate(Month = factor(month(Start_Date, label = TRUE, abbr = TRUE),
                        levels = month.abb)) %>%
  ggplot(aes(x = Month, fill = Stream)) +
  geom_bar(position = "stack") +
  scale_fill_viridis_d(option = "turbo") +
  labs(x = "Month", y = "Number of Samples", fill = "Stream",
       title = "Samples by Month (≥5-Year Streams)") +
  theme_plot() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

# scatterplot - time of year samples
inverts %>%
  mutate(
    Stream = stringr::word(LOC_NAME, 1, sep = ","),
    Year   = year(Start_Date),
    doy    = yday(Start_Date)
  ) %>%
  group_by(Stream) %>%
  filter(n_distinct(Year) >= 5) %>%
  ungroup() %>%
  ggplot(aes(x = Year, y = doy, color = Stream)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_continuous(
    breaks = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349),
    labels = month.abb
  ) +
  scale_x_continuous(breaks = seq(1980, 2025, by = 5)) +
  scale_color_viridis_d(option = "turbo") +
  labs(x = "Year", y = "Sampling Month", color = "Stream",
       title = "Sampling Dates by Year (≥5-Year Streams)") +
  theme_plot() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))


# scatterplot - time of year samples
inverts %>%
  mutate(
    Stream = stringr::word(LOC_NAME, 1, sep = ","),
    Year   = year(Start_Date),
    doy    = yday(Start_Date)
  ) %>%
  group_by(Stream) %>%
  filter(n_distinct(Year) >= 5) %>%
  ungroup() %>%
  ggplot(aes(y = Year, x = doy, color = Stream)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_continuous(
    breaks = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349),
    labels = month.abb
  ) +
  scale_y_continuous(breaks = seq(1980, 2025, by = 5)) +
  scale_color_viridis_d(option = "turbo") +
  labs(x = "Sampling Month", y = "Year", color = "Stream",
       title = "Sampling Dates by Year (≥5-Year Streams)") +
  theme_plot() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))


# ----------------------------------------------------------------------
#------- Get EPT ------------------------------------------------------
# ----------------------------------------------------------------------


EPT_ORDERS <- c("Ephemeroptera","Plecoptera","Trichoptera")

#---- Average in cases of multiple samples per site

# first average for a given sample code (usually only one but not always)
sample_level <- inverts %>%
  group_by(LOC_NAME, Year, Start_Date, Sample_Code) %>%
  summarise(
    abund_sample    = sum(Count, na.rm = TRUE),
    rich_sample     = n_distinct(Genus, na.rm = TRUE),
    EPT_rich_sample = n_distinct(Genus[Lab_Order %in% c("Ephemeroptera","Plecoptera","Trichoptera")], na.rm = TRUE),
    EPT_abund_sample= sum(Count[Lab_Order %in% EPT_ORDERS], na.rm = TRUE),
    ncbi_sample = sum(genus_ncbi * Count, na.rm = TRUE)/sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

# then average if multiple samples per day
site_day <- sample_level %>%
  group_by(LOC_NAME, Year, Start_Date) %>%
  summarise(
    abund_day    = mean(abund_sample,    na.rm = TRUE),
    rich_day     = mean(rich_sample,     na.rm = TRUE),
    EPT_rich_day = mean(EPT_rich_sample, na.rm = TRUE),
    EPT_abund_day = mean(EPT_abund_sample, na.rm = TRUE),
    ncbi_day = mean(ncbi_sample, na.rm = TRUE),
    n_samples    = n(),  # QC: how many samples contributed
    .groups = "drop"
  )
# note - outlier abundance (10x less) for Abrams creeks in 2014. Error?

# Average across days per site
site_year <- site_day %>%
  group_by(LOC_NAME, Year) %>%
  summarise(
    abund_site    = mean(abund_day,    na.rm = TRUE),
    rich_site     = mean(rich_day,     na.rm = TRUE),
    EPT_rich_site = mean(EPT_rich_day, na.rm = TRUE),
    EPT_abund_site = mean(EPT_abund_day, na.rm = TRUE),
    ncbi_site = mean(ncbi_day, na.rm = TRUE),
    n_samples    = n(),  # QC: how many samples contributed
    .groups = "drop"
  )

# average across sites per stream
inverts_stream <- site_year %>%
  mutate(Stream = word(LOC_NAME, 1, sep = ",")) %>%
  group_by(Stream, Year) %>%
  summarise(
    abund_stream    = mean(abund_site,    na.rm = TRUE),
    rich_stream     = mean(rich_site,     na.rm = TRUE),
    EPT_rich_stream = mean(EPT_rich_site, na.rm = TRUE),
    EPT_abund_stream = mean(EPT_abund_site, na.rm = TRUE),
    ncbi_stream = mean(ncbi_site, na.rm = TRUE),
    n_samples    = n(),  # QC: how many samples contributed
    .groups = "drop"
  ) %>%
  group_by(Stream) %>%
  mutate(EPT_prop_abund = EPT_abund_stream/abund_stream,
         n_years_sampled = n_distinct(Year), ) %>%
# total years that stream was sampled
  ungroup()
inverts_stream

# Save
#write_csv(inverts_stream , file.path(drive_path, "Maine/Data/Aquatics_Macroinverts/SummaryData/invert_stream_diversity_with_ncbi.csv"))


# ----------------------------------------------------------------------
# ------ Richness, EPT, NCBI Plots ------
# ----------------------------------------------------------------------

# --- Richness
p_richness <- ggplot(
  inverts_stream %>% filter(n_years_sampled >= 5, abund_stream > 0),
  aes(x = Year, y = rich_stream, color = Stream, group = Stream)) +
  geom_line() +
  geom_point(size = 2) +
  #geom_smooth(method = "lm", linetype = "dashed", se = F, size = .75) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(1980, 2020, by = 10)) +
  labs(x = "Year", y = "Genus Richness per Stream", color = "Stream") +
  theme_plot() +
  ggtitle("Genera Richness") +
  scale_color_viridis_d(option = "turbo") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right") 
p_richness

# Add linear fit
p_richness <- ggplot(
  inverts_stream %>% filter(n_years_sampled >= 5, abund_stream > 0),
  aes(x = Year, y = rich_stream, color = Stream, group = Stream)) +
  geom_line() +
  geom_point(size = 2) +
  geom_smooth(method = "lm", linetype = "dashed", se = F, size = .75) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(1980, 2025, by = 10)) +
  labs(x = "Year", y = "Genus Richness per Stream", color = "Stream") +
  theme_plot() +
  ggtitle("Genera Richness") +
  scale_color_viridis_d(option = "turbo") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right") 
p_richness


# EPT Richness
p_ept <- ggplot(
  inverts_stream %>% 
    filter(n_years_sampled >= 5, abund_stream > 0),
  aes(x = Year, y = EPT_rich_stream , color = Stream, group = Stream)) +
  geom_line() +
  geom_point(size = 2) +
  geom_smooth(method = "lm", linetype = "dashed", se = F, size = .5) +
  scale_x_continuous(
    breaks = seq(1985, 2025, by = 10),
    limits = c(1984, 2025)) +
  labs(x = "Year", y = "EPT", color = "Stream") +
  theme_plot() +
  ggtitle("EPT Richness") +
  scale_color_viridis_d(option = "turbo") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right") 
p_ept

# --- Abundance
p_abundance <- ggplot(
  inverts_stream %>% 
    filter(n_years_sampled >= 5, abund_stream > 0),
  aes(x = Year, y = abund_stream, color = Stream, group = Stream)) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = seq(1985, 2025, by = 5),
    limits = c(1984, 2025)) +
  labs(x = "Year", y = "Genera Abundance per Stream", color = "Stream") +
  theme_plot() +
  ggtitle("Stream Abundance") +
  scale_color_viridis_d(option = "turbo") +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right") 
p_abundance

# with linear fit
p_abundance <- ggplot(
  inverts_stream %>% 
    filter(n_years_sampled >= 5, abund_stream > 0),
  aes(x = Year, y = abund_stream, color = Stream, group = Stream)) +
  geom_line() +
  geom_point(size = 2) +
  stat_smooth(method = "lm", linetype = "dashed", se = F, show.legend = F, size = 0.5) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = seq(1985, 2025, by = 5),
    limits = c(1984, 2025)) +
  labs(x = "Year", y = "Abundance per Stream", color = "Stream") +
  theme_plot() +
  ggtitle("Macroinvertebrate Abundance") +
  scale_color_viridis_d(option = "turbo") +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
         legend.position = "right") +
  annotate("segment", x = 2017, xend = 2019.5, y = 3100, yend = 3100,
           color = "black", linetype = "dashed") +
  annotate("text", x = 2011.5, y = 3000, label = "Linear Fit", size = 5,
           hjust = 0, vjust = 0, color = "black")
p_abundance


# EPT abundance
p_ept_abundance <- ggplot(
  inverts_stream %>% 
    filter(n_years_sampled >= 5, abund_stream > 0),
  aes(x = Year, y = EPT_abund_stream, color = Stream, group = Stream)) +
  geom_line() +
  geom_point(size = 2) +
  stat_smooth(method = "lm", linetype = "dashed", se = F, show.legend = F, size = 0.5) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = seq(1985, 2025, by = 5),
    limits = c(1984, 2025)) +
  labs(x = "Year", y = "Abundance per Stream", color = "Stream") +
  theme_plot() +
  ggtitle("EPT Stream Abundance") +
  scale_color_viridis_d(option = "turbo") +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
         legend.position = "right") +
  annotate("segment", x = 2017, xend = 2019.5, y = 3100, yend = 3100,
           color = "black", linetype = "dashed") +
  annotate("text", x = 2011.5, y = 3000, label = "Linear Fit", size = 5,
           hjust = 0, vjust = 0, color = "black")
p_ept_abundance

# Relative EPT abundance
p_rel_ept <- ggplot(
  inverts_stream %>% 
    filter(n_years_sampled >= 5, abund_stream > 0),
  aes(x = Year, y = EPT_prop_abun, color = Stream, group = Stream)) +
  geom_line() +
  geom_point(size = 2) +
  stat_smooth(method = "lm", linetype = "dashed", se = F, show.legend = F, size = 0.5) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = seq(1985, 2025, by = 5),
    limits = c(1984, 2025)) +
  labs(x = "Year", y = "Abundance per Stream", color = "Stream") +
  theme_plot() +
  ggtitle("EPT Relative Abundance") +
  scale_color_viridis_d(option = "turbo") +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
         legend.position = "right") 
p_rel_ept

# --- NCBI
ref_lines <- data.frame(
  y     = c(4.18, 5.09, 5.91, 7.05),
  label = c("Excellent", "Good", "Good–Fair", "Poor")) 
x_left <- min(inverts_stream$Year, na.rm = TRUE)

p_ncbi <- ggplot(
  inverts_stream %>% 
    filter(n_years_sampled >= 5, abund_stream > 0),
  aes(x = Year, y = ncbi_stream, color = Stream, group = Stream)) +
  geom_line() +
  geom_point(size = 2) +
  geom_smooth(method = "lm", linetype = "dashed", se = F, size = .75) +
  scale_x_continuous(
    breaks = seq(1985, 2025, by = 5),
    limits = c(1984, 2025)) +
  labs(x = "Year", y = "NCBI", color = "Stream") +
  scale_y_continuous(name = "NCBI per Stream", limits = c(0.5, 4.4)) +
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



#########################
# ------Run over summer months ---
#-------------------------------

# -----------------------------
# 1) Build month-aware summaries
# -----------------------------

# Add Month once (ordered Jan..Dec)
inverts_m <- inverts %>%
  mutate(
    Month = factor(month(Start_Date), levels = 1:12, labels = month.abb),
    Year  = year(Start_Date)
  )

EPT_ORDERS <- c("Ephemeroptera", "Plecoptera", "Trichoptera")

# sample level (unique Sample_Code within day)
sample_level_m <- inverts_m %>%
  group_by(LOC_NAME, Year, Month, Start_Date, Sample_Code) %>%
  summarise(
    abund_sample     = sum(Count, na.rm = TRUE),
    rich_sample      = n_distinct(Genus, na.rm = TRUE),
    EPT_rich_sample  = n_distinct(Genus[Lab_Order %in% EPT_ORDERS], na.rm = TRUE),
    EPT_abund_sample = sum(Count[Lab_Order %in% EPT_ORDERS], na.rm = TRUE),
    ncbi_sample      = sum(genus_ncbi * Count, na.rm = TRUE) / sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

# day level
site_day_m <- sample_level_m %>%
  group_by(LOC_NAME, Year, Month, Start_Date) %>%
  summarise(
    abund_day     = mean(abund_sample,     na.rm = TRUE),
    rich_day      = mean(rich_sample,      na.rm = TRUE),
    EPT_rich_day  = mean(EPT_rich_sample,  na.rm = TRUE),
    EPT_abund_day = mean(EPT_abund_sample, na.rm = TRUE),
    ncbi_day      = mean(ncbi_sample,      na.rm = TRUE),
    n_samples     = n(),
    .groups = "drop"
  )

# site × year × month
site_year_month <- site_day_m %>%
  group_by(LOC_NAME, Year, Month) %>%
  summarise(
    abund_site     = mean(abund_day,     na.rm = TRUE),
    rich_site      = mean(rich_day,      na.rm = TRUE),
    EPT_rich_site  = mean(EPT_rich_day,  na.rm = TRUE),
    EPT_abund_site = mean(EPT_abund_day, na.rm = TRUE),
    ncbi_site      = mean(ncbi_day,      na.rm = TRUE),
    n_days         = n(),
    .groups = "drop"
  )

# stream × year × month
inverts_stream_m <- site_year_month %>%
  mutate(Stream = stringr::word(LOC_NAME, 1, sep = ",")) %>%
  group_by(Stream, Year, Month) %>%
  summarise(
    abund_stream      = mean(abund_site,     na.rm = TRUE),
    rich_stream       = mean(rich_site,      na.rm = TRUE),
    EPT_rich_stream   = mean(EPT_rich_site,  na.rm = TRUE),
    EPT_abund_stream  = mean(EPT_abund_site, na.rm = TRUE),
    ncbi_stream       = mean(ncbi_site,      na.rm = TRUE),
    EPT_prop_abun     = if_else(abund_stream > 0, EPT_abund_stream / abund_stream, NA_real_),
    n_sites_agg       = n(),     # how many sites contributed in that month
    .groups = "drop"
  ) %>%
  group_by(Stream, Month) %>%
  mutate(n_years_sampled = dplyr::n_distinct(Year)) %>%  # total years sampled for this Stream in THIS Month
  ungroup()

# --- Helper function: single wide summer plot ---
summer_plot <- function(df, yvar, ylab, title, logy = FALSE, ylim = NULL) {
  
  df_summer <- df %>%
    filter(Month %in% c("Jun", "Jul", "Aug")) %>%
    group_by(Stream, Month) %>%
    mutate(n_years_sampled = dplyr::n_distinct(Year)) %>%
    ungroup() %>%
    filter(n_years_sampled >= 5) %>%
    mutate(Month = forcats::fct_relevel(Month, "Jun", "Jul", "Aug"))
  
  if (isTRUE(logy)) {
    df_summer <- df_summer %>% filter(.data[[yvar]] > 0)
  }
  
  p <- ggplot(
    df_summer,
    aes(x = Year, y = .data[[yvar]], color = Stream, group = Stream)
  ) +
    geom_line() +
    geom_point(size = 2) +
    geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.7) +
    scale_x_continuous(breaks = seq(1980, 2025, by = 10), limits = c(1984, 2025)) +
    labs(x = "Year", y = ylab, color = "Stream", title = title) +
    theme_plot() +
    scale_color_viridis_d(option = "turbo") +
    facet_wrap(~ Month, ncol = 1, strip.position = "left",
               labeller = as_labeller(c(Jun = "June", Jul = "July", Aug = "August"))) +
    theme(
      strip.text       = element_text(size = 18, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position  = "right"
    )
  
  if (isTRUE(logy)) {
    p <- p + scale_y_log10()
  } else if (!is.null(ylim)) {
    p <- p + scale_y_continuous(limits = ylim)
  }
  
  p + theme(plot.margin = margin(10, 20, 10, 10),
            aspect.ratio = 0.25)
}

# --- Define metrics (one plot each) ---
metrics <- tribble(
  ~yvar,               ~ylab,                          ~title,                           ~logy,  ~ylim,
  "rich_stream",       "Genus Richness per Stream",    "Genera Richness — Summer",       TRUE,   NA,
  "abund_stream",      "Abundance per Stream",         "Stream Abundance — Summer",      TRUE,   NA,
  "EPT_rich_stream",   "EPT Richness per Stream",      "EPT Richness — Summer",          TRUE,  NA,
  "EPT_prop_abun",     "EPT Relative Abundance",       "EPT Relative Abundance — Summer",FALSE,  c(0, 1),
  "ncbi_stream",       "NCBI per Stream",              "NCBI — Summer",                  FALSE,  c(0.5, 4.4)
)

# --- Generate all plots ---
summer_plots <- metrics %>%
  mutate(plot = pmap(list(yvar, ylab, title, logy, ylim),
                     ~ summer_plot(inverts_stream_m, ..1, ..2, ..3, ..4, ..5)))

# --- Access individual plots ---
p_richness_summer   <- summer_plots$plot[[1]]
p_abundance_summer  <- summer_plots$plot[[2]]
p_ept_rich_summer   <- summer_plots$plot[[3]]
p_ept_prop_summer   <- summer_plots$plot[[4]]
p_ncbi_summer       <- summer_plots$plot[[5]]

# --- Display any plot ---
p_richness_summer
p_abundance_summer
p_ept_rich_summer
p_ept_prop_summer
p_ncbi_summer

#------ Attach geometry for plotting

# Spatial layers 
grsm_watershed = st_read(file.path(drive_path, "Maine/Spatial/Watershed/GRSM_watershed.gpkg"))[2]
grsm_streams = st_read(file.path(drive_path, "Maine/Spatial/Streams/GRSM_streams.gpkg"))
grsm_border <- sf::st_read(file.path(drive_path, "Maine/Spatial/BOUNDARY_LN/BOUNDARY_LN.shp"))[6] 

coords_sf <- inverts %>%
  st_as_sf(coords = c("LON","LAT"), crs = 4326, remove = FALSE) %>%
  mutate(stream = str_squish(str_remove(LOC_NAME, ",?\\s*Site\\s*\\w+$"))) %>%
  group_by(stream) %>%
  select(-Year) %>%
  slice_head(n = 1)  # one actual sampled point per stream

target_crs <- 26917 # EPSG:26917, which is NAD83 / UTM zone 17N
grsm_watershed <- st_transform(grsm_watershed, target_crs)
grsm_border    <- st_transform(grsm_border, target_crs)
grsm_streams   <- st_transform(grsm_streams, target_crs)
coords_sf      <- st_transform(coords_sf, target_crs)


inverts_stream_sf <- inverts_stream %>%
  left_join(coords_sf %>% dplyr::rename(Stream = stream), by = "Stream") %>%
  sf::st_as_sf()

inverts_stream_sf <- st_transform(inverts_stream_sf, st_crs(grsm_watershed))


# Minimal prep for colorizing by name
streams_named <- grsm_streams %>%
  filter(!is.na(gnis_name), gnis_name != "") %>%
  mutate(color_id = factor(gnis_name)) %>%
  rename(Stream = gnis_name)


loc_used <- inverts_stream_sf %>%
  st_drop_geometry() %>%
  dplyr::distinct(Stream) 
streams_used <- streams_named %>%
  dplyr::semi_join(loc_used, by = "Stream")


# join locations to coords; keep one point per Location
sites_pts <- coords_sf %>%
  inner_join(loc_used, by = c("stream" = "Stream")) %>%  # keep only sampled sites
  dplyr::select(stream, geometry) %>%
  sf::st_transform(sf::st_crs(grsm_watershed))

# one label point per stream name (keep all stream lines)
stream_labels <- streams_used %>%
  dplyr::group_by(Stream, color_id) %>%
  dplyr::summarise(do_union = TRUE, .groups = "drop") %>%
  sf::st_transform(26917) %>%
  dplyr::mutate(geometry = sf::st_point_on_surface(geom)) %>%  # one point near the middle
  sf::st_transform(sf::st_crs(streams_used))

watershed_labels <- st_point_on_surface(grsm_watershed)

# 4) Plot: watersheds + border + filtered streams + locations
ggplot() +
  geom_sf(data = grsm_border,  fill = NA, color = "gray50", linewidth = 1) +
  geom_sf(data = grsm_watershed, fill = NA, color = "blue", linewidth = .75) +
  geom_sf(data = streams_used, aes(color = Stream), linewidth = .5, alpha = 0.9) +
  geom_sf_text(data = stream_labels, aes(label = Stream, color = Stream),
               size = 3, check_overlap = TRUE) +
  geom_sf(data = sites_pts, aes(color = stream), size = 4) +
  geom_sf_label(data = grsm_watershed, aes(label = name),fill = "white", fontface = "bold", size = 3, color = "blue",label.size = 0,
                label.r = unit(0.15, "lines"), size = 3) +
  scale_color_viridis_d(option = "turbo", guide = "none") +

  coord_sf() +
  labs(title = "GRSM Invert Locations") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank()) 

# No grsm_watershed
ggplot() +
  geom_sf(data = grsm_border,  fill = NA, color = "gray50", linewidth = 1) +
  #geom_sf(data = grsm_watershed,    fill = NA, color = "blue", linewidth = 0.15) +
  geom_sf(data = streams_used, aes(color = Stream), linewidth = .75, alpha = 0.9) +
  geom_sf_text(data = stream_labels, aes(label = Stream, color = Stream),
               size = 3, check_overlap = TRUE) +
  geom_sf(data = sites_pts, aes(color = stream), size = 4) +
  #geom_sf_label(data = grsm_watershed, aes(label = name),fill = "white", fontface = "bold", sieze = 3, color = "blue",label.size = 0,
    #            label.r = unit(0.15, "lines"), size = 3) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  coord_sf() +
  labs(title = "GRSM Invert Locations") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank()) 




#------ Get Slopes -----


# --- Prepare invert responses ---
inverts_stream <- inverts_stream %>%
  mutate(
    log_richness       = log(rich_stream),
    log_abundance      = log(abund_stream),
    log_EPT_richness   = log(EPT_rich_stream),
    log_EPT_abundance  = log(EPT_abund_stream),
    log_rel_EPT_abundance = log(EPT_abund_stream/abund_stream)) %>%
  mutate(
    across(
      c(log_richness, log_abundance, log_EPT_richness, log_EPT_abundance, log_rel_EPT_abundance ),
      ~ ifelse(is.finite(.x), .x, NA_real_)
    )
  )

# --- Responses to analyze ---
responses <- c(
  "log_richness",
  "log_abundance",
  "log_EPT_richness",
  "log_EPT_abundance",
  "log_rel_EPT_abundance",
  "ncbi_stream"   # left unlogged
)

# Optional pretty labels for plot titles
resp_labels <- c(
  log_richness      = "Invertebrate Richness Trends",
  log_abundance     = "Invertebrate Abundance  Trends",
  log_EPT_richness  = "EPT Richness Trends",
  log_EPT_abundance = "EPT Abundance Trends",
  log_rel_EPT_abundance = "EPT Relative Abundance Trends",
  ncbi_stream       = "NCBI Index Trends"
)

# --- Trend helper ---
fit_stream_trend <- function(stream_data, response_var){
  m  <- lm(stats::reformulate("Year", response_var), data = stream_data)
  cs <- coef(summary(m))
  i  <- match("Year", rownames(cs))
  tibble(
    slope   = ifelse(is.na(i), NA_real_, cs[i, "Estimate"]),
    p_value = ifelse(is.na(i), NA_real_, cs[i, "Pr(>|t|)"]),
    r_sq    = summary(m)$r.squared
  )
}

# --- Compute per-stream trends for inverts only ---
lm_per_stream <- purrr::map_dfr(responses, function(var){
  inverts_stream %>%
    filter(n_years_sampled >= 5, is.finite(.data[[var]])) %>%
    group_by(Stream) %>%
    group_modify(~ fit_stream_trend(.x, var)) %>%
    mutate(
      response    = var,
      significant = ifelse(!is.na(p_value) & p_value < 0.05, "yes", "no")
    ) %>%
    ungroup()
})

# --- Join geometry (one geometry per Stream) ---
loc_geom <- inverts_stream_sf %>%
  group_by(Stream) %>%
  slice_head(n = 1) %>%
  dplyr::select(Stream, geometry)

lm_per_stream_sf <- lm_per_stream %>%
  left_join(loc_geom, by = "Stream") %>%
  sf::st_as_sf()

# --- Plot loop over all responses ---
gray_out <- FALSE  # TRUE = gray nonsignificant, FALSE = color all

for (var in responses) {
  
  pts <- lm_per_stream_sf %>% filter(response == var)
  if (nrow(pts) == 0) next
  
  max_abs <- max(abs(pts$slope), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  pts <- pts %>%
    mutate(fill_val = if (gray_out) ifelse(significant == "yes", slope, NA_real_) else slope)
  
  p <- ggplot() +
    geom_sf(data = grsm_border,  fill = NA, color = "gray50", linewidth = 1) +
    geom_sf(data = streams_used, color = "grey70", linewidth = 0.5, alpha = 0.7) +
    geom_sf_text(data = stream_labels, aes(label = Stream), color = "black",
                 size = 3, check_overlap = TRUE, show.legend = FALSE) +
    geom_sf(data = pts, aes(fill = fill_val),
            shape = 21, color = "black", size = 7, stroke = 0.3) +
    scale_fill_gradientn(
      colors = c("#1E90FF", "white", "red"),
      limits = c(-max_abs, max_abs),
      na.value = "grey70",
      name = "Slope"
    ) +
    coord_sf() +
    labs(title = resp_labels[[var]]) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank())
  
  print(p)
}
