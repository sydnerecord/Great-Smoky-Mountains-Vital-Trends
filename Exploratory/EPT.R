# ======================================================================
# Libraries
# ======================================================================
library(tidyverse)     # data import (readr), manipulation (dplyr), plotting (ggplot2), strings (stringr)
library(RColorBrewer)  # palettes (not strictly required below, but you load it—kept as-is)
library(readxl)        # read_xlsx() for NCBI file
library(viridisLite)   # color options
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
    plot.title        = element_text(size = 18, face = "plain", hjust = 10),  # NOTE: hjust=10 intentionally pushes title far right
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

# Base directory (adjust as needed)
base_dir <- "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine"

# File paths
inverts_path <- file.path(base_dir, "Data", "Aquatics_Macroinverts", "SummaryData", "Specimen_Data_Export.csv")
ncbi_path    <- file.path(base_dir, "Data", "Aquatics_Macroinverts", "NCBI", "NCBI_taxa.xlsx")

# Read data
inverts <- readr::read_csv(inverts_path)
ncbi    <- readxl::read_xlsx(ncbi_path)

# ----------------------------------------------------------------------
# Whitespace cleanup
# - Trim only character columns; leave numeric/integer columns unchanged
# ----------------------------------------------------------------------
inverts <- inverts %>%
  mutate(across(everything(), ~ if (is.character(.x)) trimws(.x) else .x))

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

# ======================================================================
# Site-level EPT richness, then averaged to stream-year
# - First: Location × Site × Year richness (E, P, T, EPT)
# - Second: mean across sites in a Location for the Year
# - 'n_sites' retained for visibility/weighting if needed later
# ======================================================================
ept_site_year <- inverts %>%
  group_by(Location, Site, Year) %>%
  summarize(
    E_rich   = n_distinct(Genus[Lab_Order == "Ephemeroptera"], na.rm = TRUE),
    P_rich   = n_distinct(Genus[Lab_Order == "Plecoptera"],   na.rm = TRUE),
    T_rich   = n_distinct(Genus[Lab_Order == "Trichoptera"],  na.rm = TRUE),
    EPT_rich = E_rich + P_rich + T_rich,
    .groups  = "drop"
  )
write_csv(ept_site_year, '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/NCBI_EPT/EPT_site_summary.csv')

# Average sites within stream
ept_stream_year <- ept_site_year %>%
  group_by(Location, Year) %>%
  summarize(
    E_rich   = mean(E_rich,   na.rm = TRUE),  # mean richness across sites
    P_rich   = mean(P_rich,   na.rm = TRUE),
    T_rich   = mean(T_rich,   na.rm = TRUE),
    EPT_rich = mean(EPT_rich, na.rm = TRUE),
    n_sites  = dplyr::n_distinct(Site),       # number of sites contributing
    .groups  = "drop"
  )
write_csv(ept_stream_year, '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/NCBI_EPT/EPT_summary.csv')
# ======================================================================
# Filter to streams with sufficient temporal coverage
# - n_years_min can be tuned; keeps time-series reasonably comparable
# ======================================================================
n_years_min <- 5

streams_enough_years <- ept_stream_year %>%
  group_by(Location) %>%
  filter(dplyr::n_distinct(Year) >= n_years_min) %>%
  ungroup()

n_streams      <- length(unique(streams_enough_years$Location))
stream_palette <- viridisLite::turbo(n_streams)

# ======================================================================
# Plot: EPT richness over time (log-scaled y)
# - Reference lines (35, 28, 19) with tentative labels; see Table 3, p 21 'GRSM_CESU/Maine/Literature/EPT/NCDWRMacroinvertebrate-SOP-February 2016_final.pdf'
# - Log scale emphasizes relative change; limits chosen to 18–60
#   NOTE: Points below 18 or above 60 will be dropped by the limits.
# ======================================================================
ggplot(
  streams_enough_years %>% filter(EPT_rich > 0),         # drop zeros before log scale
  aes(x = Year, y = EPT_rich, color = Location, group = Location)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = stream_palette, name = "Stream") +  # single, manual color scale
  scale_y_log10(limits = c(18, 60)) +                             # NOTE: adjust if you see warnings about removed values
  labs(x = "Year", y = "EPT Richness") +
  geom_hline(yintercept = c(35, 28, 19), color = "grey40", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 1981, y = 37, label = "Excellent?",  color = "grey30", size = 4, hjust = 0) +
  annotate("text", x = 1981, y = 29, label = "Good?",       color = "grey30", size = 4, hjust = 0) +
  annotate("text", x = 1981, y = 20, label = "Good–Fair?",  color = "grey30", size = 4, hjust = 0) +
  theme_plot()

# ======================================================================
# NCBI (NORTH CAROLINA BIOTIC INDEX) join by Genus
# - Read taxa->NCBI table, compute mean NCBI per Genus (if multiple rows/species)
# - Join to specimens by Genus; compute site-year NCBI as abundance-weighted mean
# - Track coverage: n_genera_used, prop_counts_used
# ======================================================================
ncbi$Genus <- stringr::word(ncbi$Species, 1, sep = " ")
unique(ncbi$Genus)   # quick glance at parsed genera

# Mean NCBI per Genus (handles multi-species genera in the table)
ncbi_genus <- ncbi %>%
  group_by(Genus) %>%
  summarize(genus_ncbi = mean(NCBI, na.rm = TRUE), .groups = "drop")

# Which genera in specimens lack an NCBI value after the join base?
length(unique(inverts$Genus))
missing_genera <- inverts %>%
  dplyr::distinct(Genus) %>%
  dplyr::anti_join(ncbi_genus, by = "Genus")

# Fraction of specimen genera missing NCBI (coverage diagnostic)
nrow(missing_genera) / nrow(inverts %>% distinct(Genus)) *100  #=> 57% missing

# Attach genus-level NCBI to specimens
inv_join <- inverts %>%
  dplyr::left_join(ncbi_genus, by = "Genus")

# ----------------------------------------------------------------------
# Site-year NCBI
# - Weighted mean by counts, using only rows where genus_ncbi is present
# - Denominator uses counts where genus_ncbi is non-NA (prevents downward bias)
# - Also record coverage metrics (n genera used, count coverage)
# ----------------------------------------------------------------------
ncbi_site_year <- inv_join %>%
  dplyr::group_by(Location, Site, Year) %>%
  dplyr::summarize(
    NCBI = {
      num <- sum(genus_ncbi * Count, na.rm = TRUE)          # only scored taxa contribute
      den <- sum(Count, na.rm = TRUE)                        # all individuals (scored + unscored)
      ifelse(den > 0, num / den, NA_real_)},
    n_genera_used    = dplyr::n_distinct(Genus[!is.na(genus_ncbi)]),
    counts_used      = sum(Count[!is.na(genus_ncbi)], na.rm = TRUE),
    counts_total     = sum(Count, na.rm = TRUE),
    prop_counts_used = counts_used / counts_total,   # coverage in abundance space
    .groups = "drop") %>%
  dplyr::arrange(Location, Site, Year)
ncbi_site_year$Site <- str_squish(ncbi_site_year$Site)

ncbi_site_year 
write_csv(ncbi_site_year, '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/NCBI_EPT/NCBI_site_results.csv')

# ----------------------------------------------------------------------
# Stream-year NCBI
# - Mean of site-year NCBI values (equal-weight sites within stream)
# - Keep only streams with ≥ min_years of coverage
# ----------------------------------------------------------------------
min_years <- 5

ncbi_stream_year <- ncbi_site_year %>%
  dplyr::group_by(Location, Year) %>%
  dplyr::summarize(
    NCBI_mean = mean(NCBI, na.rm = TRUE),  # equal-weight site means
    .groups   = "drop_last") %>%
  dplyr::mutate(n_years = dplyr::n_distinct(Year)) %>%  # NOTE: within grouped data this will be 1; retained for symmetry
  #dplyr::filter(n_years >= min_years) %>%               # keeps consistency with earlier filters
  dplyr::ungroup() %>%
  dplyr::arrange(Location, Year)

write_csv(ncbi_stream_year, '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/NCBI_EPT/NCBI_results.csv')
# ----------------------------------------------------------------------
# Reference lines and labels for NCBI plot
# - Thresholds: 4.18, 5.09, 5.91, 7.05 (example categories)
# - Labels placed at left edge (x_left)
# ----------------------------------------------------------------------
ref_lines <- data.frame(
  y     = c(4.18, 5.09, 5.91, 7.05),
  label = c("Excellent", "Good", "Good–Fair", "Poor"))

x_left <- min(ncbi_stream_year$Year, na.rm = TRUE)

# ======================================================================
# Plot: Stream-year mean NCBI over time
# - Manual color scale with your palette to keep consistency
# - y-limits set to c(1, 5.3) to frame typical NCBI ranges; adjust as needed
# ======================================================================
p <- ggplot(data = ncbi_stream_year,
            aes(x = Year, y = NCBI_mean, color = Location, group = Location)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(1, 5.3)) +
  scale_color_manual(values = stream_palette, name = "Stream") +
  labs(x = "Year", y = "Mean NCBI", color = "Stream") +
  theme_plot() +
  theme(legend.position  = "right",panel.grid.minor = element_blank()) +
  # Horizontal reference thresholds
  geom_hline(data = ref_lines, aes(yintercept = y), color = "grey40", 
             linetype = "dashed", linewidth = 0.8,inherit.aes = FALSE) +
  # Reference labels at left margin
  geom_text(
    data = ref_lines, aes(x = x_left, y = y + 0.1, label = label),
    color  = "grey30", size = 4, hjust = 0, inherit.aes = FALSE)

p  # draw the plot
