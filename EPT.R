
library(tidyverse)
library(RColorBrewer)
library(readxl)

theme_plot <- function(legend = TRUE) {
  ggplot2::theme(
    panel.grid        = element_blank(), 
    aspect.ratio      = .75,
    axis.text         = element_text(size = 18, color = "black"), 
    axis.ticks.length = unit(0.2, "cm"),
    axis.title        = element_text(size = 18),
    axis.title.y      = element_text(margin = margin(r = 10)),
    axis.title.x      = element_text(margin = margin(t = 10)),
    axis.title.x.top  = element_text(margin = margin(b = 5)),
    plot.title        = element_text(size = 18, face = "plain", hjust = 10),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background  = element_blank(),
    strip.background  = element_blank(),
    legend.text       = element_text(size = 15),   
    legend.title      = element_text(size = 18),   
    legend.position   = if (legend) "right" else "none",
    text              = element_text(family = "Helvetica")
  )
}


inverts = read_csv('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export.csv')

# remove white space around values
inverts <- inverts %>%
  mutate(across(everything(), ~ if(is.character(.x)) trimws(.x) else .x))

#Add helpful columns
inverts$Location = word(inverts$LOC_NAME, 1, sep = ",")
inverts$Site = word(inverts$LOC_NAME, 2, sep = ",")
inverts$Genus = word(inverts$Lab_Scientific_Name, 1, sep = " ")
inverts$Year = as.numeric(word(inverts$Start_Date , 1, sep = "-"))
str(inverts)
inverts$Year 


ept_summary <- inverts %>%
  group_by(Location, Year) %>%
  summarise(
    E_rich = n_distinct(Genus[Lab_Order == "Ephemeroptera"], na.rm = TRUE),
    P_rich = n_distinct(Genus[Lab_Order == "Plecoptera"], na.rm = TRUE),
    T_rich = n_distinct(Genus[Lab_Order == "Trichoptera"], na.rm = TRUE),
    EPT_rich = E_rich + P_rich + T_rich,
    .groups = "drop"
  )


library(dplyr)

# 1️⃣ Site-level EPT richness
ept_site_year <- inverts %>%
  group_by(Location, Site, Year) %>%
  summarise(
    E_rich = n_distinct(Genus[Lab_Order == "Ephemeroptera"], na.rm = TRUE),
    P_rich = n_distinct(Genus[Lab_Order == "Plecoptera"], na.rm = TRUE),
    T_rich = n_distinct(Genus[Lab_Order == "Trichoptera"], na.rm = TRUE),
    EPT_rich = E_rich + P_rich + T_rich,
    .groups = "drop"
  )

# 2️⃣ Average across sites within each stream-year
ept_stream_year <- ept_site_year %>%
  group_by(Location, Year) %>%
  summarise(
    E_rich = mean(E_rich, na.rm = TRUE),
    P_rich = mean(P_rich, na.rm = TRUE),
    T_rich = mean(T_rich, na.rm = TRUE),
    EPT_rich = mean(EPT_rich, na.rm = TRUE),
    n_sites = dplyr::n_distinct(Site),
    .groups = "drop"
  )


# Set the minimum number of years required
n_years_min <- 5  # change as needed

# Filter to streams with ≥ n_years_min years of data
streams_enough_years <- ept_stream_year %>%
  group_by(Location) %>%
  filter(dplyr::n_distinct(Year) >= n_years_min) %>%
  ungroup()

# Plot
library(RColorBrewer)

library(viridisLite)

# palette sized to your data
n_streams <- length(unique(streams_enough_years$Location))
stream_palette <- viridisLite::turbo(n_streams)

ggplot(streams_enough_years %>% filter(EPT_rich > 0),
       aes(x = Year, y = EPT_rich, color = Location, group = Location)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = stream_palette, name = "Stream") +   # <-- only ONE color scale
  scale_y_log10(limits = c(18, 60)) +
  labs(x = "Year", y = "EPT Richness") +
  geom_hline(yintercept = c(35, 28, 19), color = "grey40", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 1981, y = 37, label = "Excellent?",  color = "grey30", size = 4, hjust = 0) +
  annotate("text", x = 1981, y = 29, label = "Good?",       color = "grey30", size = 4, hjust = 0) +
  annotate("text", x = 1981, y = 20, label = "Good–Fair?",  color = "grey30", size = 4, hjust = 0) +
  theme_plot()


ept_ratio_stream_year <- inverts %>%
  group_by(Location, Site, Year) %>%
  summarise(
    ept_n  = sum(Count[Lab_Order %in% c("Ephemeroptera","Plecoptera","Trichoptera")], na.rm = TRUE),
    chiro_n = sum(Count[Lab_Family == "Chironomidae"], na.rm = TRUE),
    ratio_EPT_Chiro = ifelse((ept_n + chiro_n) > 0, ept_n / (ept_n + chiro_n), NA_real_),
    .groups = "drop"
  ) %>%
  group_by(Location, Year) %>%
  summarise(
    ratio_EPT_Chiro = mean(ratio_EPT_Chiro, na.rm = TRUE),
    n_sites = dplyr::n_distinct(Site),
    .groups = "drop"
  )


ept_pct_stream_year <- inverts %>%
  group_by(Location, Site, Year) %>%
  summarise(
    total_n = sum(Count, na.rm = TRUE),
    ept_n   = sum(Count[Lab_Order %in% c("Ephemeroptera","Plecoptera","Trichoptera")], na.rm = TRUE),
    pct_EPT = ifelse(total_n > 0, 100 * ept_n / total_n, NA_real_),
    .groups = "drop"
  ) %>%
  group_by(Location, Year) %>%
  summarise(
    pct_EPT = mean(pct_EPT, na.rm = TRUE),
    n_sites = dplyr::n_distinct(Site),
    .groups = "drop"
  )

ept_pooled <- inverts %>%
  group_by(Location, Year) %>%
  summarise(
    total_n = sum(Count, na.rm = TRUE),
    ept_n   = sum(Count[Lab_Order %in% c("Ephemeroptera","Plecoptera","Trichoptera")], na.rm = TRUE),
    pct_EPT_pooled = 100 * ept_n / total_n,
    .groups = "drop"
  )

# Set your minimum number of years required
n_years_min <- 5  # change as needed

# Filter to streams with at least n_years_min years of data
ept_pct_filtered <- ept_pct_stream_year %>%
  group_by(Location) %>%
  filter(dplyr::n_distinct(Year) >= n_years_min) %>%
  ungroup()

# Plot example: %EPT over time
ggplot(ept_pct_filtered, aes(x = Year, y = pct_EPT, color = Location, group = Location)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = stream_palette, name = "Stream") +   # <-- only ONE color scale
  geom_point(size = 2) +
  labs(
    x = "Year",
    y = "% EPT Abundance",
    color = "Stream"
  ) +
  theme_plot() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )



# NCBI 

ncbi = read_xlsx('/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Aquatics_Macroinverts/NCBI/NCBI_taxa.xlsx')
ncbi$Genus = word(ncbi$Species, 1, sep = " ")
unique(ncbi$Genus)

# Get mean NCBI per Genus
ncbi_genus = ncbi %>%
  group_by(Genus) %>%
  summarize(genus_ncbi = mean(NCBI, na.rm = T))

missing_genera <- inverts %>%
  dplyr::distinct(Genus) %>%
  dplyr::anti_join(ncbi_genus, by = "Genus")

#Fraction missing genera
nrow(missing_genera) / nrow(inverts %>% distinct(Genus)) #57%

inv_join <- inverts %>%
  dplyr::left_join(ncbi_genus, by = "Genus")

ncbi_site_year <- inv_join %>%
  dplyr::group_by(Location, Site, Year) %>%
  dplyr::summarize(
    NCBI = sum(genus_ncbi * Count, na.rm = TRUE) /
      sum(Count[!is.na(genus_ncbi)], na.rm = TRUE),
    n_genera_used    = dplyr::n_distinct(Genus[!is.na(genus_ncbi)]),
    counts_used      = sum(Count[!is.na(genus_ncbi)], na.rm = TRUE),
    counts_total     = sum(Count, na.rm = TRUE),
    prop_counts_used = counts_used / counts_total,
    .groups = "drop"
  ) %>%
  dplyr::arrange(Location, Site, Year)

min_years <- 5  # set your minimum number of years

ncbi_stream_year <- ncbi_site_year %>%
  dplyr::group_by(Location, Year) %>%
  dplyr::summarize(
    NCBI_mean = mean(NCBI, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  dplyr::mutate(n_years = dplyr::n_distinct(Year)) %>%
  dplyr::filter(n_years >= min_years) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Location, Year)


ref_lines <- data.frame(
  y = c(4.18, 5.09, 5.91, 7.05),
  label = c("Excellent", "Good", "Good–Fair", "Poor")
)

x_left <- min(ncbi_stream_year$Year, na.rm = TRUE)

p <- ggplot(ncbi_stream_year, aes(x = Year, y = NCBI_mean, color = Location, group = Location)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(1, 5.3)) +
  scale_color_manual(values = stream_palette, name = "Stream") +
  labs(x = "Year", y = "Mean NCBI", color = "Stream") +
  theme_plot() +
  theme(legend.position = "right", panel.grid.minor = element_blank()) +
  geom_hline(
    data = ref_lines,
    aes(yintercept = y),
    color = "grey40", linetype = "dashed", linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = ref_lines,
    aes(x = x_left, y = y + 0.1, label = label),
    color = "grey30", size = 4, hjust = 0,
    inherit.aes = FALSE
  )

p
