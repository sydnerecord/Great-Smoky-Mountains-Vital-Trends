# Atmospheric Deposition — NTN VWMs by year (core ions)

library(tidyverse)
library(lubridate)
library(janitor)
# --- theme (bigger legend) ---
theme_plot <- function(legend = TRUE) {
  theme(
    panel.grid        = element_blank(),
    panel.background  = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1),
    aspect.ratio      = .75,
    axis.text         = element_text(size = 18, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title        = element_text(size = 18),
    axis.title.y      = element_text(margin = margin(r = 10)),
    axis.title.x      = element_text(margin = margin(t = 10)),
    legend.position   = if (legend) "right" else "none",
    legend.text       = element_text(size = 16),
    legend.title      = element_text(size = 18),
    text              = element_text(family = "Helvetica")
  )
}

# --- helper: set NADP sentinel negatives to NA ---
na_fix <- function(x) { x[x %in% c(-9, -999, -9999)] <- NA; x }

# --- read data (NTN only needed for this figure) ---
ntn <- readr::read_csv(
  '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Atmospheric_Deposition/SummaryData/NTN/NTN-tn11-w-i-mg.csv',
  show_col_types = FALSE
)

# --- compute Volume-Weighted Means by year for core ions ---
core_ions <- c("SO4","NO3","NH4","Cl","Ca","Na","K","Mg")   # Br excluded on purpose

vwm_year <- ntn %>%
  mutate(year = year(dateOn),
         ppt  = na_fix(ppt)) %>%                           # precipitation depth
  mutate(across(all_of(core_ions), na_fix)) %>%
  pivot_longer(all_of(core_ions), names_to = "analyte", values_to = "conc") %>%
  filter(!is.na(conc), !is.na(ppt), ppt > 0) %>%
  group_by(year, analyte) %>%
  summarise(vwm_mgL = sum(conc * ppt) / sum(ppt), .groups = "drop")

# order legend by overall mean VWM (highest first)
ord <- vwm_year %>%
  group_by(analyte) %>%
  summarise(mu = mean(vwm_mgL, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mu)) %>%
  pull(analyte)

vwm_year <- vwm_year %>% mutate(analyte = factor(analyte, levels = ord))

# --- plot ---
p <- ggplot(vwm_year, aes(year, vwm_mgL, color = analyte)) +
  geom_line() +
  geom_point(size = 1) +
  scale_y_log10() +
  scale_color_viridis_d(option = "plasma") +  # dark, high-contrast palette
  labs(x = "Year", y = "Volume-weighted mean (mg/L)", color = "Analyte") +
  theme_plot(legend = TRUE)

print(p)

pdf("~/Downloads/deposition.pdf")
p
dev.off()
# --- AMoN (Ammonia Monitoring Network, gaseous NH3) ---
amon <- read_csv(
  '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Atmospheric_Deposition/SummaryData/AMoN/AMoN-tn01-W-i.csv',
  show_col_types = FALSE
)

# --- CASTNET (Ozone & air chemistry) ---
castnet <- read_csv(
  '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Atmospheric_Deposition/SummaryData/CASTNET/CASTNET_GRSM_export.csv',
  show_col_types = FALSE
)

# --- MDN (Mercury Deposition Network, wet Hg) ---
mdn_tn11 <- read_csv(
  '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Atmospheric_Deposition/SummaryData/MDN/MDN-tn11-W-i.csv',
  show_col_types = FALSE
)

mdn_tn12 <- read_csv(
  '/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Data/Atmospheric_Deposition/SummaryData/MDN/MDN-tn12-W-i.csv',
  show_col_types = FALSE
)


# ==========================
# AMoN (Ammonia Monitoring Network)
# ==========================


# AMoN: build amon_clean and annual means
amon_clean <- amon %>%
  # keep good QA grades (A/B are standard; tweak if you want to keep C)
  filter(QR %in% c("A","B")) %>%
  # ensure numeric concentration (AMoN CONC is µg/m³)
  mutate(conc_ugm3 = as.numeric(CONC),
         year      = year(startDate)) %>%
  # drop missing conc
  filter(!is.na(conc_ugm3))

# average A/T replicates within each deployment window
amon_deploy <- amon_clean %>%
  group_by(SITEID, startDate, endDate, year) %>%
  summarise(conc_ugm3 = mean(conc_ugm3, na.rm = TRUE),
            n_reps    = n(),
            .groups   = "drop")

# annual mean NH3 (µg/m³)
amon_year <- amon_deploy %>%
  group_by(year) %>%
  summarise(nh3_mean = mean(conc_ugm3, na.rm = TRUE),
            n_deploys = n(),
            .groups = "drop")

# quick plot
ggplot(amon_year, aes(year, nh3_mean)) +
  geom_line() +
  geom_point(size = 1) +
  labs(x = "Year", y = "Annual mean NH4 (µg/m³)") +
  theme_plot(legend = FALSE)



# ==========================
# MDN (Mercury Deposition Network)
# ==========================

mdn_year <- bind_rows(mdn_tn11 %>% mutate(site="TN11"),
                      mdn_tn12 %>% mutate(site="TN12")) %>%
  clean_names() %>%
  { dt <- .
  date_nm <- intersect(c("start_date","begin_date","sample_start_date","date_start"), names(dt))[1]
  hg_nm   <- intersect(c("total_mercury_ng_l","total_hg_ng_l","hg_ng_l",
                         "concentration_ng_l","conc_ng_l","conc"), names(dt))[1]
  if (is.na(date_nm) || is.na(hg_nm)) stop("MDN: missing date or Hg column; see names(mdn_tn11)")
  dt %>%
    mutate(datetime = parse_date_time(.data[[date_nm]], orders=c("mdy HMS","mdy HM","mdy","ymd HMS","ymd HM","ymd"), quiet=TRUE),
           year = year(datetime),
           hg_ng_l = parse_number(.data[[hg_nm]])) %>%
    filter(!is.na(year), !is.na(hg_ng_l)) %>%
    group_by(site, year) %>%
    summarise(hg_mean_ng_l = mean(hg_ng_l, na.rm=TRUE), .groups="drop")
  }


ggplot(mdn_year, aes(year, hg_mean_ng_l, color=site)) +
  geom_line() + geom_point(size=1) +
  labs(x="Year", y="Annual mean Hg (ng/L)", color="Site", head = ) +
  scale_color_brewer(palette="Set1") +
  theme_plot(legend = TRUE)
