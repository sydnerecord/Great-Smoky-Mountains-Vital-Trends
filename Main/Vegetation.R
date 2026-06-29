############################################################################
# Vegetation analyses ----
###########################################################################

# Load packages for analyses
library(lubridate)

# Set working directory
setwd("G:/Shared drives/GRSM_CESU")

# Read in tree data
trees <- read.csv("Data/Vegetation/trees.csv", header = TRUE, stringsAsFactors = FALSE)

# Understand temporal extent of data ----
# Parse Event_Date column into separate year and month columns to get a better 
# sense of the temporal coverage of the tree data
# Create column for year of observation
Year <- year(as_datetime(trees$Event_Date))
# Create a column for month of observation
Month <- month(as_datetime(trees$Event_Date))
# Append year and month columns to trees dataframe
trees <- cbind(trees, Year, Month)
# Determine months of sampling
unique(trees$Month)
# Determine years of sampling
unique(trees$Year)

# This dataset includes one year (2015) of 16 prototype plots that have not been
# resampled, one full cycle of vegetation data on the 160 plot network 
# (2016-2020) and one full cycle of soils data (2017-2021). The second cycle of 
# data is in progress with 3/5 years of vegetation data complete 
# (2021, 2023, 2024) and 2/5 years of soils data (2023, 2024). No data were 
# collected in 2022.

# Use subset_by_cycle function to subset trees by first and second cycles
# Cycle 1 vegetation (2016-2020)
cycle1_veg <- subset_by_cycle(trees, 2016, 2020)
# Cycle 2 vegetation (2021, 2023, 2024 — no 2022)
cycle2_veg <- subset_by_cycle(trees, 2021, 2024, exclude_years = 2022)

###############################################
# Basal Area Growth ----
###############################################

# Match up Cycle t and Cycle t+1 data to calculate basal area growth per 
# individual tree. To do this, we want to make sure that the t+1 cycle data do 
# not have duplicates for a unique tree tag and location name for an individual.
# If there are duplicates in the second cycle, then the mean of the cycle 2 
# values for the individual are taken for the DBH values from which basal area
# growth will be calculated. Here cycle t is cycle 1 and cycle t+1 is cycle 2

dbh_lookup <- aggregate(DBH ~ Tag + LOC_NAME, data = cycle2_veg, FUN = mean, na.rm = TRUE)
names(dbh_lookup)[names(dbh_lookup) == "DBH"] <- "DBHnext"

# Add column DBHnext to cycle t data to reflect cycle t+1 dbh
cycle1_veg <- merge(cycle1_veg, dbh_lookup, by = c("Tag", "LOC_NAME"), all.x = TRUE)

# Calculate mean Year from cycle2_veg per Tag + Loc_Name to get average remeasurement year
# Average is used to account for any trees that show up twice in a cycle
year_lookup <- aggregate(Year ~ Tag + LOC_NAME, data = cycle2_veg, FUN = mean, na.rm = TRUE)
names(year_lookup)[names(year_lookup) == "Year"] <- "Year_cycle2"

# Merge mean cycle 2 year onto cycle1_veg
cycle1_veg <- merge(cycle1_veg, year_lookup, by = c("Tag", "LOC_NAME"), all.x = TRUE)

# Calculate mean cycle 1 year per Tag + Loc_Name
year_lookup_c1 <- aggregate(Year ~ Tag + LOC_NAME, data = cycle1_veg, FUN = mean, na.rm = TRUE)
names(year_lookup_c1)[names(year_lookup_c1) == "Year"] <- "Year_cycle1"
cycle1_veg <- merge(cycle1_veg, year_lookup_c1, by = c("Tag", "LOC_NAME"), all.x = TRUE)

# Calculate number of years between measurements
cycle1_veg$Years_between <- cycle1_veg$Year_cycle2 - cycle1_veg$Year_cycle1

# Calculate basal area for cycle 1 and cycle 2 (DBH in cm, BA in cm²)
cycle1_veg$BA_cycle1 <- pi * (cycle1_veg$DBH / 2)^2
cycle1_veg$BA_cycle2 <- pi * (cycle1_veg$DBHnext / 2)^2

# Calculate annualized basal area growth per individual (cm² per year)
cycle1_veg$BA_growth_annual <- (cycle1_veg$BA_cycle2 - cycle1_veg$BA_cycle1) / cycle1_veg$Years_between

# Flag negative growth (TRUE = negative growth, likely measurement error)
cycle1_veg$Negative_growth_flag <- ifelse(!is.na(cycle1_veg$BA_growth_annual) & 
                                            cycle1_veg$BA_growth_annual < 0, TRUE, FALSE)

# Summary of results
# Check basic summary statistics to determine if outliers exist. Note some values
# with large negative growth (min) indicating possible measurement error.
summary(cycle1_veg$BA_growth_annual) 

# Exclude values with NA for growth. This indicates individuals for which there 
# was not a cycle t+1 remeasurement yet.
cycle1_veg_noNA <- cycle1_veg[!is.na(cycle1_veg$BA_growth_annual), ]
dim(cycle1_veg_noNA)[1]

# Determine percent of individuals flagged with negative growth
cycle1_veg_neg <- cycle1_veg_noNA[cycle1_veg_noNA$Negative_growth_flag==TRUE,]
cycle1_veg_neg <- cycle1_veg_neg[cycle1_veg_neg$BA_growth_annual<0,]
dim(cycle1_veg_neg)[1]/dim(cycle1_veg_noNA)[1] #21% of individuals show negative growth

# Determine range of negative growth values
range(cycle1_veg_neg$BA_growth_annual)
# Determine range of positive growth values
range(cycle1_veg_noNA[cycle1_veg_noNA$Negative_growth_flag==FALSE,]$BA_growth_annual)

# Visualize outliers in growth 
boxplot(cycle1_veg_noNA$BA_growth_annual, ylab = "Basal Area Growth (cm)")


