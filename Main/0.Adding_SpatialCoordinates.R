#1.  Setup Workspace -----
# Add locations and check for duplicates
library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(stringr)

# Base path - update to your computer 
general_path <- "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine"  

# Helper to shorten file.path calls
gp <- function(...) file.path(general_path, ...)

# Consider using R Project GitHub/Great-Smoky-Mountains-Vital-Trends/Great-Smoky-Mountains-Vital-Trends.Rproj


# 2. Load & summarize individual datasets -------
## 2.1  FISH (Three-pass data) -------


# Add Location data

# Three pass location data accessed here:
 # https://grsm-nps.opendata.arcgis.com/search?collection=dataset&q=Fisheries
 # click GRSM THREE PASS -> Download -> CSV
 # Csv saved at GRSM_CESU/Maine/Data/Aquatics_Fish/Three_Pass/Locations/GRSM_THREE_PASS.csv

# Download three pass population data and join location to it
three_pass_data <- read_xlsx(
  gp("Data/Aquatics_Fish/Three_Pass/Summary_data/GRSM_Fish_3-Pass_Summary.xlsx"),
  sheet = "Summary") %>%
  rename(STATION_NAME = Site) %>%
  distinct(.) %>% # exclude duplicated rows
  #Join locations
  left_join(
    read_csv(gp("Data/Aquatics_Fish/Three_Pass/Locations/GRSM_THREE_PASS.csv")) %>% 
      select(LAT, LON, WATERSHED, STATION_NAME, STREAMNAME) %>%
      distinct(),
    by = "STATION_NAME") %>%
  rename(LOC_NAME = STATION_NAME) %>%
  distinct(.) # exclude duplicated rows

# Write three pass data with location, watershed and stream name
write_csv(three_pass_data, gp("Data/Aquatics_Fish/Three_Pass/Summary_data/GRSM_Fish_3-Pass_Summary_with_coordinates.csv"))

#   Data Variable description
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


# To get averages...
three_pass_annual_average <- three_pass_data %>%
  mutate(Year = year(Date)) %>%
  group_by(LOC_NAME, Year, WATERSHED, STREAMNAME, LAT, LON) %>%
  summarise(
    Richness = n_distinct(Species),
    ADTDens  = mean(ADTDens, na.rm = TRUE),
    ADTBiom  = mean(ADTBiom, na.rm = TRUE),
    .groups  = "drop"
  )





## 2.2  MACROINVERTEBRATES ------


# download and prepare location file for a join to population data 

#Invertebrate location data from original GRSM folders
invert_loc_codes <- read_csv(gp("Data/Aquatics_Macroinverts/Documents/Locations.csv")) %>%
  mutate(
    station_code = STATION_NAME %>% 
      str_remove("(I&M|I\\&M)$") %>% #edit names to join to data file
      str_replace_all("[^A-Z0-9]", "")) %>%
  select(station_code, STATION_NAME, LOC_NAME, Watershed, StreamName, LAT, LON) %>%
  distinct()

# Join location to original data and remove duplicate rows
inverts_original <- read_csv(gp("Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export.csv")) 

# Response variables = Count (number of individuals per taxon)
# Categorical variables station code (sampling lcoation), LOC_NAME more descriptive but sometimes has multiple sampling codes


# Add spatial coordinates and clean up -------


inverts = inverts_original %>%
  left_join(invert_loc_codes, by = "LOC_NAME")  %>%
  distinct(.)  # remove duplicate rows

 
# Whitespace cleanup: Trim only character columns; leave numeric/integer columns unchanged
inverts <- inverts %>%
  mutate(across(where(is.character), str_squish))

# Add helpful columns - Stream location, site, Genus, and year
inverts$Location <- stringr::word(inverts$LOC_NAME, 1, sep = ",")
inverts$Site     <- stringr::word(inverts$LOC_NAME, 2, sep = ",")
inverts$Genus    <- stringr::word(inverts$Lab_Scientific_Name, 1, sep = " ")
inverts$Year     <- as.integer(stringr::word(inverts$Start_Date, 1, sep = "-"))


# compare to original data lacks lat lon, some duplicates
nrow(inverts_original)
nrow(inverts)

# show duplicates
dups <- inverts_original %>%
  mutate(Year =as.integer(format(EventDate, "%Y"))) %>%
  group_by(across(everything())) %>% # group by every column, only duplicates will appear twice
  filter(n() > 1) %>% #within each group, n() counts how many rows there are. Filter to more than 1
  ungroup() %>%
  arrange(across(everything())) #sets output so that duplicates are next to each other - note if not, its not obvious
dups #492 rows of duplicates - ie, 246 rows have a duplicate

dups %>%
  count(Year, name = "n_duplicates") %>%
  arrange(Year)
# Response variable = Count (number of individuals per taxon)
# Categorical variables sampling station code (sampling location), LOC_NAME more descriptive but sometimes has multiple sampling codes
# Sample_Code = station code plus year

#save data file with added locations
write_csv(inverts, gp("Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export_with_coordinates.csv"))

###
## 2.3  VEGETATION (trees + seedlings) ----------
###
forest_loc <- read_csv(gp("Data/Forest_Health/Locations.csv")) %>%
  select(LOC_NAME, VS_WATERSHED, LAT, LON) #note: Watershed names differ from the official USGS HUC10 (watershed) hydrologic unit names

trees_loc <- read_csv(gp("Data/Forest_Health/Trees.csv")) %>% 
  #Add location and watershed info
  left_join(forest_loc, by = "LOC_NAME")  %>% #Add year
  mutate(Year = year(Event_Date))

seedlings_loc <- read_csv(gp("Data/Forest_Health/Seedlings.csv")) %>%
  left_join(forest_loc, by = "LOC_NAME") %>%
  mutate(Year = year(Event_Date))

wooody_stems_loc <- read_csv(gp("Data/Forest_Health/Woody_Stems.csv")) %>%
  left_join(forest_loc, by = "LOC_NAME") %>%
  mutate(Year = year(Event_Date))

veg_comm_struct_loc <- read_csv(gp("Data/Forest_Health/veg_comm_struct.csv")) %>%
  left_join(forest_loc, by = "LOC_NAME") %>%
  mutate(Year = year(Event_Date))
str(veg_comm_struct_loc)

write_csv(trees_loc, gp("Data/Forest_Health/Trees_with_coordinates.csv"))
write_csv(seedlings_loc, gp("Data/Forest_Health/Seedlings_with_coordinates.csv"))
write_csv(wooody_stems_loc, gp("Data/Forest_Health/Woody_Stems_with_coordinates.csv"))
write_csv(veg_comm_struct_loc, gp("Data/Forest_Health/veg_comm_struct_with_coordinates.csv"))

###
## 2.4  SOILS ---------
###
soil_location = read_csv(gp("Data/Soil_Quality/Soils.csv"))

soils <- read_csv(gp("Data/Soil_Quality/Soils.csv")) %>%
  left_join(forest_loc, by = "LOC_NAME") 

write_csv(soils, gp("Data/Soil_Quality/Soils_with_coordinates.csv"))
