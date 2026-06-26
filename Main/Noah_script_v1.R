library(sf)

setwd("G:/Shared drives/GRSM_CESU")

##
# Water Quality Sites ----
##
NWQ <- read.csv('data/NationalWaterQualityDownload_resultphyschem.csv')
names(NWQ)
results <- NWQ[ ,c('MethodSpeciationName','CharacteristicName',
                'ResultSampleFractionText','ResultMeasureValue',
                'ResultMeasure/MeasureUnitCode')]
names(NWQ)

## compare Water Quality sites to Vital Trends Sites ----
    lat <- NWQ[,grep('Latitude',names(NWQ))]
    lon <- NWQ[,grep('Longitude',names(NWQ))]
    ID <- NWQ[,'MonitoringLocationIdentifier']
    name <- NWQ[,'MonitoringLocationName']
    
    water.sites <- data.frame(ID, name, lat, lon)
    U.water.sites <- unique(water.sites)
    dim(water.sites)
    plot(U.water.sites[,c('lon','lat')], cex=1.5, col='gray', pch=19)
    
    #read in the vital trend datasets with coordinates
    soil <- read.csv('Maine/Data/Soil_Quality/Soils_with_coordinates.csv')
    fish3pass <- read.csv('Maine/Data/Aquatics_Fish/Three_Pass/Summary_data/GRSM_Fish_3-Pass_Summary_with_coordinates.csv')
    veg <- read.csv('Maine/Data/Forest_Health/Locations.csv')
    invert <- read.csv('Maine/Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export_with_coordinates.csv')
    
    #plot them
    points(soil[,c('LON','LAT')], col='brown', pch=19)
    points(fish3pass[,c('LON', 'LAT')],col='red', pch=3)
    points(veg[,c('LON','LAT')],col='blue',pch=2)
    points(invert[,c('LON','LAT')],col='green',pch=4)
    
    #soil and veg are the same, except soil is missing 32 veg points
    sum(!(unique(soil[,'LOC_NAME']) %in% unique(veg[,'LOC_NAME'])))
    sum(!(unique(veg[,'LOC_NAME']) %in% unique(soil[,'LOC_NAME'])))
    
    length(unique(veg[,'LOC_NAME']))
    length(unique(fish3pass[,'LOC_NAME']))
    length(unique(soil[,'LOC_NAME']))
    
    fish.sites <- unique(fish3pass[,c('LOC_NAME','LON','LAT')])
    fish.sites <- na.omit(fish.sites)
    veg.sites <- unique(veg[,c('LOC_NAME','LON','LAT')])
    invert.sites <- unique(invert[,c('LOC_NAME','LON','LAT')])
    
    # Project to UTM zone 17N (crs = 32617
    points_sf_wgs84 <- st_as_sf(fish.sites, coords = c("LON", "LAT"), crs = 4326)
    fish.projected <- st_transform(points_sf_wgs84, crs = 32617)
    points_sf_wgs84 <- st_as_sf(veg.sites, coords = c("LON", "LAT"), crs = 4326)
    veg.projected <- st_transform(points_sf_wgs84, crs = 32617)
    points_sf_wgs84 <- st_as_sf(invert.sites, coords = c("LON", "LAT"), crs = 4326)
    invert.projected <- st_transform(points_sf_wgs84, crs = 32617)
    
    
    ##
# Co-Locating sites across Datasets ----
##

#read in the datasets with coordinates
soil <- read.csv('Maine/Data/Soil_Quality/Soils_with_coordinates.csv')
fish3pass <- read.csv('Maine/Data/Aquatics_Fish/Three_Pass/Summary_data/GRSM_Fish_3-Pass_Summary_with_coordinates.csv')
veg <- read.csv('Maine/Data/Forest_Health/Locations.csv')
invert <- read.csv('Maine/Data/Aquatics_Macroinverts/SummaryData/Specimen_Data_Export_with_coordinates.csv')

#plot them
plot(soil[,c('LAT','LON')])
points(fish3pass[,c('LAT','LON')],col='red')
points(veg[,c('LAT','LON')],col='blue')
points(invert[,c('LAT','LON')],col='green')

#soil and veg are the same, except soil is missing 32 veg points
sum(!(unique(soil[,'LOC_NAME']) %in% unique(veg[,'LOC_NAME'])))
sum(!(unique(veg[,'LOC_NAME']) %in% unique(soil[,'LOC_NAME'])))

length(unique(veg[,'LOC_NAME']))
length(unique(fish3pass[,'LOC_NAME']))
length(unique(soil[,'LOC_NAME']))

fish.sites <- unique(fish3pass[,c('LOC_NAME','LON','LAT')])
fish.sites <- na.omit(fish.sites)
veg.sites <- unique(veg[,c('LOC_NAME','LON','LAT')])
invert.sites <- unique(invert[,c('LOC_NAME','LON','LAT')])

# Project to UTM zone 17N (crs = 32617
points_sf_wgs84 <- st_as_sf(fish.sites, coords = c("LON", "LAT"), crs = 4326)
fish.projected <- st_transform(points_sf_wgs84, crs = 32617)
points_sf_wgs84 <- st_as_sf(veg.sites, coords = c("LON", "LAT"), crs = 4326)
veg.projected <- st_transform(points_sf_wgs84, crs = 32617)
points_sf_wgs84 <- st_as_sf(invert.sites, coords = c("LON", "LAT"), crs = 4326)
invert.projected <- st_transform(points_sf_wgs84, crs = 32617)

veg.sites[15:19,]
fish.sites[1,]
# Distance between each fish site (rows) and each veg site (columns)
dist_fish_veg <- st_distance(fish.projected, veg.projected)
nearest_veg_idx <- st_nearest_feature(fish.projected, veg.projected)
fish.veg.dist <- st_distance(fish.projected, veg.projected[nearest_veg_idx, ], by_element = TRUE)
maxdist <- 400
nearest_veg_idx[ as.numeric(fish.veg.dist) > maxdist ] <- NA

# Distance between each fish site (rows) and each invert site (columns)
dist_fish_invert <- st_distance(fish.projected, invert.projected)
nearest_invert_idx <- st_nearest_feature(fish.projected, invert.projected)
fish.invert.dist <- st_distance(fish.projected, invert.projected[nearest_invert_idx, ], by_element = TRUE)
nearest_invert_idx[ as.numeric(fish.invert.dist) > maxdist ] <- NA

fish.mat <- matrix(nrow=dim(fish.sites)[1],ncol=5,
                   dimnames=list(NULL,c('fish_name','invert_name','invert_dist',
                                 'veg_name','veg_dist')))
fish.mat[,'fish_name'] <- fish.sites[,'LOC_NAME']
fish.mat[,'veg_name'] <- veg.sites[nearest_veg_idx,'LOC_NAME']
fish.mat[,'veg_dist'] <- round(fish.veg.dist)
fish.mat[,'invert_name'] <- veg.sites[nearest_invert_idx,'LOC_NAME']
fish.mat[,'invert_dist'] <- round(fish.invert.dist)
head(fish.mat)

# Distance between each invert site (rows) and each veg site (columns)
dist_invert_veg <- st_distance(invert.projected, veg.projected)
nearest_veg_idx <- st_nearest_feature(invert.projected, veg.projected)
invert.veg.dist <- st_distance(invert.projected, veg.projected[nearest_veg_idx, ], by_element = TRUE)
maxdist <- 400
nearest_veg_idx[ as.numeric(invert.veg.dist) > maxdist ] <- NA


invert.mat <- matrix(nrow=dim(invert.sites)[1],ncol=3,
                   dimnames=list(NULL,c('invert_name',
                                 'veg_name','veg_dist')))
invert.mat[,'invert_name'] <- invert.sites[,'LOC_NAME']
invert.mat[,'veg_name'] <- veg.sites[nearest_veg_idx,'LOC_NAME']
invert.mat[,'veg_dist'] <- round(invert.veg.dist)
head(invert.mat)


write.csv(fish.mat,na = '',file='fish_co-locations_v1.csv', row.names = FALSE)
write.csv(invert.mat,na = '',file='invert_co-locations_v1.csv', row.names = FALSE)
##
#JUNK -----
##
head(fish.mat)

for(i in 1:dim(fish.mat)[1]){
    

}
    
    
        #get the coordinates of the focal fish site
    fish.coord.focal <- fish3pass[fish3pass[ 'LOC_NAME'] == fish.sites[i] , 
                                  c('LON','LAT')][1,]
        inverts[ ,'LON'] - as.numeric(fish.coord.focal['LON'])
        inverts[ ,'LAT'] - as.numeric(fish.coord.focal['LAT'])


