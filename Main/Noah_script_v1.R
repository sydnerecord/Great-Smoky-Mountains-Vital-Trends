library(sf)
library(lubridate)
setwd("G:/Shared drives/GRSM_CESU")

##
# Water Quality Sites ----
##
#looking for: 
#    pH
#    ANC
#    Discharge
#    Sulfate and Nitrate

#these data downloaded from https://www.waterqualitydata.us/#advanced=true
#using a bounding box of 35.9, 35.3, -83, -84.1, 

NWQ <- read.csv('data/NationalWaterQualityDownload_resultphyschem.csv')
names(NWQ)
results <- NWQ[ ,c('MethodSpeciationName','CharacteristicName',
                'ResultSampleFractionText','ResultMeasureValue',
                'ResultMeasure.MeasureUnitCode',
                'ResultAnalyticalMethod.MethodIdentifier',
                'ResultAnalyticalMethod.MethodIdentifierContext',
                'ResultAnalyticalMethod.MethodName',
                'ResultAnalyticalMethod.MethodUrl')]

head(results)
names(results)
results[1:20,c(2,8)]
results[1:20,c(2,3)]

nitrate <- results[results[,2]=='Nitrate', c(3,4,8)]
unique(nitrate[,c(1,3)])
summary(lm(as.numeric(nitrate[,2]) ~ factor(nitrate[,1])))

sulfate <- results[results[,2]=='Sulfate', c(3,4,8)]
unique(sulfate[,c(1,3)])
summary(lm(as.numeric(sulfate[,2]) ~ factor(sulfate[,1])))
summary(lm(as.numeric(sulfate[,2]) ~ factor(sulfate[,3])))
        
ph <- results[results[,2]=='pH', c(3,4,8)]
unique(ph[,c(1,3)])
summary(lm(as.numeric(ph[,2]) ~ factor(ph[,3])))

anc <- results[results[,2]=='Gran acid neutralizing capacity', c(3,4,8)]
unique(anc[,c(1,3)])
summary(lm(as.numeric(anc[,2]) ~ factor(anc[,1])))
summary(lm(as.numeric(anc[,2]) ~ factor(anc[,3])))


date <- ymd(NWQ[,'ActivityStartDate'])



## compare Water Quality sites to Vital Trends Sites ----
    lat <- NWQ[,grep('Latitude',names(NWQ))]
    lon <- NWQ[,grep('Longitude',names(NWQ))]
    ID <- NWQ[,'MonitoringLocationIdentifier']
    name <- NWQ[,'MonitoringLocationName']
    
    water.sites <- data.frame(ID, name, lat, lon)
    u.water.sites <- unique(water.sites)
    u.water.sites <- u.water.sites[!is.na(u.water.sites[,'lon']),]
    dim(u.water.sites)
    plot(u.water.sites[,c('lon','lat')], cex=1.3, col='gray', pch=19)
    
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

    points_sf_wgs84 <- st_as_sf(u.water.sites, coords = c("lon", "lat"), crs = 4326)
    water.projected <- st_transform(points_sf_wgs84, crs = 32617)
    
    
    # Distance between each fish site (rows) and each water site (columns)
    dist_fish_water <- st_distance(fish.projected, water.projected)
    nearest_water_idx <- st_nearest_feature(fish.projected, water.projected)
    fish.water.dist <- st_distance(fish.projected, water.projected[nearest_water_idx, ], by_element = TRUE)
    fish.mat <- matrix(nrow=dim(fish.sites)[1],ncol=3,
                        dimnames=list(NULL,c('fish_name','water_name','water_dist')))
    fish.mat[,'fish_name'] <- fish.sites[,'LOC_NAME']
    fish.mat[,'water_name'] <- water.sites[nearest_water_idx,'ID']
    fish.mat[,'water_dist'] <- round(fish.water.dist)
    fish.mat[,'water_dist']
    mean(as.numeric(fish.mat[,'water_dist']))
    
    # Distance between each ivert site (rows) and each water site (columns)
    dist_invert_water <- st_distance(invert.projected, water.projected)
    nearest_water_idx <- st_nearest_feature(invert.projected, water.projected)
    invert.water.dist <- st_distance(invert.projected, water.projected[nearest_water_idx, ], by_element = TRUE)
    invert.mat <- matrix(nrow=dim(invert.sites)[1],ncol=3,
                        dimnames=list(NULL,c('invert_name','water_name','water_dist')))
    invert.mat[,'invert_name'] <- invert.sites[,'LOC_NAME']
    invert.mat[,'water_name'] <- water.sites[nearest_water_idx,'ID']
    invert.mat[,'water_dist'] <- round(invert.water.dist)
    invert.mat[,'water_dist']
    hist(as.numeric(invert.mat[,'water_dist']))
    
### Reduce number of water sites for analysis ----
    
    #make 1000 m buffer around all fish and iverts and crop water to those
    fish.invert.proj <- rbind(fish.projected, invert.projected)
    plot(fish.invert.proj)
    buff <- st_buffer(fish.invert.proj, 1000)
    plot(buff)
    water.buff <- water.projected[buff,]
    unique(water.buff[,'ID'])
    focal.water.sites <- unique(water.buff$ID)
    f.rows <- water.sites[,'ID'] %in% focal.water.sites
    NWQ.crop <- NWQ[f.rows,]
    
    results <- NWQ.crop[ ,c('CharacteristicName',
                'ResultSampleFractionText','ResultMeasureValue',
                'ResultAnalyticalMethod.MethodName',
                'ActivityLocation.LongitudeMeasure',
                'ActivityLocation.LatitudeMeasure',
                'MonitoringLocationIdentifier',
                'MonitoringLocationName',
                'ActivityStartDate')]

## Cut down to a single method per metric -----
    #only one Nitrate Method dominates these samples
    nitrate <- results[results[,'CharacteristicName'] =='Nitrate',]
    table(nitrate[,c(2,4)])
    nitrate <- results[results[,'CharacteristicName'] =='Nitrate' &
                       results[,'ResultSampleFractionText'] =='Filterable' &
                           results[,'ResultMeasureValue'] != '',]
    table(nitrate[,c(2,4)])
    

    #only one Sulfate Method dominates these samples
    sulfate <- results[results[,'CharacteristicName'] =='Nitrate',]
    table(sulfate[,c(2,4)])
    sulfate <- results[results[,'CharacteristicName'] =='Sulfate' &
                       results[,'ResultSampleFractionText'] =='Filterable' &
                           results[,'ResultMeasureValue'] != '',]
    table(sulfate[,c(2,4)])
    
    #pH Method and Fraction text needed to reduce to one dominant method
    ph <- results[results[,'CharacteristicName'] =='pH',]
    table(ph[,c(2,4)])
    ph <- results[results[,'CharacteristicName'] =='pH' &
                       results[,'ResultSampleFractionText'] =='' &
                      results[,'ResultAnalyticalMethod.MethodName'] == 'pH'&
                           results[,'ResultMeasureValue'] != '',]
    table(ph[,c(2,4)])
    
    #anc Method and Fraction text needed to reduce to one dominant method
    anc <- results[results[,'CharacteristicName']=='Gran acid neutralizing capacity', ]
    table(anc[,c(2,4)])
    unique(anc[,2])
    anc <- results[results[,'CharacteristicName'] =='Gran acid neutralizing capacity' &
                       results[,'ResultSampleFractionText'] =='Total' &
                       results[,'ResultAnalyticalMethod.MethodName'] == 'ANC using Gran Titration'&
                           results[,'ResultMeasureValue'] != '',]
    table(anc[,c(2,4)])

# Examine Water metrics over time ----    
plot(as.numeric(ph[,'ResultMeasureValue']) ~   ymd(ph[,'ActivityStartDate']))
        

# Plot Water Metrics over space ----
    metric <- nitrate
    values <- log(as.numeric(metric[,'ResultMeasureValue']))
    scaled_vals <- (values - min(values)) / (max(values) - min(values))
    cols <- heat.colors(100)[round(scaled_vals * 99) + 1]

    plot(jitter(metric[,'ActivityLocation.LongitudeMeasure'],300),
    jitter(metric[,'ActivityLocation.LatitudeMeasure'],300),
           pch=19,
           col=cols)

    metric <- sulfate
    values <- log(as.numeric(metric[,'ResultMeasureValue']))
    scaled_vals <- (values - min(values)) / (max(values) - min(values))
    cols <- heat.colors(100)[round(scaled_vals * 99) + 1]

    plot(jitter(metric[,'ActivityLocation.LongitudeMeasure'],300),
    jitter(metric[,'ActivityLocation.LatitudeMeasure'],300),
           pch=19,
           col=cols)

    
    metric <- ph
    values <- as.numeric(metric[,'ResultMeasureValue'])
    scaled_vals <- (values - min(values)) / (max(values) - min(values))
    cols <- heat.colors(100)[round(scaled_vals * 99) + 1]

    plot(jitter(metric[,'ActivityLocation.LongitudeMeasure'],300),
    jitter(metric[,'ActivityLocation.LatitudeMeasure'],300),
           pch=19,
           col=cols)

    metric <- anc
    values <- as.numeric(metric[,'ResultMeasureValue'])
    scaled_vals <- (values - min(values)) / (max(values) - min(values))
    cols <- heat.colors(100)[round(scaled_vals * 99) + 1]

    plot(jitter(metric[,'ActivityLocation.LongitudeMeasure'],300),
    jitter(metric[,'ActivityLocation.LatitudeMeasure'],300),
           pch=19,
           col=cols)

    
    
    
    
    
    
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



    
    unique(results[,'ResultAnalyticalMethod.MethodName'])
    table(anc[,'ResultSampleFractionText'])
head(ph)
    
    
               head(nitrate)
    dim(NWQ.crop)
    
    
    
    date <- ymd(NWQ[,'ActivityStartDate'])
    
    
    water.buff <- st_intersection(water.projected, buff)
        
    plot(water.projected[1])
    plot(water.buff[1])
    u.water.istes
    