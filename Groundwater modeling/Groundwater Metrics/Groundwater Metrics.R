#==========================================================================================================================================#
# PROJECT NAME: Groundwater structures fish growth and production across a riverscape
# CREATOR: Jeffrey R. Baldock, WY Coop Unit, University of Wyoming, jbaldock@uwyo.edu 
# PURPOSE: derive metrics of groundwater influence for YOY growth sites, temperature sites, and network prediction points
#==========================================================================================================================================#

library(tidyverse)
library(sf)
library(mapview)
library(terra)
library(GGally)
library(parallel)


#------------------------------------------------------------#
# Load spring prevalence rasters, sites, and watersheds
#------------------------------------------------------------#

# spring prevalence from Maxent: complete tiff and buffered/no lake tiff
spring_full <- rast("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Groundwater Seeps/SpringTIFFs/SpringPrev_UpperSnake_BedSurf.tif")

# lakes in Upper Snake
lakes <- vect("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Spatial/Lakes/UpperSnake_Lakes.shp") 

# basin without lakes
basin <- vect("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Spatial/Basin Delineation/BasinDelineation/MajorBasins_Watersheds.shp") 
basin <- subset(basin, basin$site == "UpperSnake")
basin_nolakes <- erase(basin, lakes)

# flowline
flowline <- vect("Watershed Delineation/Flowline/UpperSnake_Flowline.shp")
flowbuff <- flowline %>% buffer(width = 100)

# remove lakes and buffer 
spring_nolakes <- mask(spring_full, basin_nolakes)
plot(spring_nolakes)
spring_nolakes_buff100 <- mask(spring_nolakes, flowbuff)
plot(spring_nolakes_buff100)
writeRaster(spring_nolakes, "G:/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Groundwater Seeps/SpringTIFFs/SpringPrev_UpperSnake_BedSurf_nolakes.tif")
writeRaster(spring_nolakes_buff100, "G:/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Groundwater Seeps/SpringTIFFs/SpringPrev_UpperSnake_BedSurf_nolakes_flowbuff100.tif")
spring_nolakes <- rast("G:/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Groundwater Seeps/SpringTIFFs/SpringPrev_UpperSnake_BedSurf_nolakes.tif")
spring_nolakes_buff100 <- rast("G:/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Groundwater Seeps/SpringTIFFs/SpringPrev_UpperSnake_BedSurf_nolakes_flowbuff100.tif")

# yoy growth sites and prediction sites
sites_yoy <- vect("Growth Section Coords/YOY_GrowthSectionCoords_DownstreamOnly.shp") 
sites_yoy$site <- paste(sites_yoy$stream, "_", sites_yoy$section, sites_yoy$updown, sep = "")
sites_pred <- vect("Watershed Delineation/Prediction Points/UpperSnake_PredPts_300m_nolake.shp")
sites_pred$site <- 1:dim(sites_pred)[1]
sites_temp <- vect("Watershed Delineation/Logger Points/UwyoTempLogger_Points_snapped.shp")
sites_temp$sitename <- sites_temp$site
sites_temp <- aggregate(sites_temp, by = "sitename")
sites_temp <- sites_temp %>% sort("sitename")
sites_temp$site <- 1:dim(sites_temp)[1]


# watershed shapefiles
sheds_yoy <- read_sf(dsn = "Watershed Delineation/Section Points", layer = "YOY_GrowthSection_Watersheds")
sheds_pred <- read_sf(dsn = "Watershed Delineation/Prediction Points", layer = "UpperSnake_PredPts_300m_nolake_Watersheds")
sheds_temp <- read_sf(dsn = "Watershed Delineation/Logger Points", layer = "UwyoTempLogger_Watersheds")
sheds_temp$sitename <- sheds_temp$site
sheds_temp <- sheds_temp[order(sheds_temp$sitename),]
sheds_temp$site <- 1:dim(sheds_temp)[1]
# sheds_pred <- vect("Watershed Delineation/Prediction Points/UpperSnake_PredPts_300m_nolake_Watersheds.shp")

# visualize
mapview(list(st_as_sf(sites_yoy), sheds_yoy), col.regions = c("blue", "white"), col = c("blue", "blue"), alpha.regions = c(1,0.2))
plot(sites_pred, pch = 16, cex = 0.5)
mapview(st_as_sf(sites_pred))


#------------------------------------------------------------#
# Calculate groundwater metrics for each basin
#------------------------------------------------------------#

### YOY Growth Section Sites

# extract average and weighted spring prevalence for each basin
sites <- sheds_yoy$site
gwlist <- list()
st <- Sys.time()
for (i in 1:length(sites)) {
  spring_mask <- mask(crop(spring_nolakes_buff100, sheds_yoy[sheds_yoy$site == sites[i],]), sheds_yoy[sheds_yoy$site == sites[i],]) # crop and mask by basin
  dist_rast <- distance(spring_mask, sites_yoy[sites_yoy$site == sites[i],]) %>% mask(spring_mask) # calculate distance between each raster cell and site location
  gwlist[[i]] <- tibble(site = sites[i],
                        basin_areasqkm = sheds_yoy$aresqkm[i],
                        springprev_point = extract(spring_nolakes_buff100, sites_yoy[sites_yoy$site == sites[i],], na.rm = TRUE)[,2],
                        springprev_basinmean = as.numeric(global(spring_mask, "mean", na.rm = T)), # extract(spring_buff, sheds_yoy[sheds_yoy$site == sites[i],], fun = mean, na.rm = TRUE)[,2],
                        springprev_iew01km = as.numeric(global(spring_mask * (1 / exp(dist_rast/1000)), "sum", na.rm = T) / global(1 / exp(dist_rast/1000), "sum", na.rm = T)),
                        springprev_iew05km = as.numeric(global(spring_mask * (1 / exp(dist_rast/5000)), "sum", na.rm = T) / global(1 / exp(dist_rast/5000), "sum", na.rm = T))
                        )
  print(i)
}
Sys.time() - st
gwmetrics_yoy <- do.call(rbind, gwlist) # bind as tibble
view(gwmetrics_yoy)

# calculate metrics of groundwater influence
gwmetrics_yoy <- gwmetrics_yoy %>% mutate(springprev_basinmean_log = log(springprev_basinmean),
                                          springprev_iew01km_log = log(springprev_iew01km),
                                          springprev_iew05km_log = log(springprev_iew05km))
write_csv(gwmetrics_yoy, "Groundwater Metrics/GroundwaterMetrics_raw_YOYsections.csv")



### Network Prediction Sites (can't get parallel to work. ~8.5 hrs in basic for loop)

# extract average and weighted spring prevalence for each basin
sites <- sheds_pred$site
gwlist <- list()
st <- Sys.time()
for (i in 1:length(sites)) {
  spring_mask <- mask(crop(spring_nolakes_buff100, sheds_pred[sheds_pred$site == sites[i],]), sheds_pred[sheds_pred$site == sites[i],]) # mask by basin
  dist_rast <- distance(spring_mask, sites_pred[sites_pred$site == sites[i],]) %>% mask(spring_mask) # calculate distance between each raster cell and site location
  gwlist[[i]] <- tibble(site = sites[i],
                        basin_areasqkm = sheds_pred$areasqkm[i],
                        springprev_point = extract(spring_nolakes_buff100, sites_pred[sites_pred$site == sites[i],], na.rm = TRUE)[,2],
                        springprev_basinmean = as.numeric(global(spring_mask, "mean", na.rm = T)), # extract(spring_buff, sheds_yoy[sheds_yoy$site == sites[i],], fun = mean, na.rm = TRUE)[,2],
                        springprev_iew01km = as.numeric(global(spring_mask * (1 / exp(dist_rast/1000)), "sum", na.rm = T) / global(1 / exp(dist_rast/1000), "sum", na.rm = T)),
                        springprev_iew05km = as.numeric(global(spring_mask * (1 / exp(dist_rast/5000)), "sum", na.rm = T) / global(1 / exp(dist_rast/5000), "sum", na.rm = T))
  )
  tmpFiles(remove = TRUE)
  print(i)
}
et <- Sys.time()
et - st
gwmetrics_pred <- do.call(rbind, gwlist) # bind as tibble
gwmetrics_pred$basin_areasqkm <- sheds_pred$areasqkm

# calculate metrics of groundwater influence
gwmetrics_pred <- gwmetrics_pred %>% mutate(springprev_basinmean_log = log(springprev_basinmean),
                                            springprev_iew01km_log = log(springprev_iew01km),
                                            springprev_iew05km_log = log(springprev_iew05km))
# ggpairs(gwmetrics_pred[,c(7:14)])
write_csv(gwmetrics_pred, "Groundwater Metrics/GroundwaterMetrics_raw_PredPoints.csv")


### Temp Logger Sites

# extract average and weighted spring prevalence for each basin
sites <- sheds_temp$site
gwlist <- list()
st <- Sys.time()
for (i in 1:length(sites)) {
  spring_mask <- mask(crop(spring_nolakes_buff100, sheds_temp[sheds_temp$site == sites[i],]), sheds_temp[sheds_temp$site == sites[i],]) # crop and mask by basin
  dist_rast <- distance(spring_mask, sites_temp[sites_temp$site == sites[i],]) %>% mask(spring_mask) # calculate distance between each raster cell and site location
  gwlist[[i]] <- tibble(site = sites[i],
                        basin_areasqkm = sheds_temp$aresqkm[i],
                        springprev_point = extract(spring_nolakes_buff100, sites_temp[sites_temp$site == sites[i],], na.rm = TRUE)[,2],
                        springprev_basinmean = as.numeric(global(spring_mask, "mean", na.rm = T)), # extract(spring_buff, sheds_yoy[sheds_yoy$site == sites[i],], fun = mean, na.rm = TRUE)[,2],
                        springprev_iew01km = as.numeric(global(spring_mask * (1 / exp(dist_rast/1000)), "sum", na.rm = T) / global(1 / exp(dist_rast/1000), "sum", na.rm = T)),
                        springprev_iew05km = as.numeric(global(spring_mask * (1 / exp(dist_rast/5000)), "sum", na.rm = T) / global(1 / exp(dist_rast/5000), "sum", na.rm = T))
  )
  print(i)
}
Sys.time() - st
gwmetrics_temp <- do.call(rbind, gwlist) # bind as tibble
view(gwmetrics_temp)

# calculate metrics of groundwater influence
gwmetrics_temp <- gwmetrics_temp %>% 
  left_join(tibble(site = sheds_temp$site, sitename = sheds_temp$sitename)) %>% 
  mutate(springprev_basinmean_log = log(springprev_basinmean),
         springprev_iew01km_log = log(springprev_iew01km),
         springprev_iew05km_log = log(springprev_iew05km))
write_csv(gwmetrics_temp, "Groundwater Metrics/GroundwaterMetrics_raw_TempLoggers.csv")


#------------------------------------------------------------#
# Normalize values to 0-1 based on min/max of pred sites
#------------------------------------------------------------#

gwmetrics_yoy <- read_csv("Groundwater Metrics/GroundwaterMetrics_raw_YOYsections.csv")
gwmetrics_pred <- read_csv("Groundwater Metrics/GroundwaterMetrics_raw_PredPoints.csv")
gwmetrics_temp <- read_csv("Groundwater Metrics/GroundwaterMetrics_raw_TempLoggers.csv")

# reorder columns
gwmetrics_pred <- gwmetrics_pred[names(gwmetrics_yoy)]

# pairs plots
ggpairs(gwmetrics_yoy[,c(4:9)])
ggpairs(gwmetrics_temp[,c(4:6,8:10)])
plot(gwinf_iew04km_log ~ gwinf_basinmean_log, gwmetrics_pred)

############################################################################################
# density plots, how well do our yoy sections cover what is available?
ggplot() + 
  geom_density(aes(x = springprev_basinmean), fill = "red", alpha = 0.2, data = gwmetrics_pred, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_basinmean), fill = "blue", alpha = 0.2, data = gwmetrics_yoy, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_basinmean), fill = "green", alpha = 0.2, data = gwmetrics_temp, inherit.aes = FALSE)
ggplot() + 
  geom_density(aes(x = springprev_iew01km), fill = "red", alpha = 0.2, data = gwmetrics_pred, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_iew01km), fill = "blue", alpha = 0.2, data = gwmetrics_yoy, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_iew01km), fill = "green", alpha = 0.2, data = gwmetrics_temp, inherit.aes = FALSE)
ggplot() + 
  geom_density(aes(x = springprev_iew05km), fill = "red", alpha = 0.2, data = gwmetrics_pred, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_iew05km), fill = "blue", alpha = 0.2, data = gwmetrics_yoy, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_iew05km), fill = "green", alpha = 0.2, data = gwmetrics_temp, inherit.aes = FALSE)
# logged
ggplot() + 
  geom_density(aes(x = springprev_basinmean_log), fill = "red", alpha = 0.2, data = gwmetrics_pred, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_basinmean_log), fill = "blue", alpha = 0.2, data = gwmetrics_yoy, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_basinmean_log), fill = "green", alpha = 0.2, data = gwmetrics_temp, inherit.aes = FALSE)
ggplot() + 
  geom_density(aes(x = springprev_iew01km_log), fill = "red", alpha = 0.2, data = gwmetrics_pred, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_iew01km_log), fill = "blue", alpha = 0.2, data = gwmetrics_yoy, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_iew01km_log), fill = "green", alpha = 0.2, data = gwmetrics_temp, inherit.aes = FALSE)
ggplot() + 
  geom_density(aes(x = springprev_iew05km_log), fill = "red", alpha = 0.2, data = gwmetrics_pred, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_iew05km_log), fill = "blue", alpha = 0.2, data = gwmetrics_yoy, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_iew05km_log), fill = "green", alpha = 0.2, data = gwmetrics_temp, inherit.aes = FALSE)


############################################################################################
# drop rows for which watersheds are improperly defined and major springprev_basinmean outliers
test <- merge(sites_pred, gwmetrics_pred)
test1 <- st_as_sf(test) %>% filter(basin_areasqkm < 0.05)
mapview(list(st_as_sf(flowline), test1))
test_drop <- test[test$basin_areasqkm > 0.05,]

boxplot(test_drop$springprev_basinmean)
abline(h = 0.4, col = "red")
boxplot(test_drop$springprev_iew05km)
abline(h = 0.4, col = "red")

# 12 "outliers": 5 are located with the the Gros Ventre diversion to Flat Creek and 4 are at the very headwaters of Spring Creek
# both of these locations have a high degree of hydrologic alteration and thus we feel comfortable dropping these. 
# Of the remaining 3 points, 1 is within a spring channel to the lower GV and 2 are at the very headwaters of the little Snake River.
test2 <- st_as_sf(test_drop) %>% filter(springprev_basinmean > 0.4)
mapview(list(st_as_sf(flowline), test2))
test2 <- st_as_sf(test_drop) %>% filter(springprev_iew05km > 0.4)
mapview(list(st_as_sf(flowline), test2))

# drop from data frames
gwmetrics_pred_drop <- gwmetrics_pred[gwmetrics_pred$basin_areasqkm > 0.05,] # drop rows for which watersheds are improperly defined
gwmetrics_pred_drop <- gwmetrics_pred_drop[gwmetrics_pred_drop$springprev_basinmean < 0.4,] # drop 12 outliers that are throwing off normalization
sheds_pred_drop <- sheds_pred[sheds_pred$site %in% gwmetrics_pred_drop$site,]
sites_pred_drop <- sites_pred[sites_pred$site %in% gwmetrics_pred_drop$site,]


############################################################################################
# summarize min and max of relevant columns in pred points for normalization
predmin <- gwmetrics_pred_drop %>% select(springprev_point, springprev_basinmean, springprev_iew01km, springprev_iew05km, springprev_basinmean_log, springprev_iew01km_log, springprev_iew05km_log) %>% summarize(across(everything(), min))
predmax <- gwmetrics_pred_drop %>% select(springprev_point, springprev_basinmean, springprev_iew01km, springprev_iew05km, springprev_basinmean_log, springprev_iew01km_log, springprev_iew05km_log) %>% summarize(across(everything(), max))
# predmin$springprev_point <- min(gwmetrics_pred_drop$springprev_point, na.rm = TRUE)
# predmax$springprev_point <- max(gwmetrics_pred_drop$springprev_point, na.rm = TRUE)

# normalization
gwmetrics_yoy_norm <- gwmetrics_yoy %>% mutate(springprev_point_norm = (springprev_point - predmin$springprev_point) / (predmax$springprev_point - predmin$springprev_point),
                                               springprev_basinmean_norm = (springprev_basinmean - predmin$springprev_basinmean) / (predmax$springprev_basinmean - predmin$springprev_basinmean),
                                               springprev_iew01km_norm = (springprev_iew01km - predmin$springprev_iew01km) / (predmax$springprev_iew01km - predmin$springprev_iew01km),
                                               springprev_iew05km_norm = (springprev_iew05km - predmin$springprev_iew05km) / (predmax$springprev_iew05km - predmin$springprev_iew05km),
                                               springprev_basinmean_log_norm = (springprev_basinmean_log - predmin$springprev_basinmean_log) / (predmax$springprev_basinmean_log - predmin$springprev_basinmean_log),
                                               springprev_iew01km_log_norm = (springprev_iew01km_log - predmin$springprev_iew01km_log) / (predmax$springprev_iew01km_log - predmin$springprev_iew01km_log),
                                               springprev_iew05km_log_norm = (springprev_iew05km_log - predmin$springprev_iew05km_log) / (predmax$springprev_iew05km_log - predmin$springprev_iew05km_log))
gwmetrics_pred_norm <- gwmetrics_pred_drop %>% mutate(springprev_point_norm = (springprev_point - predmin$springprev_point) / (predmax$springprev_point - predmin$springprev_point),
                                                      springprev_basinmean_norm = (springprev_basinmean - predmin$springprev_basinmean) / (predmax$springprev_basinmean - predmin$springprev_basinmean),
                                                      springprev_iew01km_norm = (springprev_iew01km - predmin$springprev_iew01km) / (predmax$springprev_iew01km - predmin$springprev_iew01km),
                                                      springprev_iew05km_norm = (springprev_iew05km - predmin$springprev_iew05km) / (predmax$springprev_iew05km - predmin$springprev_iew05km),
                                                      springprev_basinmean_log_norm = (springprev_basinmean_log - predmin$springprev_basinmean_log) / (predmax$springprev_basinmean_log - predmin$springprev_basinmean_log),
                                                      springprev_iew01km_log_norm = (springprev_iew01km_log - predmin$springprev_iew01km_log) / (predmax$springprev_iew01km_log - predmin$springprev_iew01km_log),
                                                      springprev_iew05km_log_norm = (springprev_iew05km_log - predmin$springprev_iew05km_log) / (predmax$springprev_iew05km_log - predmin$springprev_iew05km_log))
gwmetrics_temp_norm <- gwmetrics_temp %>% mutate(springprev_point_norm = (springprev_point - predmin$springprev_point) / (predmax$springprev_point - predmin$springprev_point),
                                                 springprev_basinmean_norm = (springprev_basinmean - predmin$springprev_basinmean) / (predmax$springprev_basinmean - predmin$springprev_basinmean),
                                                 springprev_iew01km_norm = (springprev_iew01km - predmin$springprev_iew01km) / (predmax$springprev_iew01km - predmin$springprev_iew01km),
                                                 springprev_iew05km_norm = (springprev_iew05km - predmin$springprev_iew05km) / (predmax$springprev_iew05km - predmin$springprev_iew05km),
                                                 springprev_basinmean_log_norm = (springprev_basinmean_log - predmin$springprev_basinmean_log) / (predmax$springprev_basinmean_log - predmin$springprev_basinmean_log),
                                                 springprev_iew01km_log_norm = (springprev_iew01km_log - predmin$springprev_iew01km_log) / (predmax$springprev_iew01km_log - predmin$springprev_iew01km_log),
                                                 springprev_iew05km_log_norm = (springprev_iew05km_log - predmin$springprev_iew05km_log) / (predmax$springprev_iew05km_log - predmin$springprev_iew05km_log))


# check to ensure direct relationship
plot(springprev_basinmean_norm ~ springprev_basinmean, gwmetrics_yoy_norm)
plot(springprev_basinmean_norm ~ springprev_basinmean, gwmetrics_pred_norm)
plot(springprev_iew05km_log_norm ~ springprev_iew05km_log, gwmetrics_temp_norm)

# density plots, normalized, dropped problem sheds
ggplot() + 
  geom_density(aes(x = springprev_basinmean_norm), fill = "red", alpha = 0.2, data = gwmetrics_pred_norm, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_basinmean_norm), fill = "blue", alpha = 0.2, data = gwmetrics_yoy_norm, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_basinmean_norm), fill = "green", alpha = 0.2, data = gwmetrics_temp_norm, inherit.aes = FALSE)
ggplot() + 
  geom_density(aes(x = springprev_iew01km_norm), fill = "red", alpha = 0.2, data = gwmetrics_pred_norm, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_iew01km_norm), fill = "blue", alpha = 0.2, data = gwmetrics_yoy_norm, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_iew01km_norm), fill = "green", alpha = 0.2, data = gwmetrics_temp_norm, inherit.aes = FALSE)
ggplot() + 
  geom_density(aes(x = springprev_iew05km_norm), fill = "red", alpha = 0.2, data = gwmetrics_pred_norm, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_iew05km_norm), fill = "blue", alpha = 0.2, data = gwmetrics_yoy_norm, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_iew05km_norm), fill = "green", alpha = 0.2, data = gwmetrics_temp_norm, inherit.aes = FALSE)

ggplot() + 
  geom_density(aes(x = springprev_basinmean_log_norm), fill = "red", alpha = 0.2, data = gwmetrics_pred_norm, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_basinmean_log_norm), fill = "blue", alpha = 0.2, data = gwmetrics_yoy_norm, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_basinmean_log_norm), fill = "green", alpha = 0.2, data = gwmetrics_temp_norm, inherit.aes = FALSE)
ggplot() + 
  geom_density(aes(x = springprev_iew01km_log_norm), fill = "red", alpha = 0.2, data = gwmetrics_pred_norm, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_iew01km_log_norm), fill = "blue", alpha = 0.2, data = gwmetrics_yoy_norm, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_iew01km_log_norm), fill = "green", alpha = 0.2, data = gwmetrics_temp_norm, inherit.aes = FALSE)
ggplot() + 
  geom_density(aes(x = springprev_iew05km_log_norm), fill = "red", alpha = 0.2, data = gwmetrics_pred_norm, inherit.aes = FALSE) + 
  geom_density(aes(x = springprev_iew05km_log_norm), fill = "blue", alpha = 0.2, data = gwmetrics_yoy_norm, inherit.aes = FALSE) +
  geom_density(aes(x = springprev_iew05km_log_norm), fill = "green", alpha = 0.2, data = gwmetrics_temp_norm, inherit.aes = FALSE)



#------------------------------------------------------------#
# Join gwmetrics data frames to spat vectors of points
#------------------------------------------------------------#

# join 
sites_yoy_m <- merge(sites_yoy, gwmetrics_yoy_norm)
sites_pred_m <- merge(sites_pred_drop, gwmetrics_pred_norm)
sites_temp_m <- merge(sites_temp, gwmetrics_temp_norm)

# yoy sites
mapview(st_as_sf(sites_yoy_m), zcol = "springprev_basinmean_norm")
mapview(st_as_sf(sites_yoy_m), zcol = "springprev_iew01km_norm")
mapview(st_as_sf(sites_yoy_m), zcol = "springprev_iew05km_norm")
mapview(st_as_sf(sites_yoy_m), zcol = "springprev_iew01km_log_norm")
mapview(st_as_sf(sites_yoy_m), zcol = "springprev_iew05km_log_norm")

# temp sites
mapview(st_as_sf(sites_temp_m), zcol = "springprev_basinmean_norm")
mapview(st_as_sf(sites_temp_m), zcol = "springprev_iew01km_norm")
mapview(st_as_sf(sites_temp_m), zcol = "springprev_iew05km_norm")
mapview(st_as_sf(sites_temp_m), zcol = "springprev_iew01km_log_norm")
mapview(st_as_sf(sites_temp_m), zcol = "springprev_iew05km_log_norm")

# pred sites
mapview(st_as_sf(sites_pred_m), zcol = "springprev_basinmean_norm")
mapview(st_as_sf(sites_pred_m), zcol = "springprev_iew01km_norm")
mapview(st_as_sf(sites_pred_m), zcol = "springprev_iew05km_norm")
mapview(st_as_sf(sites_pred_m), zcol = "springprev_iew01km_log_norm")
mapview(st_as_sf(sites_pred_m), zcol = "springprev_iew05km_log_norm")

# plot to compare metrics
gwmetrics_yoy_norm %>% 
  select(site, springprev_basinmean_norm, springprev_iew01km_norm, springprev_iew05km_norm) %>% 
  gather("var", "value", -site) %>% ggplot(aes(x = var, y = value, group = site, col = site)) + geom_line()
gwmetrics_yoy_norm %>% 
  select(site, springprev_basinmean_log_norm, springprev_iew01km_log_norm, springprev_iew05km_log_norm) %>% 
  gather("var", "value", -site) %>% ggplot(aes(x = var, y = value, group = site, col = site)) + geom_line()




#------------------------------------------------------------#
# Write out CSVs and shapefiles
#------------------------------------------------------------#

write_csv(gwmetrics_yoy_norm, "Groundwater Metrics/GroundwaterMetrics_Normalized_YOYsections.csv")
write_csv(gwmetrics_pred_norm, "Groundwater Metrics/GroundwaterMetrics_Normalized_PredPoints.csv")
write_csv(gwmetrics_temp_norm, "Groundwater Metrics/GroundwaterMetrics_Normalized_TempPoints.csv")

writeVector(sites_yoy_m, "Groundwater Metrics/GroundwaterMetrics_Normalized_YOYsections.shp", overwrite = TRUE)
writeVector(sites_pred_m, "Groundwater Metrics/GroundwaterMetrics_Normalized_PredPoints.shp", overwrite = TRUE)
writeVector(sites_temp_m, "Groundwater Metrics/GroundwaterMetrics_Normalized_TempPoints.shp", overwrite = TRUE)



#------------------------------------------------------------#
# Create HTML file for easy viewing
#------------------------------------------------------------#

# bring in groundwater metrics
gwmet <- read_sf(dsn = "Groundwater Metrics", layer = "GroundwaterMetrics_Normalized_PredPoints") %>% mutate(GroundwaterIndex = as.numeric(round(springpr_1, digits = 3)))
gwmet2 <- gwmet[,c("site", "basin_are0", "GroundwaterIndex")]
gwmet3 <- gwmet2[gwmet2$basin_are0 <= 500,] # drop large rivers

# cut lakes from flowline
flowline_nolake <- erase(flowline, lakes)

# create map layers
m1 <- mapview(basin, color = "black", col.regions = "white", alpha.regions = 0, lwd = 3)
m2 <- mapview(flowline_nolake, color = "dodgerblue4", lwd = 2)
m3 <- mapview(lakes, color = "dodgerblue4", lwd = 2, col.regions = "lightblue", alpha.regions = 1)
m4 <- mapview(gwmet3, zcol = "GroundwaterIndex", lwd = 0.75, col.regions = colorRampPalette(rev(viridis(12))), alpha.regions = 1)

# save map output
mapshot(m1 + m2 + m3 + m4, url = "Groundwater Metrics/GroundwaterIndex_InteractiveMap.html")




# test <- tibble(site = gwmetrics2$site, 
#                rank1 = rank(gwmetrics2$gwinf_basinmean_log),
#                rank2 = rank(gwmetrics2$gwinf_iew01km_log),
#                rank3 = rank(gwmetrics2$gwinf_iew04km_log),
#                rank4 = rank(gwmetrics2$gwinf_iew15km_log),
#                rank5 = rank(gwmetrics2$gwinf_basinmean),
#                rank6 = rank(gwmetrics2$gwinf_iew01km),
#                rank7 = rank(gwmetrics2$gwinf_iew04km),
#                rank8 = rank(gwmetrics2$gwinf_iew15km))
# test <- arrange(test, rank1)
# plot(c(1,52) ~ c(1,8), type = "n")
# for(i in 1:dim(test)[1]) {
#   lines(as.numeric(test[i,c(2:9)]) ~ c(1:8), type = "b")
# }
# 
# rank(gwmetrics[,c(4:7,9:12)])
# 
# ggparcoord(gwmetrics2, columns = 4:7, groupColumn = 1)
# ggparcoord(gwmetrics2, columns = 9:12, groupColumn = 1)
# ggparcoord(gwmetrics2, columns = 13:16, groupColumn = 1)
# 
# global(spring_mask * weight01km, "sum", na.rm = T) / length(cells(spring_mask))
