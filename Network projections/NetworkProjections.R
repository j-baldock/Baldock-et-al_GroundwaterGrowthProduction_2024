#==========================================================================================================================================#
# PROJECT NAME: Groundwater structures fish growth and production across a riverscape
# CREATOR: Jeffrey R. Baldock, WY Coop Unit, University of Wyoming, jbaldock@uwyo.edu 
# PURPOSE: Project temporal dynamics of growth and production across the river network
# NOTES: 
#==========================================================================================================================================#

library(tidyverse)
library(mapview)
library(sf)
library(terra)
library(viridis)
library(plotfunctions)
library(whitebox)
library(tidyterra)
library(gganimate)


#-----------------------------------------------------------#
# Some data
#-----------------------------------------------------------#

# prediction points with groundwater index (springprev_iew05km_norm)
# gwmet <- vect("Groundwater Metrics/GroundwaterMetrics_Normalized_PredPoints.shp")
# gwmet$gwinf <- as.numeric(round(gwmet$springpr_1, digits = 2))
# mapview(st_as_sf(gwmet), zcol = "gwinf")
# sort(unique(gwmet$gwinf))
# length(unique(gwmet$gwinf))

# extract row numbers for relevant days
doys <- c(232:309)
which(doys == 244) # Sept 1
which(doys == 274) # Oct 1
which(doys == 305) # Nov 1

# basin, lakes, flowline, and slope
basin <- vect("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Spatial/Basin Delineation/BasinDelineation/MajorBasins_Watersheds.shp") 
basin <- subset(basin, basin$site == "UpperSnake")
lakes <- vect("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Spatial/Lakes/UpperSnake_Lakes.shp")
flowline <- vect("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/Watershed Delineation/Flowline/UpperSnake_Flowline.shp") 
slope <- rast("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Groundwater Seeps/Covariates/Snake_Slope.tif") #%>% intersect(flowline)
sloflo <- terra::mask(slope, flowline)

plot(slope)
hist(slope)
hist(sloflo)
mapview(raster::raster(slope))

#-----------------------------------------------------------#
# calculate spatial metrics for which to filter pred points
#-----------------------------------------------------------#
# 
# # buffer points and calculate mean stream gradient for each prediction point
# ptsbuf <- buffer(gwmet, 150)
# gwmet$slope <- NA
# for (i in 1:dim(gwmet)[1]) {
#   gwmet$slope[i] <- extract(sloflo, ptsbuf[i,], fun = mean, na.rm = TRUE)[,2]
#   print(i)
# }
# hist(gwmet$slope)
# mapview(st_as_sf(gwmet), zcol = "slope")
# # Kruse et al (1997) did not find any YSC at sites with stream gradients >10%
# mapview(st_as_sf(gwmet) %>% filter(slope <= 15), zcol = "slope")
# mapview(st_as_sf(gwmet) %>% filter(slope <= 10), zcol = "slope")
# # dim(gwmet %>% filter(slope <= 15))
# # dim(subset(gwmet, slope <= 10))

# 
# # stream order raster
# dem <- rast("Watershed Delineation/UpperSnake_DEM.tif")
# plot(dem)
# wbt_fill_burn(dem = "Watershed Delineation/UpperSnake_DEM.tif",
#               streams = "Watershed Delineation/Flowline/UpperSnake_Flowline.shp",
#               output = "Network Projections/Working/Snake_DEM_burn.tif")
# wbt_breach_depressions_least_cost(dem = "Watershed Delineation/UpperSnake_DEM.tif", 
#                                   output = "Network Projections/Working/Snake_DEM_breached.tif", dist = 5, fill = TRUE)
# wbt_fill_depressions_wang_and_liu(dem = "Network Projections/Working/Snake_DEM_breached.tif", 
#                                   output = "Network Projections/Working/Snake_DEM_filled_breached.tif")
# wbt_d8_pointer(dem = "Network Projections/Working/Snake_DEM_filled_breached.tif", 
#                output = "Network Projections/Working/Snake_mask_flowdir_D8.tif")
# wbt_repair_stream_vector_topology(input = "Watershed Delineation/Flowline/UpperSnake_Flowline.shp",
#                                   output = "Network Projections/Working/UpperSnake_Flowline_repair.shp", dist = 100)
# wbt_rasterize_streams(streams = "Network Projections/Working/UpperSnake_Flowline_repair.shp",
#                       base = "Watershed Delineation/UpperSnake_DEM.tif",
#                       output = "Network Projections/Working/Snake_flowline_raster.tif")
# wbt_remove_short_streams(d8_pntr = "Network Projections/Working/Snake_mask_flowdir_D8.tif",
#                          streams = "Network Projections/Working/Snake_flowline_raster.tif",
#                          output = "Network Projections/Working/Snake_flowline_raster_clean.tif", min_length = 50)
# wbt_strahler_stream_order(d8_pntr = "Network Projections/Working/Snake_mask_flowdir_D8.tif",
#                           streams = "Network Projections/Working/Snake_flowline_raster.tif",
#                           output = "Network Projections/Working/Snake_StrahlerOrder.tif")
# wbt_raster_streams_to_vector(streams = "Network Projections/Working/Snake_StrahlerOrder.tif",
#                              d8_pntr = "Network Projections/Working/Snake_mask_flowdir_D8.tif",
#                              output = "Network Projections/Working/Snake_flowline_strahler.shp")
# # flrast <- rast("Network Projections/Working/Snake_flowline_raster.tif")
# # plot(flrast, xlim = c(516100, 516400), ylim = c(4822500,4822700))
# # plot(flowline, add = T)
# storder <- vect("Network Projections/Working/Snake_flowline_strahler.shp")
# mapview(st_as_sf(storder), zcol = "STRM_VAL")
# mapview(st_as_sf(flowline))
# 
# # snap points
# wbt_jenson_snap_pour_points(pour_pts = "Groundwater Metrics/GroundwaterMetrics_Normalized_PredPoints.shp",
#                             streams = "Watershed Delineation/Working/Snake_mask_flowline_raster.tif",
#                             output = "Network Projections/Working/GroundwaterMetrics_Normalized_PredPoints_snapped.shp",
#                             snap_dist = 300) 

# load prediction points shapefile
gwmet <- vect("Network Projections/Working/GroundwaterMetrics_Normalized_PredPoints_snapped.shp") #%>% st_set_crs(st_crs(basin)) %>% st_intersection(basin)
gwmet$gwinf <- as.numeric(round(gwmet$springpr_1, digits = 2))
# mapview(st_as_sf(gwmet), zcol = "gwinf")

# define basin area for each point
predsites <- vect("Watershed Delineation/Prediction Points/UpperSnake_PredPts_300m_nolake_Watersheds.shp")
gwmet <- gwmet %>% left_join(tibble(as.data.frame(predsites)))
# mapview(st_as_sf(gwmet), zcol = "areasqkm")
# mapview(st_as_sf(gwmet) %>% filter(areasqkm <= 500), zcol = "areasqkm") + mapview(st_as_sf(flowline))

# buffer points to calculate average slope
ptsbuf <- buffer(gwmet, 150)

# define slope for each point
gwmet$slope <- extract(sloflo, ptsbuf, fun = mean, na.rm = TRUE)[,2]
hist(gwmet$slope)
# mapview(st_as_sf(gwmet), zcol = "slope")


##-----------------------------------------------------------##
## Growth Capacity and Cumulative Production
##-----------------------------------------------------------##

# projected growth through time ~ gw
growth <- read_csv("Growth by Time/GrowthByTime_ProjectedLength_50perc.csv")
# add column for 0.50
colnames(growth) <- paste("gw", round(seq(from = 0, to = 1, length.out = 100), digits = 2), sep = "")
growth <- growth %>% mutate(gw0.50 = (gw0.49 + gw0.51)/2) 
growth <- growth %>% dplyr::select(order(colnames(growth)))
# view(growth)
# growth at Sept 1, Oct 1, Nov 1
growth2 <- t(growth[c(which(doys == 244), which(doys == 274), which(doys == 305)),])
growth2 <- as_tibble(growth2) %>% 
  setNames(c("LenSep1", "LenOct1", "LenNov1")) %>% 
  mutate(RelLenSep1 = (LenSep1 - min(growth2))/(max(growth2) - min(growth2)),
         RelLenOct1 = (LenOct1 - min(growth2))/(max(growth2) - min(growth2)),
         RelLenNov1 = (LenNov1 - min(growth2))/(max(growth2) - min(growth2)),
         gwinf = as.numeric(seq(from = 0, to = 1, length.out = 101)))
# growth at daily time step
growth3 <- t(growth)
growth3 <- as_tibble(growth3) %>% mutate(gwinf = as.numeric(seq(from = 0, to = 1, length.out = 101)))
growth3a <- (growth3[,-80] - min(growth3[,-80])) / (max(growth3[,-80]) - min(growth3[,-80]))
growth3a <- as_tibble(growth3a) %>% mutate(gwinf = as.numeric(seq(from = 0, to = 1, length.out = 101)))


# projected production through time ~ gw
produc <- read_csv("Production by Time/ProductionByTime_ProjectedCumulativeProduction_Median.csv")
# add column for 0.50
colnames(produc) <- paste("gw", round(seq(from = 0, to = 1, length.out = 100), digits = 2), sep = "")
produc <- produc %>% mutate(gw0.50 = (gw0.49 + gw0.51)/2) 
produc <- produc %>% dplyr::select(order(colnames(produc)))
produc[produc > 10] <- 10 # truncate production at 12
# production at Sept 1, Oct 1, Nov 1
produc2 <- t(produc[c(which(doys == 244), which(doys == 274), which(doys == 305)),])
produc2 <- as_tibble(produc2) %>% 
  setNames(c("ProSep1", "ProOct1", "ProNov1")) %>% 
  mutate(RelProSep1 = ((ProSep1) - min(produc2))/(max(produc2) - min(produc2)),
         RelProOct1 = ((ProOct1) - min(produc2))/(max(produc2) - min(produc2)),
         RelProNov1 = ((ProNov1) - min(produc2))/(max(produc2) - min(produc2)),
         gwinf = as.numeric(seq(from = 0, to = 1, length.out = 101)))
# production at a daily time step
produc3 <- t(produc)
produc3 <- as_tibble(produc3) %>% mutate(gwinf = as.numeric(seq(from = 0, to = 1, length.out = 101)))
produc3a <- (produc3[,-79] - min(produc3[,-79])) / (max(produc3[,-79]) - min(produc3[,-79]))
produc3a <- as_tibble(produc3a) %>% mutate(gwinf = as.numeric(seq(from = 0, to = 1, length.out = 101)))

# Join gw metrics with growth and production by gwinf
trim.trailing <- function (x) sub("\\s+$", "", x) 
gwmet$gwinf <- as.numeric(trim.trailing(gwmet$gwinf))
growth2$gwinf <- as.numeric(trim.trailing(growth2$gwinf))
growth3$gwinf <- as.numeric(trim.trailing(growth3$gwinf))
produc2$gwinf <- as.numeric(trim.trailing(produc2$gwinf))
produc3$gwinf <- as.numeric(trim.trailing(produc3$gwinf))
growth3a$gwinf <- as.numeric(trim.trailing(growth3a$gwinf))
produc3a$gwinf <- as.numeric(trim.trailing(produc3a$gwinf))
gwmet2 <- gwmet %>% left_join(growth2) %>% left_join(produc2)
gwmet3gro <- gwmet %>% left_join(growth3) 
gwmet3pro <- gwmet %>% left_join(produc3)
gwmet3agro <- gwmet %>% left_join(growth3a) 
gwmet3apro <- gwmet %>% left_join(produc3a)


##-----------------------------------------------------------##
## Plot Discrete Time Steps
##-----------------------------------------------------------##

# plot network projections
jpeg("Network Projections/NetworkProjections_BodySize_CumulativeProduction.jpg", width = 5.6, height = 6, units = "in", res = 1500)
par(mfrow = c(2,3), mar = c(0,0,0,0), oma = c(1,1,1,1))
gwmet2b <- gwmet2 %>% filter(slope <= 10 & areasqkm <= 500) # filter prediction points by stream slope

# colors
# vrPal <- colorRampPalette(rev(inferno(12)))
vrPal <- colorRampPalette(rev(hcl.colors(12, "Rocket")))
# vrPal <- colorRampPalette(c("#FCAD12FF", "#F78311FF", "#E65D2FFF", "#CB4149FF", "#A92E5EFF", "#85216BFF", "#60136EFF", "#3A0963FF", "#140B35FF", "#000004FF"))
bgcol <- "grey80"
stcol <- "grey55"

# Growth - Sept 1
gwmet2b$cols <- vrPal(100)[as.numeric(cut(c(0, 1, gwmet2b$RelLenSep1), breaks = 100))][-c(1,2)]
# gwmet3$lwds <- (gwmet3$RelLenSep1+0.3)/2
# gwmet3$lwds[gwmet3$lwds > 0.4] <- 0.4
plot(st_geometry(st_as_sf(basin)), col = bgcol, lwd = 0.6, xaxs = "i", yaxs = "i")
plot((flowline), col = stcol, lwd = 0.6, add = T)
plot((gwmet2b), pch = 16, cex = 0.25, col = gwmet2b$cols, add = T)
plot((lakes), border = NA, col = stcol, add = T, lwd = 0.6)
plot((basin), lwd = 0.6, add = T)
# box()
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.15, 0.975, labels = "(a)", cex = 1.25)
text(0.5, 1, labels = "Sept", xpd = NA, cex = 1.5)
text(0, 0.5, labels = "Growth", xpd = NA, cex = 1.5, srt = 90)
par(usr = usr)

# Growth - Oct 1
gwmet2b$cols <- vrPal(100)[as.numeric(cut(c(0, 1, gwmet2b$RelLenOct1), breaks = 100))][-c(1,2)]
# gwmet3$lwds <- (gwmet3$RelLenOct1+0.3)/2
# gwmet3$lwds[gwmet3$lwds > 0.4] <- 0.4
plot(st_geometry(st_as_sf(basin)), col = bgcol, lwd = 0.6, xaxs = "i", yaxs = "i")
plot((flowline), col = stcol, lwd = 0.6, add = T)
plot((gwmet2b), pch = 16, cex = 0.25, col = gwmet2b$cols, add = T)
plot((lakes), border = NA, col = stcol, add = T, lwd = 0.6)
plot((basin), lwd = 0.6, add = T)
# box()
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.15, 0.975, labels = "(b)", cex = 1.25)
text(0.5, 1, labels = "Oct", xpd = NA, cex = 1.5)
par(usr = usr)

# Growth - Nov 1
gwmet2b$cols <- vrPal(100)[as.numeric(cut(c(0, 1, gwmet2b$RelLenNov1), breaks = 100))][-c(1,2)]
# gwmet3$lwds <- (gwmet3$RelLenNov1+0.3)/2
# gwmet3$lwds[gwmet3$lwds > 0.4] <- 0.4
plot(st_geometry(st_as_sf(basin)), col = bgcol, lwd = 0.6, xaxs = "i", yaxs = "i")
plot((flowline), col = stcol, lwd = 0.6, add = T)
plot((gwmet2b), pch = 16, cex = 0.25, col = gwmet2b$cols, add = T)
plot((lakes), border = NA, col = stcol, add = T, lwd = 0.6)
plot((basin), lwd = 0.6, add = T)
# box()
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.15, 0.975, labels = "(c)", cex = 1.25)
gradientLegend(valRange = c(min(growth2[,c(1:3)]), max(growth2[,c(1:3)])), color = vrPal(100), length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.72, 0.7, 0.78, 0.92), dec = 2)
text(0.82, 0.97, "Length", xpd = NA, cex = 0.85)
text(0.5, 1, labels = "Nov", xpd = NA, cex = 1.5)
par(usr = usr)

# Production - Sept 1
gwmet2b$cols <- vrPal(100)[as.numeric(cut(c(0, 1, gwmet2b$RelProSep1), breaks = 100))][-c(1,2)]
gwmet2b$lwds <- (gwmet2b$RelProSep1+0.3)/2.2
gwmet2b$lwds[gwmet2b$lwds > 0.55] <- 0.55
plot(st_geometry(st_as_sf(basin)), col = bgcol, lwd = 0.6, xaxs = "i", yaxs = "i")
plot((flowline), col = stcol, lwd = 0.6, add = T)
plot((gwmet2b), pch = 16, cex = gwmet2b$lwds, col = gwmet2b$cols, add = T)
plot((lakes), border = NA, col = stcol, add = T, lwd = 0.6)
plot((basin), lwd = 0.6, add = T)
# box()
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.15, 0.975, labels = "(d)", cex = 1.25)
text(0, 0.5, labels = "Production", xpd = NA, cex = 1.5, srt = 90)
par(usr = usr)

# Production - Oct 1
gwmet2b$cols <- vrPal(100)[as.numeric(cut(c(0, 1, gwmet2b$RelProOct1), breaks = 100))][-c(1,2)]
gwmet2b$lwds <- (gwmet2b$RelProOct1+0.3)/2.2
gwmet2b$lwds[gwmet2b$lwds > 0.55] <- 0.55
plot(st_geometry(st_as_sf(basin)), col = bgcol, lwd = 0.6, xaxs = "i", yaxs = "i")
plot((flowline), col = stcol, lwd = 0.6, add = T)
plot((gwmet2b), pch = 16, cex = gwmet2b$lwds, col = gwmet2b$cols, add = T)
plot((lakes), border = NA, col = stcol, add = T, lwd = 0.6)
plot((basin), lwd = 0.6, add = T)
# box()
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.15, 0.975, labels = "(e)", cex = 1.25)
par(usr = usr)

# Production - Nov 1
gwmet2b$cols <- vrPal(100)[as.numeric(cut(c(0, 1, gwmet2b$RelProNov1), breaks = 100))][-c(1,2)]
gwmet2b$lwds <- (gwmet2b$RelProNov1+0.3)/2.2
gwmet2b$lwds[gwmet2b$lwds > 0.55] <- 0.55
plot(st_geometry(st_as_sf(basin)), col = bgcol, lwd = 0.6, xaxs = "i", yaxs = "i")
plot((flowline), col = stcol, lwd = 0.6, add = T)
plot((gwmet2b), pch = 16, cex = gwmet2b$lwds, col = gwmet2b$cols, add = T)
plot((lakes), border = NA, col = stcol, add = T, lwd = 0.6)
plot((basin), lwd = 0.6, add = T)
# box()
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.15, 0.975, labels = "(f)", cex = 1.25)
gradientLegend(valRange = c(min(produc2[,c(1:3)]), max(produc2[,c(1:3)])), color = vrPal(100), length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.72, 0.7, 0.78, 0.92), dec = 2)
text(0.82, 0.97, "Production", xpd = NA, cex = 0.85)
par(usr = usr)

dev.off()


##-----------------------------------------------------------##
## Animated!
##-----------------------------------------------------------##

#### Growth

# manipulate the data to long form
gwmet3gro <- gwmet3gro %>% filter(slope <= 10 & areasqkm <= 500) 
ggtib <- as_tibble(gwmet3gro) %>% select(paste("V", 1:79, sep = "")) %>% mutate(x = geom(gwmet3gro)[,3], y = geom(gwmet3gro)[,4])
ggtib2 <- ggtib %>% gather(key = "day", value = "length", 1:79) %>% mutate(day = extract_numeric(day))
# ggtib2$cols <- vrPal(100)[as.numeric(cut(c(min(ggtib2$length), max(ggtib2$length), ggtib2$length), breaks = 100))][-c(1,2)]

# create the basemap
base_map <- ggplot() + 
  geom_spatvector(data = basin, fill = bgcol) + 
  geom_spatvector(data = flowline, col = stcol) + 
  geom_spatvector(data = lakes, col = NA, fill = stcol) +
  geom_spatvector(data = basin, fill = NA) + 
  theme_void()
png("Network Projections/NetworkProjections_BaseMap.png", units = "in", width = 4, height = 5, res = 500, bg = "transparent")
base_map
dev.off()

# animate the data
map_with_data <- base_map + geom_point(data = ggtib2, aes(x = x, y = y, color = length, group = day), size = 0.1) + scale_color_gradientn(colors = vrPal(100))
map_with_anim <- map_with_data + transition_time(day) + ggtitle("Day: {frame_time}")
animate(map_with_anim, nframes = 79, units = "in", width = 4, height = 5, res = 500)
anim_save("Network Projections/NetworkProjections_BodySize_Animation.gif")


#### Production

# manipulate the data to long form
gwmet3pro <- gwmet3pro %>% filter(slope <= 10 & areasqkm <= 500) 
ggtibP <- as_tibble(gwmet3pro) %>% select(paste("V", 1:78, sep = "")) %>% mutate(x = geom(gwmet3pro)[,3], y = geom(gwmet3pro)[,4])
ggtibP2 <- ggtibP %>% gather(key = "day", value = "production", 1:78) %>% mutate(day = extract_numeric(day))

# animate the data
map_with_data <- base_map + geom_point(data = ggtibP2, aes(x = x, y = y, color = production, group = day), size = 0.1) + scale_color_gradientn(colors = vrPal(100))
map_with_anim <- map_with_data + transition_time(day) + ggtitle("Day: {frame_time}")
animate(map_with_anim, nframes = 78, units = "in", width = 4, height = 5, res = 500)
anim_save("Network Projections/NetworkProjections_Production_Animation.gif")

# just the final day
png("Network Projections/NetworkProjections_Production_Day78.png", units = "in", width = 4, height = 5, res = 500, bg = "transparent")
base_map + 
  geom_point(data = ggtibP2 %>% filter(day == 1), aes(x = x, y = y, color = production), size = 0.1) + 
  geom_point(data = ggtibP2 %>% filter(day == 78), aes(x = x, y = y, color = production), size = 0.1) + 
  scale_color_gradientn(colors = vrPal(100))
dev.off()

#----------------------------------------#
# What is the contribution of streams differing in their extent of GW input to basin-wide productivity relative to their availability?
#----------------------------------------#

# redo this b/c we capped production for plotting

# projected growth through time ~ gw
growth <- read_csv("Growth by Time/GrowthByTime_ProjectedLength_50perc.csv")
# add column for 0.50
colnames(growth) <- paste("gw", round(seq(from = 0, to = 1, length.out = 100), digits = 2), sep = "")
growth <- growth %>% mutate(gw0.50 = (gw0.49 + gw0.51)/2) 
growth <- growth %>% dplyr::select(order(colnames(growth)))
# view(growth)
# growth at Sept 1, Oct 1, Nov 1
growth2 <- t(growth[c(which(doys == 244), which(doys == 274), which(doys == 305)),])
growth2 <- as_tibble(growth2) %>% 
  setNames(c("LenSep1", "LenOct1", "LenNov1")) %>% 
  mutate(RelLenSep1 = (LenSep1 - min(growth2))/(max(growth2) - min(growth2)),
         RelLenOct1 = (LenOct1 - min(growth2))/(max(growth2) - min(growth2)),
         RelLenNov1 = (LenNov1 - min(growth2))/(max(growth2) - min(growth2)),
         gwinf = as.numeric(seq(from = 0, to = 1, length.out = 101)))

# projected production through time ~ gw
produc <- read_csv("Production by Time/ProductionByTime_ProjectedCumulativeProduction_Median.csv")
# add column for 0.50
colnames(produc) <- paste("gw", round(seq(from = 0, to = 1, length.out = 100), digits = 2), sep = "")
produc <- produc %>% mutate(gw0.50 = (gw0.49 + gw0.51)/2) 
produc <- produc %>% dplyr::select(order(colnames(produc)))
# production at Sept 1, Oct 1, Nov 1
produc2 <- t(produc[c(which(doys == 244), which(doys == 274), which(doys == 305)),])
produc2 <- as_tibble(produc2) %>% 
  setNames(c("ProSep1", "ProOct1", "ProNov1")) %>% 
  mutate(RelProSep1 = ((ProSep1) - min(produc2))/(max(produc2) - min(produc2)),
         RelProOct1 = ((ProOct1) - min(produc2))/(max(produc2) - min(produc2)),
         RelProNov1 = ((ProNov1) - min(produc2))/(max(produc2) - min(produc2)),
         gwinf = as.numeric(seq(from = 0, to = 1, length.out = 101)))

# Join gw metrics with growth and production by gwinf
trim.trailing <- function (x) sub("\\s+$", "", x) 
gwmet$gwinf <- trim.trailing(gwmet$gwinf)
growth2$gwinf <- trim.trailing(growth2$gwinf)
produc2$gwinf <- trim.trailing(produc2$gwinf)
gwmet2 <- gwmet %>% left_join(growth2) %>% left_join(produc2)


# covert to tibble and simplify gwinf to bins
gwmetdat <- as_tibble(gwmet2) %>% filter(slope <= 10 & areasqkm <= 500) %>% select(gwinf, LenSep1, LenOct1, LenNov1, RelLenSep1, RelLenOct1, RelLenNov1, ProSep1, ProOct1, ProNov1, RelProSep1, RelProOct1, RelProNov1) %>%
  mutate(gwinf = as.numeric(gwinf))
gwmetdat$gwinfbin <- ifelse(gwmetdat$gwinf >= 0 & gwmetdat$gwinf <= 0.1, 0.05,
                            ifelse(gwmetdat$gwinf > 0.1 & gwmetdat$gwinf <= 0.2, 0.15,
                                   ifelse(gwmetdat$gwinf > 0.2 & gwmetdat$gwinf <= 0.3, 0.25,
                                          ifelse(gwmetdat$gwinf > 0.3 & gwmetdat$gwinf <= 0.4, 0.35,
                                                 ifelse(gwmetdat$gwinf > 0.4 & gwmetdat$gwinf <= 0.5, 0.45,
                                                        ifelse(gwmetdat$gwinf > 0.5 & gwmetdat$gwinf <= 0.6, 0.55,
                                                               ifelse(gwmetdat$gwinf > 0.6 & gwmetdat$gwinf <= 0.7, 0.65,
                                                                      ifelse(gwmetdat$gwinf > 0.7 & gwmetdat$gwinf <= 0.8, 0.75,
                                                                             ifelse(gwmetdat$gwinf > 0.8 & gwmetdat$gwinf <= 0.9, 0.85,
                                                                                    ifelse(gwmetdat$gwinf > 0.9 & gwmetdat$gwinf <= 1.0, 0.95, 1
                                                                                    ))))))))))
gwmetdat$gwinfbincoarse <- ifelse(gwmetdat$gwinf >= 0.5, "high", "low")

# Frequency of groundwater influence - network is dominated by reaches with low levels of GW
hist(gwmetdat$gwinfbin, breaks = seq(from = 0, to = 1, by = 0.1), xlab = "Groundwater Influence", main = "")
coarsesum <- gwmetdat %>% group_by(gwinfbincoarse) %>% summarize(n = n()) %>% mutate(rel = n / dim(gwmetdat)[1])
barplot(coarsesum$n, xlab = "Groundwater Influence", ylab = "Frequency", names.arg = c("low", "high"))

# summarize relative abundance and relative contribution to total growth and production
gwmetdatsum <- gwmetdat %>% 
  group_by(gwinfbin) %>% 
  summarize(n = n(), cumprod = sum(ProNov1), cumgrow = sum(LenNov1)) %>%
  mutate(relabun = n / sum(n), 
         relcontrprod = cumprod / sum(cumprod),
         relcontrgrow = cumgrow / sum(cumgrow))
gwmetdatsum$reldiffprod <- (gwmetdatsum$relcontrprod - gwmetdatsum$relabun) / gwmetdatsum$relabun
gwmetdatsum$reldiffgrow <- (gwmetdatsum$relcontrgrow - gwmetdatsum$relabun) / gwmetdatsum$relabun
# plot
prodmat <- t(as.matrix(gwmetdatsum %>% select(relabun, relcontrprod)))
barplot(prodmat, beside = TRUE) # barplot

###### Chi square test
chisq.test(x = gwmetdatsum$cumprod, p = gwmetdatsum$relabun)


# sam but with coarse GW binning
gwmetdatsumcoarse <- gwmetdat %>% 
  group_by(gwinfbincoarse) %>% 
  summarize(n = n(), cumprod = sum(ProNov1), cumgrow = sum(LenNov1)) %>%
  mutate(relabun = n / sum(n), 
         relcontrprod = cumprod / sum(cumprod),
         relcontrgrow = cumgrow / sum(cumgrow))
gwmetdatsumcoarse$reldiffprod <- (gwmetdatsumcoarse$relcontrprod - gwmetdatsumcoarse$relabun) / gwmetdatsumcoarse$relabun
gwmetdatsumcoarse$reldiffgrow <- (gwmetdatsumcoarse$relcontrgrow - gwmetdatsumcoarse$relabun) / gwmetdatsumcoarse$relabun
# plot
barplot(t(as.matrix(gwmetdatsumcoarse %>% select(relabun, relcontrprod))), beside = TRUE) # barplot





gwmetdatsum$relcont <- gwmetdatsum$cumprod / gwmetdatsum$n
plot(reldiff ~ gwinfbin, gwmetdatsum) # relative difference between relative contribution and relative abundance 

plot(relabun ~ gwinfbin, gwmetdatsum)
plot(reldiffprod ~ gwinfbin, gwmetdatsum) # this is basically the same as cumprod ~ gw
plot(reldiffgrow ~ gwinfbin, gwmetdatsum) # this is basically the same as cumprod ~ gw


plot(reabun ~ gwinfbin, gwmetdatsum, lwd = 2)


#-----------------------------------------------------------#
# Growth Rates and Productivity Rates
#-----------------------------------------------------------#

# model mcmc samples
grmcmc <- read.csv("Growth by Time/GrowthByTime_Length_mcmcsamps.csv")
prmcmc <- read.csv("Production by Time/ProductionByTime_Length_mcmcsamps.csv")

# mean and sds for z-transformation
grzscore <- read_csv("Growth by Time/GrowthByTime_z-score_table.csv")
przscore <- read_csv("Production by Time/ProductionByTime_z-score_table.csv")

# doys and gw vectors
doys <- c(232:309)
gw <- seq(from = 0, to = 1, length.out = 101)
# z-score for predictions
zdoyg <- (doys - grzscore$means[2]) / grzscore$sds[2]
zgwg <- (gw - grzscore$means[1]) / grzscore$sds[1]
zdoyp <- (doys - przscore$means[2]) / przscore$sds[2]
zgwp <- (gw - przscore$means[1]) / przscore$sds[1]

# median growth rates ~ doy + gw
grmat <- matrix(NA, nrow = length(doys), ncol = length(gw))
for (i in 1:length(doys)) { grmat[i,] <- median(grmcmc[,"alpha.adj"]) + median(grmcmc[,"beta.1."])*zdoyg[i] + median(grmcmc[,"beta.2."])*(zdoyg[i]^2) + median(grmcmc[,"beta.3."])*(zgwg) + median(grmcmc[,"beta.4."])*(zgwg)*zdoyg[i] + median(grmcmc[,"beta.5."])*(zgwg)*(zdoyg[i]^2) }
grmat2 <- t(grmat[c(which(doys == 244), which(doys == 274), which(doys == 305)),])
grmat2 <- as_tibble(grmat2) %>% setNames(c("GRSep1", "GROct1", "GRNov1")) %>% mutate(gwinf = as.numeric(seq(from = 0, to = 1, length.out = 101)))

# median productivity rates ~ doy + gw
prmat <- matrix(NA, nrow = length(doys), ncol = length(gw))
for (i in 1:length(doys)) { prmat[i,] <- (median(prmcmc[,"alpha.adj"]) + median(prmcmc[,"beta.1."])*zdoyp[i] + median(prmcmc[,"beta.2."])*(zdoyp[i]^2) + median(prmcmc[,"beta.3."])*(zgwp) + median(prmcmc[,"beta.4."])*(zgwp)*zdoyp[i] + median(prmcmc[,"beta.5."])*(zgwp)*(zdoyp[i]^2))}
prmat2 <- t(prmat[c(which(doys == 244), which(doys == 274), which(doys == 305)),])
prmat2 <- as_tibble(prmat2) %>% setNames(c("PRSep1", "PROct1", "PRNov1")) %>% mutate(gwinf = as.numeric(seq(from = 0, to = 1, length.out = 101)))

# Join gw metrics with growth and production by gwinf
trim.trailing <- function (x) sub("\\s+$", "", x) 
# gwmet$gwinf <- trim.trailing(gwmet$gwinf)
grmat2$gwinf <- trim.trailing(grmat2$gwinf)
prmat2$gwinf <- trim.trailing(prmat2$gwinf)
gwmet3 <- gwmet %>% left_join(grmat2) %>% left_join(prmat2)


# plot network projections 
# scaled to min/max rates across all time steps
vrPal <- colorRampPalette((inferno(12)))

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(min(grmat2[,c(1:3)]), max(grmat2[,c(1:3)]), gwmet3$GRSep1), breaks = 100))][-c(1,2)]
plot(gwmet3["GRSep1"], pch = 16, cex = 0.5, col = gwmet3$cols)

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(min(grmat2[,c(1:3)]), max(grmat2[,c(1:3)]), gwmet3$GROct1), breaks = 100))][-c(1,2)]
plot(gwmet3["GROct1"], pch = 16, cex = 0.5, col = gwmet3$cols)

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(min(grmat2[,c(1:3)]), max(grmat2[,c(1:3)]), gwmet3$GRNov1), breaks = 100))][-c(1,2)]
plot(gwmet3["GRNov1"], pch = 16, cex = 0.5, col = gwmet3$cols)

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(min(prmat2[,c(1:3)]), max(prmat2[,c(1:3)]), gwmet3$PRSep1), breaks = 100))][-c(1,2)]
plot(gwmet3["PRSep1"], pch = 16, cex = 0.5, col = gwmet3$cols)

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(min(prmat2[,c(1:3)]), max(prmat2[,c(1:3)]), gwmet3$PROct1), breaks = 100))][-c(1,2)]
plot(gwmet3["PROct1"], pch = 16, cex = 0.5, col = gwmet3$cols)

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(min(prmat2[,c(1:3)]), max(prmat2[,c(1:3)]), gwmet3$PRNov1), breaks = 100))][-c(1,2)]
plot(gwmet3["PRNov1"], pch = 16, cex = 0.5, col = gwmet3$cols)


# plot network projections 
# scaled to min/max rates within each time step
vrPal <- colorRampPalette((inferno(12)))

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(gwmet3$GRSep1), breaks = 100))]
plot(gwmet3["GRSep1"], pch = 16, cex = 0.5, col = gwmet3$cols)

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(gwmet3$GROct1), breaks = 100))]
plot(gwmet3["GROct1"], pch = 16, cex = 0.5, col = gwmet3$cols)

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(gwmet3$GRNov1), breaks = 100))]
plot(gwmet3["GRNov1"], pch = 16, cex = 0.5, col = gwmet3$cols)

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(gwmet3$PRSep1), breaks = 100))]
plot(gwmet3["PRSep1"], pch = 16, cex = 0.5, col = gwmet3$cols)

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(gwmet3$PROct1), breaks = 100))]
plot(gwmet3["PROct1"], pch = 16, cex = 0.5, col = gwmet3$cols)

gwmet3$cols <- vrPal(100)[as.numeric(cut(c(gwmet3$PRNov1), breaks = 100))]
plot(gwmet3["PRNov1"], pch = 16, cex = 0.5, col = gwmet3$cols)


