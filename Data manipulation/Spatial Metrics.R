
library(SSN)
library(tidyverse)
library(sf)
library(riverdist)
library(raster)
library(elevatr)
library(nhdplusTools)
library(nhdR)

# section coordinates, mid-points
coords <- read_tsv("YOY_growth_sections_coords.txt") %>% 
  dplyr::select(name, latitude, longitude) 
coords <- coords %>% separate(name, into = c("stream", "section"), sep = "(?<=[A-Za-z])(?=[0-9])")
coords <- coords %>% separate(section, into = c("section", "ref"), sep = 1)
coords <- coords %>% pivot_wider(id_cols = c(stream, section), names_from = ref, values_from = c("latitude", "longitude"))
coords <- coords %>% mutate(lat = rowMeans(dplyr::select(coords, latitude_A, latitude_B)), lon = rowMeans(dplyr::select(coords, longitude_A, longitude_B)))
coords2 <- coords %>% dplyr::select(stream, section, lat, lon)

coordsp <- SpatialPointsDataFrame(coords2[,c('lon','lat')], data = coords2, proj4string = crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')) # create spatial points df
coordsp <- spTransform(coordsp, crs('+proj=utm +zone=12 +ellps=WGS84')) # reproject to UTM

# import stream network with springs
snake <- importSSN("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Spatial/lsn6.ssn", predpts = NULL)
snake <- as.SpatialLines(snake)

# get elevation
coordspel <- get_elev_point(coordsp)

# format data frmae
df <- as_tibble(data.frame(coordspel))
df <- df %>% filter(stream != "Cottonwood" ) %>% filter(stream != "Deadmans") %>% filter(stream != "Flat")
df <- df %>% mutate(stream = recode(stream, Upper = "Upper Bar BC", Lower = "Lower Bar BC", FlatUp = "Flat", ThreeChan = "3 Channel"))
df <- df %>% rename(elev = elevation, utm_e = coords.x1, utm_n = coords.x2)
df2 <- df %>% dplyr::select(stream, section, elev, lat, lon, utm_n, utm_e)

write_csv(df2, "YOY_Stream_SpatialMetrics.csv")



##################
# was attempting to calculate watershed area

# Direct download of NHD basin
nwissite <- list(featureSource = "nwissite", 
                 featureID = "USGS-13022500")
basin <- get_nldi_basin(nwissite)
basin.utm <- toUTM(basin) # reproject basin to UTM

plot(st_geometry(basin.utm))
plot(snake, add = T)

##### hoping to calculate watershed area...but not sure it's necessary????

# DEM for elevation, 
dem <- raster("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Spatial/DEM/Snake_Watershed_utm.tif")
dem.crop <- crop(dem, extent(basin.utm))
dem.mask <- mask(dem.crop, basin.utm)

plot(dem.crop)
lines(snake)

# view
plot(snake)
points(coordsp, col = "blue")

coordsp_riv <- xy2convert(x = coordsp$lon)


