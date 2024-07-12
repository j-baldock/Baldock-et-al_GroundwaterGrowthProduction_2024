#==========================================================================================================================================#
# PROJECT NAME: Groundwater structures fish growth and production across a riverscape
# CREATOR: Jeffrey R. Baldock, WY Coop Unit, University of Wyoming, jbaldock@uwyo.edu 
# PURPOSE: Merge data 
# NOTES: 
#==========================================================================================================================================#

library(tidyverse)
library(lubridate)
library(GGally)
library(AICcmodavg)
library(RColorBrewer)
library(irr)

#------------------------------------------------------#
# Load relevant datasets
#------------------------------------------------------#

# temperature
temp <- read_csv("YOY_Stream_TempMetrics.csv") %>% dplyr::select(stream, tempsite, year, date, priordate, section, meantempper, Amedian)

# redd counts
redds <- read_csv("YOY_ReddCounts_Summary.csv") %>% mutate(year = year(date)) %>% dplyr::select(stream, section, year, avg) %>% rename(avgredds = avg)

# groundwater metrics (including basin area)
gwmet <- read_csv("Groundwater Metrics/GroundwaterMetrics_Normalized_YOYsections.csv") %>% 
  separate(site, into = c("stream", "section"), sep = "_", remove = FALSE) %>% 
  separate(section, into = c("section", "updn"), sep = 1, remove = TRUE) %>% 
  mutate(section = as.numeric(section)) %>%
  dplyr::select(stream, section, basin_areasqkm, springprev_point_norm, springprev_basinmean_norm, springprev_iew01km_norm, springprev_iew05km_norm, springprev_iew01km_log_norm, springprev_iew05km_log_norm)
ggpairs(gwmet %>% dplyr::select(basin_areasqkm, springprev_point_norm, springprev_basinmean_norm, springprev_iew01km_norm, springprev_iew05km_norm, springprev_iew01km_log_norm, springprev_iew05km_log_norm))

# density
density <- read_csv("YOY_CatchDensity_2021-2022.csv") #%>% mutate(effden_permeter_log = log(effden_permeter))

# raw length/weight
lenwt <- read_csv("YOY_LengthWeight_Summary.csv") #%>% mutate(year = year(date))


#------------------------------------------------------#
# join covariate data
#------------------------------------------------------#

dat <- temp %>% 
  left_join(redds) %>% #, by = c("stream", "year", "section")) %>% #
  left_join(gwmet) %>% #, by = c("stream", "section")) %>%
  left_join(density) %>% #, by = c("stream", "date", "section")) %>%
  left_join(lenwt) #, by = c("stream", "date", "section")) %>% 
  )

# calculate prior catch, weight, length
ppdat <- dat %>% 
  dplyr::select(stream, year, date, section, catch, lenmean, wtmean) %>% 
  rename(priordate = "date", priorcatch = "catch", priorlenmean = "lenmean", priorwtmean = "wtmean")
dat <- dat %>% left_join(ppdat)

# generate misc variables
dat$stream <- as.factor(dat$stream)
dat$doy <- yday(dat$date)
dat$density_log <- log(dat$density + 0.01)
dat$cpue_log <- log(dat$cpue + 0.01)
dat$effden_permeter_log <- log(dat$effden_permeter + 0.01)
dat$avgredds_log <- log(dat$avgredds + 0.01)

# recode stream, section, year for jags
dat$strid <- as.numeric(dat$stream)
dat$yrid <- as.numeric(as.factor(dat$year))
dat$sectid <- as.factor(paste(dat$strid, dat$section, sep = "."))

# write out and read
write_csv(dat, "Growth_DataTable_WithCovariates.csv")

