
library(tidyverse)
library(lubridate)
library(RColorBrewer)

getwd()
setwd("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth")

# load density data to get min and max sampling dates per stream/year
mmdates <- read_csv("YOY_CatchDensity_2021-2022.csv") %>% 
  mutate(year = year(date)) %>%
  group_by(stream, year) %>%
  summarize(mindate = min(date), maxdate = max(date))

# load relevant temp data and filter
t21 <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Environmental/Hydrology/3_Daily temp data/Daily temp_1020-1121.csv")
t22 <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Environmental/Hydrology/3_Daily temp data/Daily temp_1121-1122.csv")
temp <- rbind(t21, t22)
temp <- temp %>% filter(stream == "blackrock" | stream == "blacktail" | stream == "cliff" | stream == "crystal" | stream == "fall" | stream == "fish" | 
                        stream == "flat" | stream == "lowerbarbc" | stream == "mosquito" | stream == "spread" | stream == "threechannel" | stream == "upperbarbc" | stream == "willow")
temp <- temp %>% mutate(year = year(date))

# visualize temp and GDD trajectories for all streams and sites
streams <- unique(temp$stream)
years <- 2021:2022
cols <- c(brewer.pal(9, "Set1"), "black")

for(i in 1:length(streams)) {
  par(mfrow = c(2,2))
  for(j in 1:length(years)) {
    d <- temp %>% filter(stream == streams[i] & year == years[j])
    sites <- unique(d$site)
    plot(meanT ~ date, d, col = NA, main = paste(streams[i], years[j], sep = " "), 
         xlim = c(date(paste(years[j], "-08-01", sep = "")), date(paste(years[j], "-11-10", sep = ""))), ylim = c(0,20),
         xlab = "", ylab = "Mean Daily Temp")
    legend("bottomleft", legend = sites, lwd = 2, col = cols[1:length(sites)])
    for(k in 1:length(sites)) {
      d1 <- d %>% filter(site == sites[k] & date >= paste(years[j], "-08-01", sep = ""))
      lines(meanT ~ date, d1, lwd = 2, col = cols[k])
      }
    plot(seq(0, 1200, length.out = 100) ~ seq.Date(date(paste(years[j], "-08-01", sep = "")), date(paste(years[j], "-11-10", sep = "")), length.out = 100), col = NA, main = paste(streams[i], years[j], sep = " "),
         xlab = "", ylab = "Cumulative Growing Degree Days")
    for(k in 1:length(sites)) {
      d1 <- d %>% filter(site == sites[k] & date >= paste(years[j], "-08-01", sep = ""))
      d1$T0 <- ifelse(d1$meanT < 0, 0, d1$meanT)
      d1$GDD <- cumsum(d1$T0)
      lines(GDD ~ date, d1, lwd = 2, col = cols[k])
      }}}


# Simplifying assumption of thermal homogeneity within streams, no longitudinal change
# Except for multi-thread, branching streams, for which we have multple loggers
# BUT, preliminary inspection of these plots suggests very little difference within streams. Thus, 1x stream
temp2 <- temp %>% filter(site == "blackrock_1A" | 
                           site == "blacktail_1A" | 
                           # site == "blacktail_4B" | 
                           site == "cliff_1A" | 
                           site == "crystal_1A" | 
                           site == "fall_4B" | 
                           site == "fish_downstream" | 
                           site == "flat_4B" | 
                           site == "lowerbarbc_SG" | 
                           site == "mosquito_1A" | 
                           site == "spread_1A" | 
                           site == "threechannel_1A" | 
                           # site == "threechannel_4B" | 
                           site == "upperbarbc_SG" | 
                           # site == "upperbarbc_3_DSRightFk" |
                           site == "willow_4B")

temp3 <- temp2 %>% filter(date >= "2021-08-01" & date <= "2021-11-05")
temp4 <- temp2 %>% filter(date >= "2022-08-01" & date <= "2022-11-05")

# identify missing starting data
temp3 %>% group_by(site) %>% summarise(mindate = min(date))
temp4 %>% group_by(site) %>% summarise(mindate = min(date))

# fix missing data
temp3b <- temp3[temp3$site == "flat_4B",][rep(1,6),]
temp3b$date <- temp3b$date - c(6:1)
temp3 <- rbind(temp3, temp3b)

temp3b <- temp3[temp3$site == "fall_4B",][rep(1,2),]
temp3b$date <- temp3b$date - c(2:1)
temp3 <- rbind(temp3, temp3b)

temp4b <- temp4[temp4$site == "threechannel_1A",][rep(1,2),]
temp4b$date <- temp4b$date - c(2:1)
temp4 <- rbind(temp4, temp4b)

# order temp data frames
temp3 <- temp3 %>% group_by(site) %>% arrange(date)
temp4 <- temp4 %>% group_by(site) %>% arrange(date)

# bind fixed data
temp2 <- rbind(temp3, temp4)

#----------------------------------------------------------------#
# visualize temp and GDD trajectories for focal streams/sites
#----------------------------------------------------------------#


sites <- unique(temp3$site)
# cols <- hcl.colors(length(sites), "Spectral")
coldf <- read_csv("StreamColors.csv")
coldf$stream <- c("blackrock", "mosquito", "spread", "fall", "cliff", "willow", "crystal", "threechannel", "upperbarbc", "fish", "flat", "blacktail", "lowerbarbc")
temp3 <- temp3 %>% left_join(coldf)
temp4 <- temp4 %>% left_join(coldf)


png("Figures/TempData.png", units = "in", width = 8.5, height = 6, res = 1000)
par(mfrow = c(2,2), mar = c(2,4,1,1), mgp = c(2.3,0.8,0), oma = c(0,0,0,6))

# 2021, mean daily temp
plot(meanT ~ date, temp3, col = NA, xlim = c(date("2021-08-01"), date("2021-11-05")), xlab = "", ylab = "Mean Daily Temp")
# rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey80")
for(i in 1:length(sites)) {
  d <- temp3 %>% filter(site == sites[i])
  lines(meanT ~ date, d, lwd = 2, col = d$cols)
}

# 2021, cumulative degree days, cutoff = 0 deg
plot(seq(0, 1200, length.out = 100) ~ seq.Date(date("2021-08-01"), date("2021-11-05"), length.out = 100), 
     col = NA, xlab = "", ylab = "Cumulative GDD (0)")
# rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey80")
for(k in 1:length(sites)) {
  d <- temp3 %>% filter(site == sites[k])
  d$T0 <- ifelse(d$meanT < 0, 0, d$meanT) # DD above 0
  d$GDD <- cumsum(d$T0)
  lines(GDD ~ date, d, lwd = 2, col = d$cols)
}

# # 2021, cumulative degree days, cutoff = 5 deg
# plot(seq(0, 700, length.out = 100) ~ seq.Date(date("2021-08-01"), date("2021-11-05"), length.out = 100), 
#      col = NA, xlab = "", ylab = "Cumulative Growing Degree Days (5)")
# rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey80")
# for(k in 1:length(sites)) {
#   d <- temp3 %>% filter(site == sites[k])
#   d$T0 <- ifelse(d$meanT < 5, 0, d$meanT-5) # DD above 5Per Huntsman et al (2018), Coleman and Fausch (2007): DD calculated based on temp above which growth typically occurs in CT
#   d$GDD <- cumsum(d$T0)
#   lines(GDD ~ date, d, lwd = 2, col = cols[k])
# }

# 2022, mean daily temp
plot(meanT ~ date, temp4, col = NA, xlim = c(date("2022-08-01"), date("2022-11-05")), xlab = "", ylab = "Mean Daily Temp")
# rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey80")
for(i in 1:length(sites)) {
  d <- temp4 %>% filter(site == sites[i])
  lines(meanT ~ date, d, lwd = 2, col = d$cols)
}

# 2022, cumulative degree days, cutoff = 0 deg
plot(seq(0, 1200, length.out = 100) ~ seq.Date(date("2022-08-01"), date("2022-11-05"), length.out = 100), 
     col = NA, xlab = "", ylab = "Cumulative GDD (0)")
# rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey80")
for(k in 1:length(sites)) {
  d <- temp4 %>% filter(site == sites[k])
  d$T0 <- ifelse(d$meanT < 0, 0, d$meanT) # DD above 0
  d$GDD <- cumsum(d$T0)
  lines(GDD ~ date, d, lwd = 2, col = d$cols)
}

# # 2022, cumulative degree days, cutoff = 5 deg
# plot(seq(0, 700, length.out = 100) ~ seq.Date(date("2022-08-01"), date("2022-11-05"), length.out = 100), 
#      col = NA, xlab = "", ylab = "Cumulative Growing Degree Days (5)")
# rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey80")
# for(k in 1:length(sites)) {
#   d <- temp4 %>% filter(site == sites[k])
#   d$T0 <- ifelse(d$meanT < 5, 0, d$meanT-5) # DD above 5Per Huntsman et al (2018), Coleman and Fausch (2007): DD calculated based on temp above which growth typically occurs in CT
#   d$GDD <- cumsum(d$T0)
#   lines(GDD ~ date, d, lwd = 2, col = cols[k])
# }
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", legend = coldf$stream, fill = coldf$cols, bty = "n", xpd = NA)

dev.off()


# write site-year temperature data as separaate files for bioenergetics simulations
sites <- unique(temp3$site)
for (i in 1:length(sites)) {
  ddd <- temp3 %>% filter(site == sites[i]) %>% mutate(day = yday(date)) %>% dplyr::select(day, meanT) %>% rename(temperature = meanT)
  ddd <- ddd[,c(2:3)] %>% add_row(day = c(1,365), temperature = c(1,1)) %>% arrange(day)
  write_csv(ddd, paste("Bioenergetics/Temperature Files/Temp_", sites[i], "_2021.csv", sep = ""))
}
sites <- unique(temp4$site)
for (i in 1:length(sites)) {
  ddd <- temp4 %>% filter(site == sites[i]) %>% mutate(day = yday(date)) %>% dplyr::select(day, meanT) %>% rename(temperature = meanT)
  ddd <- ddd[,c(2:3)] %>% add_row(day = c(1,365), temperature = c(1,1)) %>% arrange(day)
  write_csv(ddd, paste("Bioenergetics/Temperature Files/Temp_", sites[i], "_2022.csv", sep = ""))
}


#----------------------------------------------------------------#
# calculate cumulative and inter-sampling period GDD, mean temp
#----------------------------------------------------------------#

# load template for all streams, years, dates, and reaches
dat <- read_csv("YOY_SpecificGrowthRate_Summary.csv") %>% 
  dplyr::select(stream, year, date, priordate, reach) %>% rename(section = reach)
dat2 <- read_csv("streamsectiontempsite.csv")
dat <- dat %>% left_join(dat2, by = c("stream", "section"))

# fudge sampling dates for Mosquito Creek in 2021 to enable calculation of temp metrics
dat$date[dat$stream == "Mosquito" & dat$year == 2021] <- c(rep("2021-08-08", times = 4), rep("2021-08-23", times = 4), rep("2021-09-13", times = 4), rep("2021-10-07", times = 4), rep("2021-11-05", times = 4))
dat$priordate[dat$stream == "Mosquito" & dat$year == 2021] <- c(rep(NA, times = 4), rep("2021-08-08", times = 4), rep("2021-08-23", times = 4), rep("2021-09-13", times = 4), rep("2021-10-07", times = 4))

# fudge sampling dates for Fish Creek in 2022 to enable calc of temp metrics for first sampling period, originally on Aug 1, 2022
dat$date[dat$stream == "Fish" & dat$date == "2022-08-01"] <- rep("2022-08-02", times = 4)
dat$priordate[dat$stream == "Fish" & dat$date == "2022-08-26"] <- rep("2022-08-02", times = 4)

# arbitrary emergence date of August 1
dat$priordate[is.na(dat$priordate)] <- date(paste(dat$year[is.na(dat$priordate)], "-08-01", sep = ""))

# degree days above 0/5 degrees C
temp2$DD0 <- ifelse(temp2$meanT < 0, 0, temp2$meanT) # DD above 0
temp2$DD5 <- ifelse(temp2$meanT < 5, 0, temp2$meanT-5) # DD above 0

dat$GDDcum0 <- NA
dat$GDDper0 <- NA
dat$GDDcum5 <- NA
dat$GDDper5 <- NA
dat$meantempper <- NA

sites <- unique(temp2$site)

for(i in 1:dim(dat)[1]) {
  d <- temp2 %>% filter(site == dat$tempsite[i] & date >= date(paste(dat$year[i], "-08-01", sep = "")) & date <= date(paste(dat$year[i], "-11-05", sep = "")))
  d1 <- d %>% filter(date > date(paste(dat$year[i], "-08-01", sep = "")) & date <= dat$date[i])
  dat$GDDcum0[i] <- sum(d1$DD0)
  dat$GDDcum5[i] <- sum(d1$DD5)
  d1 <- d %>% filter(date > dat$priordate[i] & date <= dat$date[i])
  dat$GDDper0[i] <- sum(d1$DD0)
  dat$GDDper5[i] <- sum(d1$DD5)
  dat$meantempper[i] <- mean(d1$meanT)
  print(i)
}
view(dat)
unique(is.na(dat))

# load percent groundwater estimates from sine-wave regression model and join to data table
gw <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Temperature/GroundwaterDynamics_Summary.csv") %>%
  dplyr::select(site, Amedian, Alower, Aupper, pmedian, plower, pupper, ASD, pSD)
dat <- dat %>% left_join(gw, by = c("tempsite" = "site"))


#----------------------------------------------------------------#
# write csv
#----------------------------------------------------------------#

write_csv(dat, "YOY_Stream_TempMetrics.csv")



#----------------------------------------------------------------#
# Groundwater effects on temperature and temperature sensitivity
#----------------------------------------------------------------#

# groundwater influence
gwmet <- read_csv("Groundwater Metrics/GroundwaterMetrics_Normalized_YOYsections.csv") %>% 
  separate(site, into = c("stream", "section"), sep = "_", remove = FALSE) %>% 
  separate(section, into = c("section", "updn"), sep = 1, remove = TRUE) %>% 
  mutate(section = as.numeric(section)) %>%
  dplyr::select(stream, section, springprev_basinmean_norm) %>% 
  rename(gw = springprev_basinmean_norm) %>% 
  add_row(stream = c("dum1","dum2"), section = c(NA,NA), gw = c(0,1)) %>% arrange(gw)

# assign color based on gw influence: purple = snowmelt, yellow = spring-fed
vrPal <- colorRampPalette(viridis(12))
gwmet <- gwmet %>% group_by(stream) %>% summarize(gw = mean (gw)) %>% arrange(gw) %>% mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))])
gwmet <- gwmet %>% filter(!stream %in% c("dum1","dum2"))

# view barplot...does ordering of gw influence make sense?
par(mar = c(5,10,2,2))
barplot(gwmet$gw, col = gwmet$cols, names.arg = gwmet$stream, horiz = TRUE, las = 2)
dev.off()

# recode stream names
gwmet <- gwmet %>% mutate(stream = recode(stream, 
                                          "Crystal" = "crystal", 
                                          "Cliff" = "cliff",
                                          "Spread" = "spread",
                                          "Willow" = "willow",
                                          "Mosquito" = "mosquito",
                                          "Blackrock" = "blackrock", 
                                          "Fall" = "fall",
                                          "Flat" = "flat",
                                          "3 Channel" = "threechannel",
                                          "Fish" = "fish",
                                          "Blacktail" = "blacktail",
                                          "Upper Bar BC" = "upperbarbc",
                                          "Lower Bar BC" = "lowerbarbc")) 

# join colors and gw to temp 
temp3 <- temp3 %>% left_join(gwmet)
temp4 <- temp4 %>% left_join(gwmet)

# temp sensitivity metrics
tempsensdf <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Temperature/TempSensitivityWithYearlySiteMetrics.csv") %>%
  select(site, year, slo.est) %>% rename(ts = slo.est)

# join to gw
temp3ts <- temp3 %>% group_by(site, year(date)) %>% summarize(gw = unique(gw), cols = unique(cols)) %>% rename(year = "year(date)") %>% left_join(tempsensdf)
temp4ts <- temp4 %>% group_by(site, year(date)) %>% summarize(gw = unique(gw), cols = unique(cols)) %>% rename(year = "year(date)") %>% left_join(tempsensdf)
tsdf <- rbind(temp3ts, temp4ts) %>% mutate(yrid = ifelse(year == 2021, 21, 24))



###### PLOT
streams <- unique(gwmet$stream)
gradpal <- tibble(gw = seq(from = 0, to = 1, length.out = 100)) %>%
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>% 
  filter(gw >= min(tsdf$gw) & gw <= max(tsdf$gw))

jpeg(filename = "Temperature/TempTimeSeriesAndTempSensitivityByGroundwater.jpg", height = 2.25, width = 7, units = "in", res = 1500)
par(mfrow = c(1,3), mar = c(3.5,3.5,0.5,0.5), mgp = c(2.2,0.8,0))

plot(meanT ~ date, temp3, col = NA, xlim = c(date("2021-08-01"), date("2021-11-05")), xlab = "", ylab = "Mean Daily Temp")
for (i in 1:length(streams)) {
  d <- temp3 %>% filter(stream == streams[i])
  lines(meanT ~ date, d, lwd = 2, col = d$cols)
}
legend("topright", legend = 2021, pch = NA, bty = "n")
gradientLegend(valRange = c(min(gradpal$gw), max(gradpal$gw)), color = gradpal$cols, length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.05, 0.05, 0.12, 0.4), dec = 2)

plot(meanT ~ date, temp4, col = NA, xlim = c(date("2022-08-01"), date("2022-11-05")), xlab = "", ylab = "Mean Daily Temp")
for (i in 1:length(streams)) {
  d <- temp4 %>% filter(stream == streams[i])
  lines(meanT ~ date, d, lwd = 2, col = d$cols)
}
legend("topright", legend = 2022, pch = NA, bty = "n")

plot(log(ts) ~ log(gw), tsdf, pch = NA, xlab = "log(Groundwater Influence)", ylab = "log(Temperature Sensitivity)")
tsmod <- lm(log(ts) ~ log(gw), tsdf)
gwvec <- seq(from = min(tsdf$gw), to = max(tsdf$gw), length.out = 100)
preds <- predict(tsmod, newdata = list(gw = gwvec), interval = "confidence")
polygon(x = c(log(gwvec), rev(log(gwvec))), y = c(preds[,2], rev(preds[,3])), col = alpha("black", 0.2), lty = 0)
lines(preds[,1] ~ log(gwvec), lwd = 2)
points(log(ts) ~ log(gw), tsdf, pch = yrid, bg = tsdf$cols)
legend("topright", legend = c(2021,2022), pch = c(21,24), bty = "n")

dev.off()






