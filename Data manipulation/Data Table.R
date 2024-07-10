# MERGE DATA and EXPLORE VISUALLY


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

# site level metrics (excluding basin area)
spatial <- read_csv("YOY_Stream_SpatialMetrics.csv") 

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

# growth
growth <- read_csv("YOY_SpecificGrowthRate_Summary.csv") %>% 
  dplyr::select(stream, date, reach, growthwtmean, growthwtsd, growthlenmean, growthlensd) %>%
  rename(section = reach)

# production
production <- read_csv("YOY_Production_Summary.csv") %>% rename(section = reach) %>% dplyr::select(stream, date, section, prodwtmean, prodwtsd, prodlenmean, prodlensd)

# temp sensitivity metrics
tempsensdf <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Temperature/TempSensitivityWithYearlySiteMetrics.csv") %>%
  group_by(site) %>% summarize(# gw_strbuff = unique(gw_strbuff), 
                               # gwindex = unique(gwindex),
                               # basin_areasqkm = unique(basin_areasqkm),
                               slope_basin = unique(slope_basin),
                               slope_stream = unique(slope_stream),
                               forest_strbuff = unique(forest_strbuff),
                               loglakearea = unique(loglakearea),
                               ts = mean(slo.est)) %>% ungroup()
ggpairs(tempsensdf %>% select(slope_basin, slope_stream, forest_strbuff, loglakearea, ts))

# consumption
consump <- read_csv("Bioenergetics/YOYGrowth_BioE_DesignFile_Means.csv")
consump2 <- read_csv("Bioenergetics/FB4_Log_File_BaldockYOYGrowthMeans.csv") %>% dplyr::select(Run_Name, 'p-value') %>% rename(pvalue = 'p-value')
consump <- consump %>% left_join(consump2) %>% dplyr::select(stream, date, section, year, pvalue)


#------------------------------------------------------#
# join covariate data
#------------------------------------------------------#

dat <- temp %>% 
  left_join(redds) %>% #, by = c("stream", "year", "section")) %>% #
  left_join(spatial) %>% #, by = c("stream", "section")) %>% 
  left_join(gwmet) %>% #, by = c("stream", "section")) %>%
  left_join(density) %>% #, by = c("stream", "date", "section")) %>%
  left_join(production) %>% #, by = c("stream", "date", "section")) %>% 
  left_join(lenwt) %>% #, by = c("stream", "date", "section")) %>% 
  left_join(growth) %>% #, by = c("stream", "date", "section")) %>%
  left_join(consump) %>% #, by = c("stream", "date", "section")) %>% 
  left_join(tempsensdf, by = c("tempsite" = "site")) # %>%
  # mutate(gwindex = gw_strbuff / log(basin_areasqkm),
  #        loggwindex = log(gwindex),
  #        relloggwindex = 1 - (loggwindex/min(loggwindex))
  #        loggwindex2_scaled = 1 - (abs(loggwindex2) / max(abs(loggwindex2)))
  #        loggwindex2_scaled = (abs(min(loggwindex2)) + loggwindex2) / max(abs(min(loggwindex2)) + loggwindex2)
  

# calculate prior catch, weight, length
ppdat <- dat %>% 
  dplyr::select(stream, year, date, section, catch, lenmean, wtmean) %>% 
  rename(priordate = "date", priorcatch = "catch", priorlenmean = "lenmean", priorwtmean = "wtmean")
dat <- dat %>% left_join(ppdat)

# calculate measures of growth based on length or weight: absolute, instantaneous, relative, relative (logged), specific
dat <- dat %>%
  mutate(agr_len = (lenmean - priorlenmean) / (yday(date) - yday(priordate)),
         igr_len = (log(lenmean) - log(priorlenmean)) / (yday(date) - yday(priordate)),
         rgr_len = (lenmean - priorlenmean) / ((yday(date) - yday(priordate))*priorlenmean),
         rgrl_len = (log(lenmean) - log(priorlenmean)) / ((yday(date) - yday(priordate))*log(priorlenmean)),
         sgr_len = 100*(exp(igr_len) - 1),
         
         agr_wt = (wtmean - priorwtmean) / (yday(date) - yday(priordate)),
         igr_wt = (log(wtmean) - log(priorwtmean)) / (yday(date) - yday(priordate)),
         rgr_wt = (wtmean - priorwtmean) / ((yday(date) - yday(priordate))*priorwtmean),
         rgrl_wt = (log(wtmean) - log(priorwtmean)) / ((yday(date) - yday(priordate))*log(priorwtmean)),
         sgr_wt = 100*(exp(igr_wt) - 1))

ggpairs(dat %>% dplyr::select(priorlenmean, agr_len, igr_len, rgr_len, rgrl_len, sgr_len))
ggpairs(dat %>% dplyr::select(priorwtmean, agr_wt, igr_wt, rgr_wt, rgrl_wt, sgr_wt))

plot(rgr_len ~ igr_len, dat, col = rbPal(100)[as.numeric(cut(dat$priorlenmean, breaks = 100))])
plot(rgrl_len ~ rgr_len, dat, col = rbPal(100)[as.numeric(cut(dat$priorlenmean, breaks = 100))])
plot(rgr_len ~ agr_len, dat, col = rbPal(10)[as.numeric(cut(dat$priorlenmean, breaks = 10))])

summary(lm(rgr_len ~ meantempper*effden_permeter_log, dat))
summary(lm(rgrl_len ~ meantempper*effden_permeter_log, dat))
summary(lm(igr_len ~ meantempper*effden_permeter_log, dat))

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
# dat <- read_csv("Growth_DataTable_WithCovariates.csv")
# head(dat)
# dat$stream <- as.factor(dat$stream)
# 
# tempsensdf <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Temperature/TempSensitivityWithYearlySiteMetrics.csv") %>%
#   group_by(site) %>% summarize(gw_strbuff = unique(gw_strbuff), 
#                                # gwindex = unique(gwindex),
#                                basin_areasqkm = unique(basin_areasqkm),
#                                slope_basin = unique(slope_basin),
#                                slope_stream = unique(slope_stream),
#                                forest_strbuff = unique(forest_strbuff),
#                                loglakearea = unique(loglakearea),
#                                ts = mean(slo.est)) %>% ungroup()
# 
# dat2 <- dat %>% left_join(tempsensdf, by = c("tempsite" = "site")) %>%
#   mutate(gwindex = gw_strbuff / log(basin_areasqkm),
#          loggwindex = log(gwindex),
#          relloggwindex = 1 - (loggwindex/min(loggwindex))
#          # loggwindex2_scaled = 1 - (abs(loggwindex2) / max(abs(loggwindex2)))
#          # loggwindex2_scaled = (abs(min(loggwindex2)) + loggwindex2) / max(abs(min(loggwindex2)) + loggwindex2)
#          )
# boxplot(unique(dat2$relloggwindex))
# boxplot(unique(dat2$gw_strbuff))
# view(dat2 %>% group_by(stream) %>% summarize(gwi = unique(loggwindex), amed = unique(Amedian)))
# ggpairs(dat2 %>% select(Amedian, ts, gw_strbuff, gwindex, relloggwindex))
# view(dat2)
# 
# 
# summary(lm(growthwtmean ~ doy + I(doy^2) + doy*Amedian, dat))
# summary(lm(prodwtmean ~ doy + I(doy^2) + doy*Amedian, dat))
# summary(lm(prodwtmean ~ doy + I(doy^2) + doy*Amedian, dat))
# 
# write_csv(dat2, "Growth_DataTable_WithCovariates.csv")

#### add consumption
# dat <- read_csv("Growth_DataTable_WithCovariates.csv")
# consump <- read_csv("Bioenergetics/YOYGrowth_BioE_DesignFile_Means.csv")
# consump2 <- read_csv("Bioenergetics/FB4_Log_File_BaldockYOYGrowthMeans.csv") %>% select(Run_Name, 'p-value') %>% rename(pvalue = 'p-value')
# consump <- consump %>% left_join(consump2) %>% select(stream, date, section, year, pvalue)
# dat <- dat %>% left_join(consump)
# write_csv(dat, "Growth_DataTable_WithCovariates.csv")
# 
# 
# dat <- read_csv("Growth_DataTable_WithCovariates.csv") %>% group_by(stream) %>% summarize(relloggwindex = unique(relloggwindex), ts = unique(ts))
# plot(ts ~ relloggwindex, dat, xlab = "Groundwater Influence", ylab = "Temperature Sensitivity")
# dat <- dat %>% arrange(relloggwindex)
# 
# par(mar = c(4,8,1,1))
# barplot(dat$relloggwindex, horiz = T, names.arg = dat$stream, las = 2, col = viridis::viridis(13), xlim = c(0,0.8))

# compare metrics of groundwater influence
ggpairs(dat %>% dplyr::select(basin_areasqkm, ts, Amedian, springprev_point_norm, springprev_basinmean_norm, springprev_iew01km_norm) %>%
          mutate(x = springprev_basinmean_norm / log(basin_areasqkm), 
                 lx = log(springprev_basinmean_norm / log(basin_areasqkm)), 
                 lxx = log(springprev_basinmean_norm)))
# springprev_basinmean_norm gives essentially the same information as previously used metrics of groundwater influence, such as that divided by log basin area, on a log scale
# springprev_basinmean_norm is negatively related to basin area (neg log relationship)
# springprev_basinmean_norm is negatively related to temp sensitivity and Amedian (neg log relationship)


#------------------------------------------------------#
# Bind data - Total growth
#------------------------------------------------------#

lenwt <- read_csv("YOY_LengthWeight_Summary.csv") #%>% mutate(year = year(date))
growthtot <- read_csv("YOY_SpecificGrowthRate_SummerTotal_Summary.csv") %>% 
  dplyr::select(stream, date, priordate, reach, year, growthwtmean, growthwtsd) %>% 
  rename(section = reach) %>% mutate(date = as_date(date))
dat3 <- dat %>% group_by(stream) %>% summarize(tempsite = unique(tempsite))

growthtot <- growthtot %>% left_join(lenwt, by = c("stream", "date", "section")) %>% left_join(dat3)
growthtot <- growthtot %>% left_join(tempsensdf, by = c("tempsite" = "site")) %>%
  mutate(gwindex = gw_strbuff / log(basin_areasqkm),
         loggwindex = log(gwindex),
         relloggwindex = 1 - (loggwindex/min(loggwindex))
  )

write_csv(growthtot, "Growth_SummerTotal_DataTable_WithCovariates.csv")
growthtot <- read_csv("Growth_SummerTotal_DataTable_WithCovariates.csv")

# plot and model effect of gwindex on total summer growth rate
summary(lm(growthwtmean ~ relloggwindex, growthtot))
plot(lm(growthwtmean ~ relloggwindex, growthtot))
ggplot(growthtot, aes(x = relloggwindex, y = growthwtmean)) + geom_point() + geom_smooth(method = "lm")

# plot and model effect of gwindex on of end of season body weight
summary(lm(wtmean ~ relloggwindex, growthtot))
plot(lm(wtmean ~ relloggwindex, growthtot))
ggplot(growthtot, aes(x = relloggwindex, y = wtmean)) + geom_point() + geom_smooth(method = "lm")

# plot and model effect of gwindex on of end of season body length
summary(lm(lenmean ~ relloggwindex, growthtot))
plot(lm(lenmean ~ relloggwindex, growthtot))
ggplot(growthtot, aes(x = relloggwindex, y = lenmean)) + geom_point() + geom_smooth(method = "lm")


#------------------------------------------------------#
# Exploratory plotting
#------------------------------------------------------#

rbPal <- colorRampPalette(c('red','blue'))
# dat$Col <- rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))]

# pairs plot of predictor variables
ggpairs(dat %>% dplyr::select(density, density_log, cpue, cpue_log, effden_permeter, effden_permeter_log, doy, priorlenmean, meantempper, springprev_basinmean_norm))
# pairs plot of response variables
ggpairs(dat %>% dplyr::select(growthwtmean, prodwtmean, growthlenmean, prodlenmean, pvalue))

dd <- dat %>% dplyr::select(specgrowthlen, specgrowthwt, doy, springprev_basinmean_norm, priorlenmean)
dd <- dd[complete.cases(dd),]
plot(specgrowthwt ~ specgrowthlen, dd, col = rbPal(100)[as.numeric(cut(dd$priorlenmean, breaks = 100))])
m1 <- (lm(specgrowthwt ~ specgrowthlen + priorlenmean, dat))
plot(m1$residuals ~ dd$doy, col = rbPal(100)[as.numeric(cut(dd$springprev_basinmean_norm, breaks = 100))])
abline(lm(m1$residuals ~ dd$doy))
plot(m1)

plot(wtmean ~ density, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(growthwtmean ~ density, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(prodwtmean ~ density, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])

plot(wtmean ~ logdens, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(growthwtmean ~ logdens, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(prodwtmean ~ logdens, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])

plot(wtmean ~ doy, dat, col = rbPal(100)[as.numeric(cut(dat$Amedian, breaks = 100))])
plot(growthwtmean ~ doy, dat, col = rbPal(100)[as.numeric(cut(dat$Amedian, breaks = 100))])
plot(prodwtmean ~ doy, dat, col = rbPal(100)[as.numeric(cut(dat$Amedian, breaks = 100))])
plot(density ~ doy, dat, col = rbPal(100)[as.numeric(cut(dat$Amedian, breaks = 100))])
plot(logdens ~ doy, dat, col = rbPal(100)[as.numeric(cut(dat$Amedian, breaks = 100))])

plot(wtmean ~ Amedian, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(growthwtmean ~ Amedian, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(growthwtmean ~ Amedian, dat, col = rbPal(100)[as.numeric(cut(dat$GDDper0, breaks = 100))])
plot(prodwtmean ~ Amedian, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(density ~ Amedian, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(logdens ~ Amedian, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])

plot(wtmean ~ Amedian, dat, col = rbPal(100)[as.numeric(cut(dat$logdens, breaks = 100))])
plot(growthwtmean ~ Amedian, dat, col = rbPal(100)[as.numeric(cut(dat$logdens, breaks = 100))])
plot(prodwtmean ~ Amedian, dat, col = rbPal(100)[as.numeric(cut(dat$logdens, breaks = 100))])

plot(wtmean ~ GDDcum0, dat, col = rbPal(100)[as.numeric(cut(dat$logdens, breaks = 100))])
plot(growthwtmean ~ GDDper0, dat, col = rbPal(100)[as.numeric(cut(dat$logdens, breaks = 100))])
plot(prodwtmean ~ GDDper0, dat, col = rbPal(100)[as.numeric(cut(dat$logdens, breaks = 100))])
plot(logdens ~ GDDper0, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])

plot(wtmean ~ avgredds, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(growthwtmean ~ avgredds, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(prodwtmean ~ avgredds, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(logdens ~ avgredds, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])

plot(wtmean ~ logredds, dat, col = rbPal(100)[as.numeric(cut(dat$Amedian, breaks = 100))])
plot(growthwtmean ~ logredds, dat, col = rbPal(100)[as.numeric(cut(dat$Amedian, breaks = 100))])
plot(prodwtmean ~ logredds, dat, col = rbPal(100)[as.numeric(cut(dat$Amedian, breaks = 100))])
plot(logdens ~ logredds, dat, col = rbPal(100)[as.numeric(cut(dat$Amedian, breaks = 100))])


plot(log(prodwtmean) ~ Amedian, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(log(prodwtmean) ~ GDDper0, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(log(prodwtmean) ~ doy, dat, col = rbPal(100)[as.numeric(cut(dat$Amedian, breaks = 100))])


plot(log(catch+1) ~ logredds, dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
plot(logdens ~ I(logredds/seclen), dat, col = rbPal(100)[as.numeric(cut(dat$doy, breaks = 100))])
summary(lm(I(log(catch+1)) ~ logredds, dat))
summary(lm(logdens ~ I(logredds/seclen), dat))


summary(lm(pvalue ~ effden_permeter_log, dat))

summary(lm(pvalue ~ doy + I(doy^2) + springprev_basinmean_norm + doy*springprev_basinmean_norm + I(doy^2)*springprev_basinmean_norm, dat))
summary(lm(specgrowthlen ~ doy + I(doy^2) + springprev_basinmean_norm + doy*springprev_basinmean_norm + I(doy^2)*springprev_basinmean_norm, dat))
summary(lm(specgrowthlen ~ pvalue * springprev_basinmean_norm, dat))


car::vif(lm(pvalue ~ effden_permeter_log + doy + springprev_basinmean_norm, dat))


#------------------------------------------------------#
# Seasonal variation in growth and production
#------------------------------------------------------#

cvdat <- dat %>% group_by(stream, year, section, Amedian, dum) %>% 
  summarize(gn = sum(!is.na(growthwtmean)),
            pn = sum(!is.na(prodwtmean)),
            gmean = mean(growthwtmean, na.rm = T), 
            gsd = sd(growthwtmean, na.rm = T), 
            pmean = mean(prodwtmean, na.rm = T), 
            psd = sd(prodwtmean, na.rm = T)) %>%
  mutate(gcv = gsd / gmean, pcv = psd / pmean) %>% ungroup() %>%
  filter(gn >= 3)

plot(gcv ~ Amedian, cvdat, bg = cols[cvdat$dum], pch = cvdat$year - 2000)
legend("topleft", legend = levels(unique(dat$stream)), fill = cols)
summary(lm(gcv ~ Amedian, cvdat))
abline(lm(gcv ~ Amedian, cvdat))
plot(lm(gcv ~ Amedian, cvdat))

plot(pcv ~ Amedian, cvdat, bg = cols[cvdat$dum], pch = cvdat$year - 2000)
legend("topleft", legend = levels(unique(dat$stream)), fill = cols)
summary(lm(pcv ~ Amedian, cvdat))
abline(lm(pcv ~ Amedian, cvdat))
plot(lm(gcv ~ Amedian, cvdat))

#------------------------------------------------------#
# Intra-class correlation plots
#------------------------------------------------------#

# streams <- levels(unique(dat$stream))
# years <- c(2021:2022)
# 
# growthicc <- tibble(year = c(rep(2021, 13), rep(2022, 13)), stream = rep(streams, 2), icc = NA)
# for (i in 1:length(years)) {
#   for (j in 1:length(streams)) {
#     d <- dat %>% filter(year == years[i] & stream == streams[j]) %>% dplyr::select(section, growthwtmean) %>% pivot_wider(., names_from = section, values_from = growthwtmean, values_fn = list) %>% unnest(cols = everything())
#     ic <- icc(d, model = "twoway", type = "agreement", unit = "single")
#     growthicc$icc[i*j] <- ic$value
#   }
# }
dat %>% ggplot() + geom_point(aes(x = section, y = growthwtmean, col = date)) + facet_wrap(~ stream + year, scales = "free")
dat %>% ggplot() + geom_point(aes(x = section, y = prodwtmean, col = date)) + facet_wrap(~ stream + year, scales = "free")


#------------------------------------------------------#
# Exploratory modeling
#------------------------------------------------------#

summary(m1 <- lm(growthwtmean ~ Amedian + GDDper0 + logdens + doy + logdens*doy + Amedian*doy, dat))
summary(m2 <- lm(growthwtmean ~ Amedian + GDDper0 + logdens + doy + logdens*doy + Amedian*doy + GDDper0*doy, dat))
summary(m2a <- lm(growthwtmean ~ Amedian + GDDper0 + logdens + doy + logdens*doy + Amedian*doy + GDDper0*doy + I(doy^2), dat))
summary(m2b <- lm(growthwtmean ~ Amedian + GDDper0 + logdens + doy + logdens*doy + Amedian*doy + GDDper0*doy + I(doy^2) + logdens*I(doy^2), dat))
summary(m2c <- lm(growthwtmean ~ Amedian + GDDper5 + logdens + doy + logdens*doy + Amedian*doy + GDDper5*doy, dat))

summary(m3 <- lm(growthwtmean ~ Amedian + GDDper0 + logdens + doy + logdens*doy + Amedian*doy + GDDper0*doy + Amedian*logdens, dat))
summary(m4 <- lm(growthwtmean ~ Amedian + logdens + doy + logdens*doy + Amedian*doy, dat))
summary(m5 <- lm(growthwtmean ~ Amedian + logdens + doy + meantempper, dat))

aictab(cand.set = list(m1, m2, m2a, m3, m4, m5), modnames = c("m1", "m2", "m2a", "m3", "m4", "m5"))
aictab(cand.set = list(m2, m2c), modnames = c("m2", "m2c"))


summary(m6 <- lm(prodwtmean ~ Amedian + GDDper0 + doy, dat))
summary(m7 <- lm(prodwtmean ~ Amedian + GDDper0 + doy + Amedian*doy, dat))
summary(m8 <- lm(prodwtmean ~ Amedian + GDDper0 + doy + GDDper0*doy, dat))
summary(m9 <- lm(prodwtmean ~ Amedian + GDDper0 + doy + Amedian*doy + GDDper0*doy, dat))
summary(m10 <- lm(prodwtmean ~ Amedian + GDDper0, dat))
summary(m11 <- lm(prodwtmean ~ Amedian, dat))

aictab(cand.set = list(m6, m7, m8, m9, m10, m11), modnames = c("m6", "m7", "m8", "m9", "m10", "m11"))


summary(m12 <- lm(growthwtmean ~ Amedian + doy + I(doy^2) + doy*Amedian, dat))
summary(m13 <- lm(growthwtmean ~ Amedian + doy + doy*Amedian, dat))
aictab(cand.set = list(m12, m13), modnames = c("m12", "m13"))

summary(lm(growthwtmean ~ Amedian + GDDper0 + doy + elev + elev*doy, dat))
