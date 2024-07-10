#==========================================================================================================================================#
# PROJECT NAME: Groundwater structures fish growth and production across a riverscape
# CREATOR: Jeffrey R. Baldock, WY Coop Unit, University of Wyoming, jbaldock@uwyo.edu 
# PURPOSE: Estimate effects of cathment characteristics (including groundwater) on stream temperature sensitivity to air temperature
# NOTES:
#==========================================================================================================================================#



library(GGally)
library(tidyverse)
library(lubridate)
library(plotfunctions)
library(RColorBrewer)
library(R2jags)
library(runjags)
library(MCMCvis)
library(loo)
library(HDInterval)
library(scales)
library(viridis)
library(ggmcmc)
library(mapview)
library(plotfunctions)
library(daymetr)


##---------------------------------------------------------------------------------------##
## Temperature data and Watershed metrics
##---------------------------------------------------------------------------------------##

# read in raw stream temp data
dat <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Environmental/Hydrology/StreamTempDatabase_05232023.csv") %>%
  mutate(doy = yday(date), year = year(date)) %>%
  group_by(stream, location, site, lat, long, date, year, doy) %>%
  summarise(meanT = mean(meanT), maxT = mean(maxT), minT = mean(minT)) %>% ungroup() %>% 
  filter(date >= "2021-01-01" & date <= "2022-12-31") %>% 
  filter(site %in% c("blackrock_1A", "blacktail_1A", "cabin_NA", "christian_NA", "cliff_1A", "coburn_NA", "cody_Moment",
                     "cottonwoodNPS_lower", "crystal_1A", "deadmans_1A", "dell_NA", "ditch_upper", "dog_NA", 
                     "fish_Upper_GV", "fish_downstream", "flat_4B", "goosewing_NA", "granite_NA", "grosventre_upper",
                     "hoback_upper", "lava_NA", "lowerbarbc_SG", "mosquito_1A", "nowlin_NA", "pacific_NA", "shoal_NA",
                     "slate_NA", "spread_1A", "spring_JLD", "spring_TSS", "srsidechannel_spawning", "threechannel_1A",
                     "threechannel_3A", "unkspring_Ford", "upperbarbc_SG", "willow_4B"))

# summarize site information
siteinfo <- dat %>% group_by(site) %>% 
  summarise(lat = unique(lat), long = unique(long), startd = min(date), endd = max(date), dayss = n()) %>% 
  ungroup() 
siteinfo <- siteinfo %>% mutate(siteid = 1:dim(siteinfo)[1])
dat <- dat %>% left_join(siteinfo %>% select(site, siteid))
dat$yearid <- ifelse(dat$year == 2021, 0, 1)

# download point location Daymet air temp data and join to stream temp data
sites <- siteinfo$site
airlist <- vector("list", length = length(sites))
for (i in 1:length(sites)) {
  si <- siteinfo %>% filter(site == sites[i])
  air <- download_daymet(site = si$site, lat = si$lat, lon = si$long, start = year(si$startd), end = year(si$endd), internal = T)
  air$data <- air$data %>% mutate(tmean = (tmax..deg.c. + tmin..deg.c.)/2, date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"))
  air2 <- as_tibble(air$data[,c(11,10)]) %>% rename(airt = tmean)
  air2$site <- air$site
  airlist[[i]] <- air2
}
airdf <- do.call(rbind, airlist)
dat <- dat %>% left_join(airdf)

#
dat <- dat %>% filter(year %in% c(2021,2022) & doy >= 213 & doy <= 319)
dat_sum <- dat %>% group_by(site, year) %>% summarize(ndays = n(), nperc = n()/107) %>% ungroup() %>% filter(nperc >= 0.7) %>% mutate(site_year = paste(site, year, sep = "_"))
dat <- dat %>% mutate(site_year = paste(site, year, sep = "_")) %>% filter(site_year %in% dat_sum$site_year)
# dat_study <- dat_study %>% left_join(dat_study_sum_fil[,c(5:6)])

# create lagged variables
dat <- dat %>% group_by(site, year) %>% mutate(meanTlag1 = lag(meanT, n = 1L), airtlag1 = lag(airt, n = 1L)) %>% ungroup()

### watershed metrics
watmet <- read_csv("Temperature/TempSensPCAdata.csv") %>% select(site, elev, slope, forest, logarea, loglake) %>% 
  group_by(site) %>% summarize(elev = unique(elev), slope = unique(slope), forest = unique(forest), logarea = unique(logarea), loglake = unique(loglake))
gwmet <- read_csv("Groundwater Metrics/GroundwaterMetrics_Normalized_TempPoints.csv") %>%
  select(sitename, springprev_basinmean_norm, springprev_iew01km_norm, springprev_iew05km_norm, springprev_basinmean_log_norm, springprev_iew01km_log_norm, springprev_iew05km_log_norm) %>%
  rename(site = sitename)
watmet <- watmet %>% left_join(gwmet %>% select(site, springprev_iew05km_norm)) %>% rename(gw = springprev_iew05km_norm)
watmet <- watmet %>% left_join(siteinfo %>% select(site, siteid))
watmet$loggw <- log(watmet$gw)
# z-score basin characteristics
watmet <- watmet %>% mutate(slopez = scale(slope, center = TRUE, scale = T),
                            forestz = scale(forest, center = TRUE, scale = T),
                            logareaz = scale(logarea, center = TRUE, scale = T),
                            loglakez = scale(loglake, center = TRUE, scale = T),
                            loggwz = scale(loggw, center = TRUE, scale = T))
# ggpairs(watmet %>% select(slopez, forestz, logareaz, loglakez, loggwz))


sitess <- dat %>% group_by(site) %>% summarize(siteid = unique(siteid))

write_csv(watmet %>% select(site, elev, slope, forest, logarea, loglake, gw) %>% mutate(area = exp(logarea), lake = exp(loglake)), "Temperature/TempSens_BasinMetrics.csv")

##---------------------------------------------------------------------------------------##
## Fit JAGS Model
##---------------------------------------------------------------------------------------##

Covs <- model.matrix(~ slopez + forestz + logareaz + loglakez + loggwz + slopez*loggwz + forestz*loggwz + logareaz*loggwz + loglakez*loggwz, data = watmet, na.rm = FALSE)
jags.data <- list("nObs" = dim(dat)[1], "yearid" = dat$yearid, "siteid" = as.factor(dat$siteid), "nSites" = length(unique(dat$siteid)), 
                  "airt" = dat$airt, "Covs" = Covs[,-1], "nCovs" = dim(Covs)[2]-1, "temp" = dat$meanT)
params <- c("alpha", "beta.year", "beta.air", "alpha.ts", "beta.ts", "sigma", "sigma.ts", "loglik", "mu.t", "temp","res", "var_fit", "var_res", "R2")
tsmod <- jags(jags.data, inits = NULL, parameters.to.save = params, model.file = "Temperature/tsmod2.txt", 
              n.chains = 3, n.thin = 20, n.burnin = 2000, n.iter = 22000, DIC = TRUE)
tsmod
tsmod$BUGSoutput$summary[,8][tsmod$BUGSoutput$summary[,8] > 1.01] # any R-hat values too high? 
MCMCtrace(tsmod, ind = TRUE, params = c("alpha", "beta.year", "beta.air", "alpha.ts", "beta.ts", "sigma", "sigma.ts", "R2"), 
          filename = "Temperature/TemperatureSensitivity_traceplots.pdf") # traceplots

saveRDS(tsmod, "Temperature/TempModelFit.rds")
tsmod <- readRDS("Temperature/TempModelFit.rds")

##---------------------------------------------------------------------------------------##
## MODEL DIAGNOSTICS/SELECTION - calculate and plot LOO
##---------------------------------------------------------------------------------------##

# LOO-CV to evaluate appropriate model specification
loglik <- tsmod$BUGSoutput$sims.list$loglik
reff <- relative_eff(exp(loglik), chain_id = c(rep(1,1000),rep(2,1000),rep(3,1000)))
loo_ts <- loo(loglik, r_eff = reff)
plot(loo_ts, label_points = TRUE)

# set the top model based on LOO
top_mod <- tsmod


##---------------------------------------------------------------------------------------##
## Save MCMC samples from top model (and second best model, top + cubic temp)
##---------------------------------------------------------------------------------------##

# Save MCMC samples from top model
top_mod.mcmc <- as.data.frame(as.matrix(as.mcmc(top_mod)))
modelout <- top_mod$BUGSoutput
# generate MCMC samples and store as an array
McmcList <- vector("list", length = dim(modelout$sims.array)[2])
for(i in 1:length(McmcList)) {
  McmcList[[i]] = as.mcmc(modelout$sims.array[,i,])
}
# rbind MCMC samples from 3 chains and save as object
Mcmcdat <- rbind(McmcList[[1]], McmcList[[2]], McmcList[[3]])
param.summary <- modelout$summary

# save model output
# write_csv(as.data.frame(Mcmcdat), "Temperature/TemperatureSensitivity_mcmcsamps.csv")
# Mcmcdat <- read.csv("Temperature/TemperatureSensitivity_mcmcsamps.csv")
# write.csv(as.data.frame(param.summary), "Temperature/TemperatureSensitivity_ParameterSummary.csv", row.names = T)
# param.summary <- read.csv("Temperature/TemperatureSensitivity_ParameterSummary.csv", row.names = 1)


##---------------------------------------------------------------------------------------##
## MODEL DIAGNOSTICS
##---------------------------------------------------------------------------------------##

Mcmcdat <- as_tibble(Mcmcdat)

# subset expected and observed MCMC samples
ppdat_exp <- as.matrix(Mcmcdat[,startsWith(names(Mcmcdat), "mu.t[")])
ppdat_obs <- as.matrix(Mcmcdat[,startsWith(names(Mcmcdat), "temp[")])

# Bayesian p-value
sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2])

# Posterior Predictive Check: plot median posterior expected length ~ observed length with Bayesian p-value
jpeg("Temperature/Figures/TemperatureSensitivity_PPCheck.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(x = seq(from = -5, to = 20, length.out = 100), y = seq(from = -5, to = 20, length.out = 100), pch = NA, xlab = expression(paste("Observed Temperature ("^"o", "C)", sep = "")), ylab = expression(paste("Median Posterior Expected Temperature ("^"o", "C)", sep = "")))
for (i in 1:length(dat$meanT)) { points(median(unlist(Mcmcdat[,paste("mu.t[",i,"]", sep = "")])) ~ dat$meanT[i])}
legend("topleft", bty = "n", legend = paste("Bayesian p-value = ", round(sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2]), digits = 3), sep = ""))
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()

# histogram of residuals
jpeg("Temperature/Figures/TemperatureSensitivity_ResidHist.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist((ppdat_obs - ppdat_exp), main = "", xlab = "Observed - Expected")
legend("topright", bty = "n", legend = expression(paste("Median = 0.045 "^"o", "C", sep = "")))
abline(v = median(unlist(ppdat_obs - ppdat_exp)), col = "red", lwd = 2)
dev.off()

# RMSE
sqrt(mean((ppdat_obs - ppdat_exp)^2)) # point estimate
rmse <- vector("numeric", 3000L)
for (j in 1:dim(Mcmcdat)[1]) {
  res <- vector("numeric", 258L)
  for (i in 1:length(dat$meanT)) {
    res[i] <- (ppdat_obs[j,i] - ppdat_exp[j,i])^2
  }
  rmse[j] <- sqrt(mean(res))
  print(j)
}
jpeg("Temperature/Figures/TemperatureSensitivity_RMSE.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist(rmse, xlab = expression(paste("Root Mean Square Error ("^"o", "C)", sep = "")), ylab = "Posterior Frequency", main = "")
legend("topright", bty = "n", legend = expression(paste("Median = 1.163 "^"o", "C", sep = "")))
abline(v = median(rmse), col = "red", lwd = 2)
dev.off()


# Bayesian R^2 per Gelman et al. (2018), Nakagawa and Shielzeth (2013)
round(median(unlist(Mcmcdat[,"R2"])), digits = 3)
round(hdi(Mcmcdat[,"R2"], credMass = 0.95)[1], 3)
round(hdi(Mcmcdat[,"R2"], credMass = 0.95)[2], 3)

jpeg("Temperature/Figures/TemperatureSensitivity_BayesianR2.jpg", units = "in", width = 5.5, height = 5, res = 1500)
par(mar = c(4,4,0.5,0.5), mgp = c(2.5,1,0))
hist(unlist(Mcmcdat[,"R2"]), main = NA, xlab = expression(Bayesian ~ R^2), ylab = "Posterior Frequency", xlim = c(0.898,0.908))
legend("topright", bty = "n", cex = 0.9, legend = "Median = 0.902 (0.900, 0.903)")
abline(v = median(unlist(Mcmcdat[,"R2"])), col = "red", lwd = 2)
dev.off()


##### Combined Figure
jpeg("Temperature/Figures/TemperatureSensitivity_ModelDiagnosticsCombined.jpg", units = "in", width = 7, height = 7, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), mfrow = c(2,2))

# PP Check
plot(x = seq(from = -5, to = 20, length.out = 100), y = seq(from = -5, to = 20, length.out = 100), bty = "l", pch = NA, xlab = expression(paste("Observed Temperature ("^"o", "C)", sep = "")), ylab = expression(paste("Median Posterior Expected Temperature ("^"o", "C)", sep = "")))
for (i in 1:length(dat$meanT)) { points(median(unlist(Mcmcdat[,paste("mu.t[",i,"]", sep = "")])) ~ dat$meanT[i])}
legend("top", bty = "n", legend = paste("Bayesian p-value = ", round(sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2]), digits = 3), sep = ""))
abline(a = 0, b = 1, col = "red", lwd = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(a)")
par(usr = usr)

# RMSE
hist(rmse, xlab = expression(paste("Root Mean Square Error ("^"o", "C)", sep = "")), ylab = "Posterior Frequency", main = "")
box(bty = "l")
legend("topright", bty = "n", legend = expression(paste("Median = 1.163 "^"o", "C", sep = "")))
abline(v = median(rmse), col = "red", lwd = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(b)")
par(usr = usr)

# Marg/Cond R2
hist(unlist(Mcmcdat[,"R2"]), main = NA, xlab = expression(Bayesian ~ R^2), ylab = "Posterior Frequency")
legend("topright", bty = "n", cex = 0.9, legend = c("Median = 0.902", "(0.900, 0.903)"))
abline(v = median(unlist(Mcmcdat[,"R2"])), col = "red", lwd = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(c)")
par(usr = usr)

dev.off()

##---------------------------------------------------------------------------------------##
## PLOTTING
##---------------------------------------------------------------------------------------##

# Parameter (beta) estimates - dot plot
mod.gg <- ggs(as.mcmc(top_mod))

# Temp Sensitivity - slopes
jpeg("Temperature/Figures/TemperatureSensitivity_ParameterDotPlot_BetasTempSens.jpg", units = "in", width = 5, height = 9, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^beta.air", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = sitess$site)
dev.off()

# Temp Sensitivity - slopes
jpeg("Temperature/Figures/TemperatureSensitivity_ParameterDotPlot_BetasWatershedEffects.jpg", units = "in", width = 5, height = 5, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^beta.ts", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = colnames(Covs)[-1])
dev.off()



# set up gradient legend
vrPal <- colorRampPalette(rev(viridis(12)))
gradpal <- tibble(gw = seq(from = 0, to = 1, length.out = 100)) %>%
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>% 
  filter(gw >= min(watmet$gw) & gw <= max(watmet$gw))
# number of values for prediction
nvalues <- 100
# generate sequences of predictor variables and scales axes
x.gw <- seq(from = min(watmet$loggwz), to = max(watmet$loggwz), length.out = nvalues)
# manually set up and transform axes to original scale of the data
range(watmet$gw)
x.axis.gw <- c(0.05, 0.1, 0.2, 0.4, 0.8) #seq(from = 0.05, to = 1.01, length.out = 5)
x.axis.gw.log <- log(x.axis.gw)
x.scaled.gw <- (x.axis.gw.log - mean(watmet$loggw)) / sd(watmet$loggw)

# tidy table
outdf <- tibble(site = watmet$site,
                loggwz = watmet$loggwz,
                tsmean = param.summary[c(34:64),1],
                tssd = param.summary[c(34:64),2],
                ts025 = param.summary[c(34:64),3],
                ts975 = param.summary[c(34:64),7])
outdf$coldum <- ifelse(outdf$site %in% c("blackrock_1A", "blacktail_1A", "cliff_1A", "crystal_1A", "fish_downstream", "flat_4B", "lowerbarbc_SG", 
                                         "mosquito_1A", "spread_1A", "threechannel_1A", "threechannel_3A", "upperbarbc_SG", "willow_4B"), "black", "white")

# predictions
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- unlist(Mcmcdat[i,"alpha.ts"]) + unlist(Mcmcdat[i,"beta.ts[5]"])*x.gw }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)

# plot
jpeg("Temperature/Figures/TemperatureSensitivity_GroundwaterEffect.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(tsmean ~ loggwz, outdf, pch = NA, bty = "l", ylim = c(0,1), xlab = "Groundwater index (log scale)", ylab = "Temperature sensitivity", axes = FALSE)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
polygon(c(x.gw, rev(x.gw)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty = 0)
lines(pred_median ~ x.gw, col = "black", lwd = 2)
arrows(x0 = outdf$loggwz, x1 = outdf$loggwz, y0 = outdf$ts025, y1 = outdf$ts975, length = 0)
points(tsmean ~ loggwz, outdf, pch = 21, bg = outdf$coldum)
legend("topright", pch = 21, pt.bg = c("black", "white"), legend = c("Focal streams", "Other streams"), bty = "n")
dev.off()



###### PLOT

# color palettes
vrPal <- colorRampPalette(rev(viridis(12)))
gradpal <- tibble(gw = seq(from = 0, to = 1, length.out = 100)) %>%
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>% 
  filter(gw >= min(watmet$gw) & gw <= max(watmet$gw))

# assign color based on gw influence: purple = snowmelt, yellow = spring-fed
watmet2 <- watmet %>% add_row(site = c("dum1","dum2"), gw = c(0,1)) %>% arrange(gw)
watmet2 <- watmet2 %>% group_by(site) %>% summarize(gw = mean(gw)) %>% arrange(gw) %>% mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))])
watmet2 <- watmet2 %>% filter(!site %in% c("dum1","dum2"))


jpeg(filename = "Temperature/Figures/TempTimeSeriesAndTempSensitivityByGroundwater_B.jpg", height = 2.25, width = 7, units = "in", res = 1500)
par(mfrow = c(1,3), mar = c(3.5,3.5,0.5,1), mgp = c(2.2,0.8,0))

sites <- unique(dat$site)
dat2 <- dat %>% left_join(watmet2)
dat2$cols <- ifelse(dat2$site %in% c("blackrock_1A", "blacktail_1A", "cliff_1A", "crystal_1A", "fish_downstream", "flat_4B", "lowerbarbc_SG", "mosquito_1A", "spread_1A", "threechannel_1A", "threechannel_3A", "upperbarbc_SG", "willow_4B"), dat2$cols, "grey70")
dat2$lwds <- ifelse(dat2$site %in% c("blackrock_1A", "blacktail_1A", "cliff_1A", "crystal_1A", "fish_downstream", "flat_4B", "lowerbarbc_SG", "mosquito_1A", "spread_1A", "threechannel_1A", "threechannel_3A", "upperbarbc_SG", "willow_4B"), 2, 1)
dat2a <- dat2 %>% filter(site %in% c("blackrock_1A", "blacktail_1A", "cliff_1A", "crystal_1A", "fish_downstream", "flat_4B", "lowerbarbc_SG", "mosquito_1A", "spread_1A", "threechannel_1A", "threechannel_3A", "upperbarbc_SG", "willow_4B"))
dat2b <- dat2 %>% filter(!site %in% c("blackrock_1A", "blacktail_1A", "cliff_1A", "crystal_1A", "fish_downstream", "flat_4B", "lowerbarbc_SG", "mosquito_1A", "spread_1A", "threechannel_1A", "threechannel_3A", "upperbarbc_SG", "willow_4B"))

# 2021
plot(meanT ~ date, dat2, col = NA, xlim = c(date("2021-08-01"), date("2021-11-05")), xlab = "Time", ylab = expression(paste("Mean daily temperature ("^"o","C)")), bty = "l")
sites <- unique(dat2b$site)
for (i in 1:length(sites)) {
  d <- dat2b %>% filter(site == sites[i] & year == 2021 & date <= date("2021-11-05"))
  lines(meanT ~ date, d, col = d$cols, lwd = d$lwds)
}
sites <- unique(dat2a$site)
for (i in 1:length(sites)) {
  d <- dat2a %>% filter(site == sites[i] & year == 2021 & date <= date("2021-11-05"))
  lines(meanT ~ date, d, col = d$cols, lwd = d$lwds)
}
legend("topright", legend = c("(a)"), pch = NA, bty = "n")
gradientLegend(valRange = c(min(gradpal$gw), max(gradpal$gw)), color = gradpal$cols, length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.05, 0.05, 0.12, 0.4), dec = 2)

# 2022
plot(meanT ~ date, dat2, col = NA, xlim = c(date("2022-08-01"), date("2022-11-05")), xlab = "Time", ylab = expression(paste("Mean daily temperature ("^"o","C)")), bty = "l")
sites <- unique(dat2b$site)
for (i in 1:length(sites)) {
  d <- dat2b %>% filter(site == sites[i] & year == 2022 & date <= date("2022-11-05"))
  lines(meanT ~ date, d, col = d$cols, lwd = d$lwds)
}
sites <- unique(dat2a$site)
for (i in 1:length(sites)) {
  d <- dat2a %>% filter(site == sites[i] & year == 2022 & date <= date("2022-11-05"))
  lines(meanT ~ date, d, col = d$cols, lwd = d$lwds)
}
legend("topright", legend = c("(b)"), pch = NA, bty = "n")


### TS ~ Groundwater Index

# pred values
nvalues <- 100 # number of values for prediction
x.gw <- seq(from = min(watmet$loggwz), to = max(watmet$loggwz), length.out = nvalues) # generate sequences of predictor variables and scales axes
range(watmet$gw) # manually set up and transform axes to original scale of the data
x.axis.gw <- c(0.05, 0.1, 0.2, 0.4, 0.8) #seq(from = 0.05, to = 1.01, length.out = 5)
x.axis.gw.log <- log(x.axis.gw)
x.scaled.gw <- (x.axis.gw.log - mean(watmet$loggw)) / sd(watmet$loggw)

# tidy table
outdf <- tibble(site = watmet$site,
                loggwz = watmet$loggwz,
                tsmean = param.summary[c(34:64),1],
                tssd = param.summary[c(34:64),2],
                ts025 = param.summary[c(34:64),3],
                ts975 = param.summary[c(34:64),7])
outdf <- outdf %>% left_join(watmet2)
outdf$cols <- ifelse(outdf$site %in% c("blackrock_1A", "blacktail_1A", "cliff_1A", "crystal_1A", "fish_downstream", "flat_4B", "lowerbarbc_SG", 
                                       "mosquito_1A", "spread_1A", "threechannel_1A", "threechannel_3A", "upperbarbc_SG", "willow_4B"), outdf$cols, "white")

# predictions
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- unlist(Mcmcdat[i,"alpha.ts"]) + unlist(Mcmcdat[i,"beta.ts[5]"])*x.gw }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)

# plot
# jpeg("Temperature/Figures/TemperatureSensitivity_GroundwaterEffect.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
# par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(tsmean ~ loggwz, outdf, pch = NA, bty = "l", ylim = c(0,1), xlab = "Groundwater index (log scale)", ylab = "Temperature sensitivity", axes = FALSE)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
polygon(c(x.gw, rev(x.gw)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty = 0)
lines(pred_median ~ x.gw, col = "black", lwd = 2)
arrows(x0 = outdf$loggwz, x1 = outdf$loggwz, y0 = outdf$ts025, y1 = outdf$ts975, length = 0)
points(tsmean ~ loggwz, outdf, pch = 21, bg = outdf$cols)
legend("topright", legend = c("(c)"), pch = NA, bty = "n")
# legend("topright", pch = c(NA,21,21), pt.bg = c(NA,"black", "white"), legend = c("(c)", "Focal streams", "Other streams"), bty = "n")

dev.off()

