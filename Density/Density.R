#==========================================================================================================================================#
# PROJECT NAME: Groundwater structures fish growth and production across a riverscape
# CREATOR: Jeffrey R. Baldock, WY Coop Unit, University of Wyoming, jbaldock@uwyo.edu 
# PURPOSE: Manipulate raw growth data to get density estimates
# NOTES:
# - currently density is YOY per meter of stream...
#   - may want to consider est. density per square meter
#   - potentially measure stream width using NAIP imagery??
# - YOY almost exclusively use bank habitat...
#   - per meter may be more relevant than per square meter?
# - First bit of code is data manipulation of the raw catch data. Jump to line 230
# - Lines 285-522 attempt to model density by groundwater using various polynomial functions. These do not fit the data well
# - Lines 530+ describe association between groundwater and density using LOESS, as in the manuscript
#==========================================================================================================================================#

library(tidyverse)
library(lubridate)
library(RColorBrewer)
library(ggpubr)
library(GGally)
library(R2jags)
library(runjags)
library(MCMCvis)
library(loo)
library(HDInterval)
library(scales)

setwd("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth")
getwd()

dat <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/GrowthData_working_withAges_YOYonly.csv") %>%
  mutate(section = factor(section)) %>%
  mutate(stream = factor(stream)) %>%
  mutate(streamnum = as.numeric(stream)) %>%
  mutate(year = year(date)) %>%
  mutate(doy = yday(date)) %>%
  mutate(elap = as.numeric(elap))

# calculate average elapsed sampling time to impute missing values
meanelap <- dat %>% group_by(stream, section) %>% summarize(meanelap = mean(elap, na.rm = T), sdelap = sd(elap, na.rm = T)) %>% ungroup() %>% mutate(cvelap = sdelap / meanelap)
boxplot(meanelap$cvelap)
plot(sdelap ~ meanelap, meanelap)

# section lengths: quick and dirty...could snap coordinates to stream network and calc stream dist...but may not work for spring streams
xx <- tibble(stream = rep(unique(dat$stream), each = 4), section = rep(c(1:4), times = length(unique(dat$stream))), seclen = 150)
xx$section <- as_factor(xx$section)
# fix short sections
xx$seclen[2] <- 111 # UBBC 4
xx$seclen[6] <- 111 # Blacktail 2

# complete set of streams, dates, and sections
nsec <- dat %>% group_by(stream, date) %>% summarize(nsec = n_distinct(section)) %>% ungroup()
xxx <- tibble(stream = rep(nsec$stream, each = 4), date = rep(nsec$date, each = 4), section = rep(c(1:4), times = 125))
xxx$section <- as_factor(xxx$section)
xxx <- xxx %>% left_join(xx, by = c("stream", "section"))

# Adjust incomplete sampling on some sections
xxx$seclen[xxx$stream == "Upper Bar BC" & xxx$date == "2022-08-25" & xxx$section == 1] <- 106  # UBBC, 8/25/22, section 1...stopped short due to backpack issue
xxx$seclen[xxx$stream == "Spread" & xxx$date == "2021-10-04" & xxx$section == 4] <- 80         # Spread, 10/4/21, section 4...battery died mid section
xxx$seclen[xxx$stream == "3 Channel" & xxx$date == "2022-10-05" & xxx$section == 4] <- 103     # 3 Channel, 10/5/22, section 5...stopped at 43.55562, -110.78717

# drop sections not sampled on certain dates
xxx <- xxx[!(xxx$stream == "Fall" & xxx$date == "2022-08-02" & xxx$section == 4),]
xxx <- xxx[!(xxx$stream == "Lower Bar BC" & xxx$date == "2021-08-21" & xxx$section == 4),]

# # add rows for stream/dates sampled, but no fish caught inn any section
# xxx <- xxx %>% 
#   add_row(stream = "Spread", date = date("2022-10-26"), section = as_factor(1), seclen = 150) %>%
#   add_row(stream = "Spread", date = date("2022-10-26"), section = as_factor(2), seclen = 150) %>%
#   add_row(stream = "Spread", date = date("2022-10-26"), section = as_factor(3), seclen = 150) %>%
#   add_row(stream = "Blackrock", date = date("2021-09-16"), section = as_factor(1), seclen = 150) %>%
#   add_row(stream = "Blackrock", date = date("2021-09-16"), section = as_factor(2), seclen = 150) %>%
#   add_row(stream = "Blackrock", date = date("2021-09-16"), section = as_factor(3), seclen = 150) %>%
#   add_row(stream = "Blackrock", date = date("2021-09-16"), section = as_factor(4), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-08-12"), section = as_factor(1), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-08-12"), section = as_factor(2), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-08-12"), section = as_factor(3), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-08-12"), section = as_factor(4), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-09-01"), section = as_factor(1), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-09-01"), section = as_factor(2), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-09-01"), section = as_factor(3), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-09-01"), section = as_factor(4), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-09-17"), section = as_factor(1), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-09-17"), section = as_factor(2), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-09-17"), section = as_factor(3), seclen = 150) %>%
#   add_row(stream = "Willow", date = date("2021-09-17"), section = as_factor(4), seclen = 150) 

# summarize catch data
abun <- dat %>% group_by(stream, date, year, section, elap) %>% summarize(catch = sum(species == "SRC")) %>% filter(section %in% c(1:4)) %>% ungroup() %>% left_join(dplyr::select(meanelap, stream, section, meanelap))
abun$elap <- ifelse(is.na(abun$elap), abun$meanelap, abun$elap)

databun <- xxx %>% left_join(dplyr::select(abun, stream, date, section, elap, catch), by = c("stream", "date", "section"))
databun$catch[is.na(databun$catch)] <- 0
databun$density <- databun$catch / databun$seclen  # YOY density, per meter
databun$cpue <- databun$catch / databun$elap # CPUE: YOY catch per minute
databun$cpue[is.na(databun$cpue)] <- 0

#-------------#
# plot
#-------------#

databun <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/YOY_CatchDensity_2021-2022.csv")
databun$section <- as.factor(databun$section)
elapp <- dat %>% group_by(stream, date, section) %>% summarize(elap = mean(elap)) %>% ungroup()
databun <- databun %>% left_join(elapp) %>% left_join(meanelap %>% select(stream, section, meanelap))
databun$elap <- ifelse(is.na(databun$elap), databun$meanelap, databun$elap)
databun <- databun %>% select(stream, date, section, seclen, catch, density, elap) %>% mutate(cpue = catch / elap)


# catch
databun %>% filter(date < "2022-01-01") %>%
  ggplot(aes(x = date, y = catch, color = section)) +
  geom_line() + geom_point() +
  facet_wrap(~ stream, scales = "free")

databun %>% filter(date > "2022-01-01") %>%
  ggplot(aes(x = date, y = catch, color = section)) +
  geom_line() + geom_point() +
  facet_wrap(~ stream, scales = "free")

# density
databun %>% filter(date < "2022-01-01") %>%
  ggplot(aes(x = date, y = density, color = section)) +
  geom_line() + geom_point() +
  facet_wrap(~ stream, scales = "free")

databun %>% filter(date > "2022-01-01") %>%
  ggplot(aes(x = date, y = density, color = section)) +
  geom_line() + geom_point() +
  facet_wrap(~ stream, scales = "free")

# logged density - modeling growth (or whatever) as a function of logged density may be most appropriate given 
# well established negative power relationship between density and body size/growth (see Imre and Grant papers and related)
# density-dependence strongest at relatively low densities
databun %>% filter(date < "2022-01-01") %>%
  ggplot(aes(x = date, y = log(density+0.00001), color = section)) +
  geom_line() + geom_point() +
  facet_wrap(~ stream)

databun %>% filter(date > "2022-01-01") %>%
  ggplot(aes(x = date, y = log(density+0.00001), color = section)) +
  geom_line() + geom_point() +
  facet_wrap(~ stream)

# relationship between density and cpue per RA comment
databun2 <- databun %>% mutate(logcpue = log(cpue), logdensity = log(density)) 
databun2$logcpue[is.infinite(databun2$logcpue)] <- NA
databun2$logdensity[is.infinite(databun2$logdensity)] <- NA

jpeg("DensityCpueRegression.jpg", width = 4.75, height = 4, units = "in", res = 1000)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(logcpue ~ logdensity, databun2, xlab = "log(Density, YOY per meter)", ylab = "log(CPUE, YOY per minute)")
dpmod <- (lm(logcpue ~ logdensity, databun2, na.action = na.exclude))
legend("bottomright", bty = "n", legend = paste("y = ", round(coefficients(dpmod)[1], digits = 2), " + ", round(coefficients(dpmod)[2], digits = 2), "*x", sep = ""))
abline(dpmod, col = "red", lwd = 3)
dev.off()

plot(cpue ~ density, databun2, xlab = "Density, YOY per meter", ylab = "CPUE, YOY per minute")
dpmod <- (lm(logcpue ~ logdensity, databun2, na.action = na.exclude))
legend("bottomright", bty = "n", legend = paste("y = ", round(coefficients(dpmod)[1], digits = 2), " + ", round(coefficients(dpmod)[2], digits = 2), "*x", sep = ""))


#-------------#
# save file
#-------------#

write_csv(databun, "YOY_CatchDensity_2021-2022.csv")


#-------------#
# Effective Density per Walters and Post (1993), Huntsman et al (2021)
#-------------#


streams <- unique(dat$stream)
years <- unique(dat$year)
lst1 <- list()
lst2 <- list()
lst3 <- list()
lst4 <- list()

# dat2 <- dat %>% filter(stream == "3 Channel", date == "2022-10-05")

for(m in 1:length(years)) {
  for(l in 1:length(streams)) {
    dat2 <- dat %>% filter(stream == streams[l] & year == years[m])
    dates <- unique(dat2$date)
    if(dim(dat2)[1] == 0) next
    for(i in 1:length(dates)) {
      for(k in 1:4) {
        dat3 <- dat2 %>% filter(date == dates[i] & section == k)
        if(dim(dat3)[1] == 0) next
        lens <- dat3$lengthmm[complete.cases(dat3$lengthmm)]
        dat3$effden <- dat3$lengthmm^2
        dat3$effden <- ifelse(is.na(dat3$effden), sample(lens, 1)^2, dat3$effden) # bootstrap unobserved size
        lst1[[k]] <- dat3
      }
      lst2[[i]] <- do.call(rbind, lst1)
    }
    lst3[[l]] <- do.call(rbind, lst2)
  }
  lst4[[m]] <- do.call(rbind, lst3)
}
dat4 <- do.call(rbind, lst4)

# view effective density output
range(dat4$lengthmm, na.rm = T)
sqrt(range(dat4$effden, na.rm = T))
plot(effden ~ lengthmm, dat4)
view(dat4)
summary(dat4)

# summarize by stream/section and calculate effective abundance per sampling event by summing effective abundance per fish
effabun <- dat4 %>% group_by(stream, date, year, section) %>% summarize(effabun = sum(effden)) %>% ungroup() 

# join to density/cpue data
databun <- databun %>% left_join(effabun) 
databun$effabun <- ifelse(is.na(databun$effabun), 0, databun$effabun)
databun$effden_permeter = databun$effabun / databun$seclen

# visualize
boxplot(databun$effden_permeter)
boxplot(log(databun$effden_permeter+0.0001))
ggpairs(dplyr::select(databun, density, cpue, effden_permeter) %>% mutate(logdensity = log(density), logcpue = log(cpue), logeffden = log(effden_permeter)))


write_csv(databun, "YOY_CatchDensity_2021-2022.csv")



##---------------------------------------------------------------------------------------##
# Density dynamics by groundwater influence
##---------------------------------------------------------------------------------------##

dat <- read_csv("YOY_CatchDensity_2021-2022.csv") 


# groundwater influence
gwmet <- read_csv("Groundwater Metrics/GroundwaterMetrics_Normalized_YOYsections.csv") %>% 
  separate(site, into = c("stream", "section"), sep = "_", remove = FALSE) %>% 
  separate(section, into = c("section", "updn"), sep = 1, remove = TRUE) %>% 
  mutate(section = as.numeric(section)) %>%
  dplyr::select(stream, section, springprev_iew05km_norm) %>% 
  rename(gw = springprev_iew05km_norm) %>% 
  add_row(stream = c("dum1","dum2"), section = c(NA,NA), gw = c(0,1)) %>% arrange(gw)

# assign color based on gw influence: purple = snowmelt, yellow = spring-fed
vrPal <- colorRampPalette(rev(viridis(12)))
gwmet <- gwmet %>% #group_by(stream) %>% summarize(gw = mean (gw)) %>% 
  arrange(gw) %>% mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))])
gwmet <- gwmet %>% filter(!stream %in% c("dum1","dum2"))

# view barplot...does ordering of gw influence make sense?
par(mar = c(5,10,2,2))
barplot(gwmet$gw, col = gwmet$cols, names.arg = gwmet$stream, horiz = TRUE, las = 2)
dev.off()

# join to density data
dat <- dat %>% left_join(gwmet)
dat$yrid <- ifelse(year(dat$date) == 2021, 21, 24)
dat$density_log <- log(dat$density + 0.01)

# maximum density by gw...probably the better approach
dat2 <- dat %>% filter(year %in% c(2021:2022)) %>% 
  group_by(stream, section, gw, cols, year, yrid) %>% 
  summarize(maxden = max(density), sumden = sum(density)) %>% 
  mutate(density_log = log(maxden + 0.01),
         gw_log = log(gw))

# plot(maxden ~ gw, dat2, ylim = c(0,2))
# plot(log(maxden) ~ gw, dat2)
# plot(maxden ~ log(gw), dat2)
# plot(log(maxden) ~ log(gw), dat2)
# 
# plot(density_log ~ gw, dat2, pch = NA, xlab = "Groundwater Influence", ylab = "log(Max. Density, fish per meter)")
# points(density_log ~ gw, dat2, pch = dat2$yrid, bg = dat2$cols)
# denmod <- lm(density_log ~ gw + I(gw^2), dat2)
# summary(denmod)
# gwvec <- seq(from = min(dat2$gw), to = max(dat2$gw), length.out = 100)
# preds <- predict(denmod, newdata = list(gw = gwvec), interval = "confidence")
# polygon(x = c(gwvec, rev(gwvec)), y = c(preds[,2], rev(preds[,3])), col = alpha("black", 0.2), lty = 0)
# lines(preds[,1] ~ gwvec, lwd = 2)
# legend("topleft", legend = c(2021,2022), pch = c(21,24), bty = "n")


##---------------------------------------------------------------------------------------##
## Misc JAGS fields
##---------------------------------------------------------------------------------------##

# grouping variables as factors for JAGS model
dat2$stream <- as.factor(dat2$stream)
dat2$strid <- as.factor(as.numeric(dat2$stream))
dat2$yrid <- as.factor(as.numeric(as.factor(dat2$year)))
dat2$sectid <- as.factor(paste(dat2$strid, dat2$section, sep = "."))

# tibble linking each section to each stream for JAGS model with section nested within stream
unq_sect_str <- dat2 %>% group_by(strid) %>% reframe(sectid = unique(sectid)) %>% ungroup()

# some additional data for JAGS
nYears <- length(unique(dat2$yrid))
nSections <- length(unique(dat2$sectid))
nStreams <- length(unique(dat2$strid))

StreamNames <- dat2 %>% group_by(stream) %>% summarize(strid = unique(strid), gw = mean(gw)) %>% arrange(gw) %>% ungroup() #%>% mutate(cols = (viridis::viridis(13)), strid2 = 1:13)
StreamNames <- StreamNames %>% arrange(strid)

# Sort 
dat2 <- dat2 %>% arrange(yrid, strid, sectid)
dat2 <- dat2 %>% ungroup() %>% mutate(gwz = scale(gw, center = TRUE, scale = TRUE))


##---------------------------------------------------------------------------------------##
## Fit JAGS Model
##---------------------------------------------------------------------------------------##

summary(lm(density_log ~ gwz + I(gwz^2), data = dat2))

# linear
Covs <- model.matrix(~ gwz, data = dat2, na.rm = FALSE)
jags.data <- list("nObs" = dim(dat2)[1], "dens" = dat2$density_log, "Covs" = matrix(Covs[,-1]), "nCovs" = dim(Covs)[2]-1,                                       
                  "strid" = dat2$strid, "yrid" = dat2$yrid, "sectid" = dat2$sectid, "sect_per_str" = unq_sect_str$strid,                                                                        
                  "nYears" = nYears, "nStreams" = nStreams, "nSections" = nSections)
params <- c("alpha.adj", "alpha.year.adj", "alpha.stream.adj", "alpha.sect.adj", "beta", "sigma", "sigma.yr", "sigma.str", "sigma.sect", "sigma.a", "sigma.b", "loglik", "dens", "mu")
denmod1 <- jags(jags.data, inits = NULL, parameters.to.save = params, model.file = "Density/densitymod.txt", 
               n.chains = 3, n.thin = 50, n.burnin = 50000, n.iter = 200000, DIC = TRUE)
denmod1
denmod1$BUGSoutput$summary[,8][denmod1$BUGSoutput$summary[,8] > 1.01] # any R-hat values too high? 

# quadratic
Covs <- model.matrix(~ gwz + I(gwz^2), data = dat2, na.rm = FALSE)
jags.data <- list("nObs" = dim(dat2)[1], "dens" = dat2$density_log, "Covs" = Covs[,-1], "nCovs" = dim(Covs)[2]-1,                                       
                  "strid" = dat2$strid, "yrid" = dat2$yrid, "sectid" = dat2$sectid, "sect_per_str" = unq_sect_str$strid,                                                                        
                  "nYears" = nYears, "nStreams" = nStreams, "nSections" = nSections)
params <- c("alpha.adj", "alpha.year.adj", "alpha.stream.adj", "alpha.sect.adj", "beta", "sigma", "sigma.yr", "sigma.str", "sigma.sect", "sigma.a", "sigma.b", "loglik", "dens", "mu")
denmod2 <- jags(jags.data, inits = NULL, parameters.to.save = params, model.file = "Density/densitymod.txt", 
                n.chains = 3, n.thin = 50, n.burnin = 50000, n.iter = 200000, DIC = TRUE)
denmod2
denmod2$BUGSoutput$summary[,8][denmod2$BUGSoutput$summary[,8] > 1.01] # any R-hat values too high? 

# cubic
Covs <- model.matrix(~ gwz + I(gwz^2) + I(gwz^3), data = dat2, na.rm = FALSE)
jags.data <- list("nObs" = dim(dat2)[1], "dens" = dat2$density_log, "Covs" = Covs[,-1], "nCovs" = dim(Covs)[2]-1,                                       
                  "strid" = dat2$strid, "yrid" = dat2$yrid, "sectid" = dat2$sectid, "sect_per_str" = unq_sect_str$strid,                                                                        
                  "nYears" = nYears, "nStreams" = nStreams, "nSections" = nSections)
params <- c("alpha.adj", "alpha.year.adj", "alpha.stream.adj", "alpha.sect.adj", "beta", "sigma", "sigma.yr", "sigma.str", "sigma.sect", "sigma.a", "sigma.b", "loglik", "dens", "mu")
denmod3 <- jags(jags.data, inits = NULL, parameters.to.save = params, model.file = "Density/densitymod.txt", 
                n.chains = 3, n.thin = 50, n.burnin = 50000, n.iter = 200000, DIC = TRUE)
denmod3
denmod3$BUGSoutput$summary[,8][denmod3$BUGSoutput$summary[,8] > 1.01] # any R-hat values too high? 

# quartic
Covs <- model.matrix(~ gwz + I(gwz^2) + I(gwz^3) + I(gwz^4), data = dat2, na.rm = FALSE)
jags.data <- list("nObs" = dim(dat2)[1], "dens" = dat2$density_log, "Covs" = Covs[,-1], "nCovs" = dim(Covs)[2]-1,                                       
                  "strid" = dat2$strid, "yrid" = dat2$yrid, "sectid" = dat2$sectid, "sect_per_str" = unq_sect_str$strid,                                                                        
                  "nYears" = nYears, "nStreams" = nStreams, "nSections" = nSections)
params <- c("alpha.adj", "alpha.year.adj", "alpha.stream.adj", "alpha.sect.adj", "beta", "sigma", "sigma.yr", "sigma.str", "sigma.sect", "sigma.a", "sigma.b", "loglik", "dens", "mu")
denmod4 <- jags(jags.data, inits = NULL, parameters.to.save = params, model.file = "Density/densitymod.txt", 
                n.chains = 3, n.thin = 50, n.burnin = 50000, n.iter = 200000, DIC = TRUE)
denmod4
denmod4$BUGSoutput$summary[,8][denmod4$BUGSoutput$summary[,8] > 1.01] # any R-hat values too high? 

# combine into list
fitlist <- list(denmod1, denmod2, denmod3, denmod4)


##---------------------------------------------------------------------------------------##
## MODEL DIAGNOSTICS/SELECTION - calculate and plot LOO
##---------------------------------------------------------------------------------------##

loo_list <- list()
for (i in 1:length(fitlist)) {
  m <- fitlist[[i]]
  loglik <- m$BUGSoutput$sims.list$loglik
  reff <- relative_eff(exp(loglik), chain_id = c(rep(1,3000),rep(2,3000),rep(3,3000)))
  loo_list[[i]] <- loo(loglik, r_eff = reff)
  print(i)
}
lc <- loo_compare(loo_list)
print(lc, simplify = FALSE, digits = 2)
plot(loo_list[[3]])
write.csv(data.frame(lc), "/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/Density/DensityByGroundwater_LOO_output.csv")

# set the top model based on LOO
top_mod <- Fit_list_gr[[1]]
Covs <- Covs_list[[1]]

# set the top model based on LOO
top_mod <- denmod3

MCMCtrace(denmod3, ind = TRUE, params = c("alpha.adj", "alpha.year.adj", "alpha.stream.adj", "alpha.sect.adj", "beta", "sigma", "sigma.yr", "sigma.str", "sigma.sect", "sigma.a", "sigma.b"), 
          filename = "Density/DensityGW_traceplots.pdf") # traceplots

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
write_csv(as.data.frame(Mcmcdat), "Density/DensityGroundwater_mcmcsamps.csv")
Mcmcdat <- read.csv("Density/DensityGroundwater_mcmcsamps.csv")
write.csv(as.data.frame(param.summary), "Density/DensityGroundwater_ParameterSummary.csv", row.names = T)
param.summary <- read.csv("Density/DensityGroundwater_ParameterSummary.csv", row.names = 1)


##---------------------------------------------------------------------------------------##
## Model Diagnostics
##---------------------------------------------------------------------------------------##

# subset expected and observed MCMC samples
ppdat_exp <- as.matrix(Mcmcdat[,startsWith(names(Mcmcdat), "mu.")])
ppdat_obs <- as.matrix(Mcmcdat[,startsWith(names(Mcmcdat), "dens.")])

# Bayesian p-value
sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2])

# Posterior Predictive Check: plot median posterior expected length ~ observed length with Bayesian p-value
jpeg("Density/Figures/DensityByGroundwater_PPCheck.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(x = seq(from = -4, to = 2, length.out = 100), y = seq(from = -4, to = 2, length.out = 100), pch = NA, xlab = "Observed Length (mm)", ylab = "Median Posterior Expected Length (mm)")
for (i in 1:dim(Covs)[1]) { points(median(Mcmcdat[,paste("mu.",i,".", sep = "")]) ~ median(Mcmcdat[,paste("dens.",i,".", sep = "")]))}
legend("topleft", bty = "n", legend = paste("Bayesian p-value = ", round(sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2]), digits = 3), sep = ""))
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()

# histogram of residuals
jpeg("Density/Figures/DensityByGroundwater_ResidHist.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist((ppdat_obs - ppdat_exp), main = "", xlab = "Observed - Expected")
legend("topright", bty = "n", legend = paste("Median Resid. = ", round(median((ppdat_obs - ppdat_exp)), digits = 2), sep = ""))
abline(v = median(unlist(ppdat_obs - ppdat_exp)), col = "red", lwd = 2)
dev.off()

# RMSE
sqrt(mean((ppdat_obs - ppdat_exp)^2)) # point estimate
rmse <- vector("numeric", 9000L)
for (j in 1:dim(Mcmcdat)[1]) {
  res <- vector("numeric", 97L)
  for (i in 1:dim(Covs)[1]) {
    res[i] <- (ppdat_obs[j,i] - ppdat_exp[j,i])^2
  }
  rmse[j] <- sqrt(mean(res))
  print(j)
}
jpeg("Density/Figures/DensityByGroundwater_RMSE.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist(rmse, xlab = "Root Mean Square Error", ylab = "Posterior Frequency", main = "")
legend("topright", bty = "n", legend = paste("Median RMSE = ", round(median(rmse), digits = 2), sep = ""))
abline(v = median(rmse), col = "red", lwd = 2)
dev.off()


##---------------------------------------------------------------------------------------##
## PARAMETER DOT PLOTS
##---------------------------------------------------------------------------------------##

# Parameter (beta) estimates - dot plot
mod.gg <- ggs(as.mcmc(top_mod))

# slopes
jpeg("Density/Figures/DensityGroundwater_ParameterDotPlot_Betas.jpg", units = "in", width = 5, height = 4, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^beta", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = colnames(Covs[,-1]))
dev.off()

# alpha - stream offset
jpeg("Density/Figures/DensityGroundwater_ParameterDotPlot_AlphaStream.jpg", units = "in", width = 5, height = 9, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^alpha.stream.adj", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = StreamNames$stream)
dev.off()

# alpha - section offset
jpeg("Density/Figures/DensityGroundwater_ParameterDotPlot_AlphaSection.jpg", units = "in", width = 5, height = 9, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^alpha.sect.adj", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = unq_sect_str$sectid)
dev.off()


##---------------------------------------------------------------------------------------##
## MARGINAL EFFECTS PLOT
##---------------------------------------------------------------------------------------##


nvalues <- 100 # number of values for prediction
x.gw <- seq(from = min(dat2$gwz), to = max(dat2$gwz), length.out = nvalues) # generate sequences of predictor variables and scales axes
# manually set up and transform axes to original scale of the data
range(dat2$gw)
range(dat2$gwz)
x.axis.gw <- seq(from = 0, to = 1, length.out = 6)
x.scaled.gw <- (x.axis.gw - mean(dat2$gw)) / sd(dat2$gw)

# predictions
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*x.gw + Mcmcdat[i,"beta.2."]*(x.gw^2) + Mcmcdat[i,"beta.3."]*(x.gw^3) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)

# plot
jpeg(filename = "Density/Figures/DensityByGroundwater_MarginalEffect.jpg", height = 3.75, width = 3.75, units = "in", res = 1500)
par(mar = c(3.5,3.5,0.5,1), mgp = c(2.2,0.8,0))
plot(density_log ~ gwz, dat2, pch = NA, xlab = "Groundwater index", ylab = "ln(Peak density, fish/meter)", bty = "l", axes = F)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
polygon(c(x.gw, rev(x.gw)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty = 0)
lines(pred_median ~ x.gw, col = "black", lwd = 2)
points(density_log ~ gwz, dat2, pch = as.numeric(dat2$yrid)+20, bg = dat2$cols)
legend("bottomright", legend = c(2021,2022), pch = c(1,0), bty = "n")
dev.off()


##---------------------------------------------------------------------------------------##
## LOESS method
##---------------------------------------------------------------------------------------##
# polynomials really do not fit the data well, thus just use LOESS smoother b/c this is really not a primary objective


# fit the loess
denlo <- loess(density_log ~ gwz, dat2)

# order the data
j <- order(dat2$gwz)

# predict (with standard errors) from the LOESS object 
denlo_pred <- predict(denlo, se = TRUE)

# plot
jpeg(filename = "Density/Figures/DensityByGroundwater_LOESS.jpg", height = 3.75, width = 3.75, units = "in", res = 1500)
par(mar = c(3.5,3.5,0.5,1), mgp = c(2.2,0.8,0))
plot(density_log ~ gwz, dat2, pch = NA, xlab = "Groundwater index", ylab = "ln(Peak density, fish/meter)", bty = "l", axes = F)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
polygon(c(dat2$gwz[j], rev(dat2$gwz[j])), 
        c(denlo_pred$fit[j] - qt(0.975, denlo_pred$df)*denlo_pred$se.fit[j], rev(denlo_pred$fit[j] + qt(0.975, denlo_pred$df)*denlo_pred$se.fit[j])), 
        col = scales::alpha("black", 0.2), lty = 0)
lines(denlo_pred$fit[j] ~ dat2$gwz[j], col = "black", lwd = 2)
points(density_log ~ gwz, dat2, pch = as.numeric(dat2$yrid)+20, bg = dat2$cols)
legend("bottomright", legend = c(2021,2022), pch = c(1,0), bty = "n")
dev.off()

