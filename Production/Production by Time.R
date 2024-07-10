#==========================================================================================================================================#
# PROJECT NAME: Groundwater structures fish growth and production across a riverscape
# CREATOR: Jeffrey R. Baldock, WY Coop Unit, University of Wyoming, jbaldock@uwyo.edu 
# PURPOSE: Estimate effects of groundwater on temporal trends in YOY production, 
# NOTES: 
#==========================================================================================================================================#

library(lme4)
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


##---------------------------------------------------------------------------------------##
## Data
##---------------------------------------------------------------------------------------##

dat <- read_csv("Growth_DataTable_WithCovariates.csv") %>% mutate(gw = springprev_iew05km_norm) %>% 
  add_row(stream = c("dum1","dum2"), gw = c(0,1)) # add dummy rows to ensure right color scheme

# assign color based on gw influence: purple = snowmelt, yellow = spring-fed
vrPal <- colorRampPalette(rev(viridis(12)))
coldf <- dat %>% group_by(stream, section) %>% 
  summarize(gw = unique(gw)) %>% ungroup() %>%
  arrange(gw) %>%
  mutate(site = paste(stream, section, sep = "_"),
         cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))])
coldf <- coldf %>% filter(!stream %in% c("dum1","dum2"))

# view barplot...does ordering of gw influence make sense?
par(mar = c(5,10,2,2))
barplot(coldf$gw, col = coldf$cols, names.arg = coldf$site, horiz = TRUE, las = 2)
dev.off()

# join colors to data
dat <- dat %>% left_join(coldf)

# days elapsed between sampling intervals
dat$elapdays <- yday(dat$date) - yday(dat$priordate)

# calculate length and weight-based production
dat$prodlen <- ((dat$lenmean - dat$priorlenmean) * dat$density) / dat$elapdays
dat$prodwt <- ((dat$wtmean - dat$priorwtmean) * dat$density) / dat$elapdays

boxplot(dat$prodlen)
abline(h = -0.1)
boxplot(dat$prodwt)
abline(h = -0.01)

# censor to remove negative outliers. Do not censor based on catch b/c low catches are important!
# dat <- dat %>% mutate(prodlen_censor = ifelse(catch >= 5 & priorcatch >= 5, prodlen, NA),
#                       prodwt_censor = ifelse(catch >= 5 & priorcatch >= 5, prodwt, NA))
dat <- dat %>% mutate(prodlen_censor = ifelse(prodlen > -0.013, prodlen, NA),
                      prodwt_censor = ifelse(prodwt > -0.01, prodwt, NA))

sum(!is.na(dat$prodlen_censor))
sum(!is.na(dat$prodwt_censor))
min(dat$prodlen_censor, na.rm = T)
plot(log(prodlen_censor + 0.015) ~ doy, dat, col = dat$cols)
plot(log(prodlen_censor + 0.015) ~ gw, dat, pch = 16, col = dat$cols)

# drop missing data z-score covariates 
dat2 <- dat %>% filter(!is.na(prodlen_censor)) %>% mutate(zgw = c(scale(gw, center = T, scale = T)),
                                                          zdoy = c(scale(doy, center = T, scale = T)))

# create table of covariate means and stdevs for back calculating/plotting
scaletbl <- data.frame(var = c("gw", "doy"),
                       means = c(mean(dat2$gw), mean(dat2$doy)),
                       sds = c(sd(dat2$gw), sd(dat2$doy)))
write_csv(scaletbl, "Production by Time/ProductionByTime_z-score_table.csv")


# pairs plots
boxplot(dat2 %>% dplyr::select(zgw, zdoy))
ggpairs(dat2 %>% dplyr::select(zgw, zdoy))

mmmm <- lm(log(prodlen_censor + 0.015) ~ zdoy + I(zdoy^2) + zgw + zdoy*zgw + I(zdoy^2)*zgw, data = dat2)
plot(mmmm)
plot(log(prodlen_censor + 0.015) ~ zdoy, dat2, col = dat2$cols, pch = 16)

##------------------------------------------------##
## Misc JAGS fields
##------------------------------------------------##

# grouping variables as factors for JAGS model
dat2$stream <- as.factor(dat2$stream)
dat2$sectid <- as.factor(dat2$sectid)
dat2$strid <- as.factor(dat2$strid)
dat2$yrid <- as.factor(dat2$yrid)

# tibble linking each section to each stream for JAGS model with section nested within stream
unq_sect_str <- dat2 %>% group_by(strid) %>% reframe(sectid = unique(sectid)) %>% ungroup()

# some additional data for JAGS
nYears <- length(unique(dat2$yrid))
nSections <- length(unique(dat2$sectid))
nStreams <- length(unique(dat2$strid))

#
StreamNames <- dat2 %>% group_by(stream) %>% summarize(strid = unique(strid), gw = mean(gw)) %>% arrange(gw) %>% ungroup() %>% mutate(strid2 = 1:13) %>% 
  add_row(stream = c("dum1","dum2"), gw = c(0,1)) %>% 
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>%
  filter(!stream %in% c("dum1","dum2"))

# Sort 
dat2 <- dat2 %>% arrange(yrid, strid, sectid)

# drop outlier as it affects model convergence and results in LOO diagnostics indicating model mispecification
# dat2 <- dat2[-140,]

##---------------------------------------------------------------------------------------##
## Fit JAGS Model
##---------------------------------------------------------------------------------------##

# Given our study design (sampling shortly after/during emergence, short sampling intervals and thus overlapping size distributions among sampling events),
# much of the (bootstrapped) variation in growth rates likely stems from variability in emergence timing rather than individual variation in growth rates.
# B/c of overlapping distributions, bootstrapped measures of uncertainty are very large. Therefore, accounting for bootstrapped measures of 
# uncertainty with a state-space model formulation likely results in overly conservative estimates of the effects of covariates on growth, which may limit
# our ability to draw inferences about covariate effects AND explore residual variation in growth among streams that is not already account for in the model.
# Thus, running a model on the means (i.e., no state-space structure) is a justifiable path forward. 

# Length-based production by time model

Covs <- model.matrix(~ zdoy + I(zdoy^2) + zgw + zdoy*zgw + I(zdoy^2)*zgw, data = dat2, na.rm = FALSE)
Covs_Sigma <- model.matrix(~ zgw, data = dat2, na.rm = FALSE) 
jags.data <- list("nObs" = dim(Covs)[1], "P" = log(dat2$prodlen_censor + 0.015), "elap" = dat2$elapdays,                                                                       
                  "strid" = dat2$strid, "yrid" = dat2$yrid, "sectid" = dat2$sectid, "sect_per_str" = unq_sect_str$strid,                                                                        
                  "nYears" = length(unique(dat2$yrid)), "nStreams" = length(unique(dat2$strid)), "nSections" = length(unique(dat2$sectid)),                                                
                  "Covs" = Covs[,-1], "nCovs" = dim(Covs)[2]-1, 
                  "CovsMu" = apply(Covs[,-1], MARGIN = 2, FUN = mean, na.rm = T), "CovsSD" = apply(Covs[,-1], MARGIN = 2, FUN = sd, na.rm = T),
                  "Covs_Sigma" = matrix(Covs_Sigma[,-1]), "nSigmaCovs" = dim(Covs_Sigma)[2]-1)
params <- c("alpha.adj", "alpha.year.adj", "alpha.stream.adj", "alpha.sect.adj", "beta", "mu", "P",
            "sigma.pr", "log.sigma.alpha", "sigma.alpha", "log.sigma.beta", "sigma.sb", 
            "sigma.yr", "sigma.str", "sigma.sect", "sigma.a", "sigma.b", 
            "loglik", "res", "var_fit", "var_fix", "var_res", "var_yr", "var_sect", "var_str", "margR2", "condR2")
prodmod <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "JAGS Models/ProductionModel.txt", 
                         n.chains = 3, n.thin = 20, n.burnin = 50000, n.iter = 150000, DIC = TRUE)
MCMCtrace(prodmod, ind = TRUE, filename = "Production by Time/Traceplots/Production_Length_m_traceplots.pdf") # traceplots
prodmod$BUGSoutput$summary[,8][prodmod$BUGSoutput$summary[,8] > 1.01] # any R-hat values too high?  


##---------------------------------------------------------------------------------------##
## MODEL DIAGNOSTICS/SELECTION - calculate and plot LOO
##---------------------------------------------------------------------------------------##

# LOO-CV to evaluate appropriate model specification
loglik.pm <- prodmod$BUGSoutput$sims.list$loglik
reff.pm <- relative_eff(exp(loglik.pm), chain_id = c(rep(1,5000),rep(2,5000),rep(3,5000)))
loo_pm <- loo(loglik.pm, r_eff = reff.pm)
plot(loo_pm, label_points = TRUE)

# set the top model based on LOO
top_mod <- prodmod


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
write_csv(as.data.frame(Mcmcdat), "Production by Time/ProductionByTime_Length_mcmcsamps.csv")
Mcmcdat <- read.csv("Production by Time/ProductionByTime_Length_mcmcsamps.csv")
write.csv(as.data.frame(param.summary), "Production by Time/ProductionByTime_Length_ParameterSummary.csv", row.names = T)
param.summary_len <- read.csv("Production by Time/ProductionByTime_Length_ParameterSummary.csv", row.names = 1)


##---------------------------------------------------------------------------------------##
## MODEL DIAGNOSTICS
##---------------------------------------------------------------------------------------##

# subset expected and observed MCMC samples
ppdat_exp <- as.matrix(Mcmcdat[,startsWith(names(Mcmcdat), "mu.")])
ppdat_obs <- as.matrix(Mcmcdat[,startsWith(names(Mcmcdat), "P.")])

# Bayesian p-value
sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2])

# Posterior Predictive Check: plot median posterior expected length ~ observed length with Bayesian p-value
jpeg("Production by Time/Figures/ProductionByTime_PPCheck.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(x = seq(from = -6, to = 0.5, length.out = 100), y = seq(from = -6, to = 0.5, length.out = 100), pch = NA, xlab = "Observed Production (mm per m per day)", ylab = "Median Posterior Expected Production")
for (i in 1:dim(Covs)[1]) { points(median(Mcmcdat[,paste("mu.",i,".", sep = "")]) ~ log(dat2$prodlen_censor[i] + 0.015))}
legend("topleft", bty = "n", legend = paste("Bayesian p-value = ", round(sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2]), digits = 3), sep = ""))
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()

# histogram of residuals
jpeg("Production by Time/Figures/ProductionByTime_ResidHist.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist((ppdat_obs - ppdat_exp), main = "", xlab = "Observed - Expected")
legend("topleft", bty = "n", legend = paste("Median Resid. = ", round(median((ppdat_obs - ppdat_exp)), digits = 2), sep = ""))
abline(v = median(unlist(ppdat_obs - ppdat_exp)), col = "red", lwd = 2)
dev.off()

# RMSE
sqrt(mean((ppdat_obs - ppdat_exp)^2)) # point estimate
rmse <- vector("numeric", 15000L)
for (j in 1:dim(Mcmcdat)[1]) {
  res <- vector("numeric", 258L)
  for (i in 1:dim(Covs)[1]) {
    res[i] <- (ppdat_obs[j,i] - ppdat_exp[j,i])^2
  }
  rmse[j] <- sqrt(mean(res))
  print(j)
}
jpeg("Production by Time/Figures/ProductionByTime_RMSE.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist(rmse, xlab = "Root Mean Square Error", ylab = "Posterior Frequency", main = "")
legend("topright", bty = "n", legend = paste("Median RMSE = ", round(median(rmse), digits = 2), " mm", sep = ""))
abline(v = median(rmse), col = "red", lwd = 2)
dev.off()


# Calculate Bayesian R^2 per Gelman et al. (2018), Nakagawa and Shielzeth (2013)

var_fix <- Mcmcdat[,"var_fix"]
var_fit <- Mcmcdat[,"var_fit"]
var_str <- Mcmcdat[,"var_str"]
var_yr <- Mcmcdat[,"var_yr"]
var_sect <- Mcmcdat[,"var_sect"]
var_res <- Mcmcdat[,"var_res"]

# model derived marginal and condition R2
hist(Mcmcdat[,"margR2"])
hist(Mcmcdat[,"condR2"])

# conditional R2 calculated in two ways
hist((var_fix + var_yr + var_str + var_sect) / (var_fix + var_res + var_yr + var_str + var_sect)) # "manually"
hist(var_fit / (var_fit + var_res)) # per Gelman et al. (2019). Note: this is ~identical to ouput using BRMS and performance packages, function = r2_posterior()
median((var_fix + var_yr + var_str + var_sect) / (var_fix + var_res + var_yr + var_str + var_sect))
median(var_fit / (var_fit + var_res))

# dropping year component
hist((var_fix + var_str + var_sect) / (var_fix + var_res + var_str + var_sect))
median((var_fix + var_str + var_sect) / (var_fix + var_res + var_str + var_sect))
hist((var_fix) / (var_fix + var_res + var_str + var_sect))
median((var_fix) / (var_fix + var_res + var_str + var_sect))

# random components each
hist((var_str) / (var_fix + var_str + var_yr + var_sect + var_res))
hist((var_sect) / (var_fix + var_str + var_yr + var_sect + var_res))
hist((var_yr) / (var_fix + var_str + var_yr + var_sect + var_res))
### Very broad distribution of variation explained by year component likely driven by few (2) factors associated with that hierarchical component.
### Multiple references (Gelman and Hill???) note that hierarchical variance may be poorly estimated (overly broad) when number of factors is small. 
### This results in distributions of marginal and conditional R^2 that are also very broad. Thus, should rely on modes for inference, as the mode 
### (more so than the median) reduces the influence of overly broad year component. Modes align well with preliminary analysis using frequentist
### approach (lmer), which yields marginal R2 = 0.401 and conditional R2 = 0.495


# Marginal/conditional plot
# Note, using BRMS and r2_posterior() in "performance" package yields the same results as those presented here
jpeg("Production by Time/Figures/ProductionByTime_MargCondR2.jpg", units = "in", width = 5.5, height = 5, res = 1500)
par(mar = c(4,4,0.5,0.5), mgp = c(2.5,1,0))
hist((var_fix) / (var_fix + var_res + var_str + var_sect), col = alpha("darkorange", 0.5), xlim = c(0,1), ylim = c(0,5500), breaks = seq(from = 0, to = 1, by = 0.02), main = NA, xlab = expression(Bayesian ~ R^2), ylab = "Posterior Frequency")
hist(var_fit / (var_fit + var_res), col = alpha("forestgreen", 0.5), xlim = c(0,1), add = TRUE, breaks = seq(from = 0, to = 1, by = 0.02))
# hist((var_fix + var_str + var_sect) / (var_fix + var_res + var_str + var_sect), col = alpha("dodgerblue", 0.5), xlim = c(0,1), add = TRUE, breaks = seq(from = 0, to = 1, by = 0.02))
legend("topleft", fill = c(alpha("darkorange", 0.5), alpha("forestgreen", 0.5)), bty = "n", cex = 0.9,
       legend = c(paste("Marginal = ", round(median((var_fix) / (var_fix + var_res + var_str + var_sect)), 3), " (", round(hdi((var_fix) / (var_fix + var_res + var_str + var_sect), credMass = 0.95)[1], 3), " - ", round(hdi((var_fix) / (var_fix + var_res + var_str + var_sect), credMass = 0.95)[2], 3), ")", sep = ""), 
                  paste("Conditional = ", round(median(var_fit / (var_fit + var_res)), 3), " (", round(hdi(var_fit / (var_fit + var_res), credMass = 0.95)[1], 3), " - ", round(hdi(var_fit / (var_fit + var_res), credMass = 0.95)[2], 3), ")", sep = "")))
dev.off()


##### Combined Figure
jpeg("Production by Time/Figures/ProductionByTime_ModelDiagnosticsCombined.jpg", units = "in", width = 7, height = 7, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), mfrow = c(2,2))

# PP Check
plot(x = seq(from = -6, to = 0.5, length.out = 100), y = seq(from = -6, to = 0.5, length.out = 100), pch = NA, xlab = expression(paste("Observed Production (mm m"^"-1", "day"^"-1", ", logged)", sep = "")), ylab = "Median Posterior Expected Production", bty = "l")
for (i in 1:dim(Covs)[1]) { points(median(Mcmcdat[,paste("mu.",i,".", sep = "")]) ~ log(dat2$prodlen_censor[i] + 0.015))}
legend("top", bty = "n", legend = paste("Bayesian p-value = ", round(sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2]), digits = 3), sep = ""))
abline(a = 0, b = 1, col = "red", lwd = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(a)")
par(usr = usr)

# RMSE
hist(rmse, xlab = "Root Mean Square Error", ylab = "Posterior Frequency", main = "")
box(bty = "l")
legend("topright", bty = "n", legend = paste("Median = ", round(median(rmse), digits = 3), sep = ""))
abline(v = median(rmse), col = "red", lwd = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(b)")
par(usr = usr)

# Marg/Cond R2
hist((var_fix) / (var_fix + var_res + var_str + var_sect), col = alpha("darkorange", 0.5), xlim = c(0,1), ylim = c(0,6000), breaks = seq(from = 0, to = 1, by = 0.02), main = NA, xlab = expression(Bayesian ~ R^2), ylab = "Posterior Frequency")
hist(var_fit / (var_fit + var_res), col = alpha("forestgreen", 0.5), xlim = c(0,1), add = TRUE, breaks = seq(from = 0, to = 1, by = 0.02))
legend("topright", fill = c(alpha("darkorange", 0.5), alpha("forestgreen", 0.5)), bty = "n", cex = 0.9,
       legend = c(paste("Marginal = ", round(median((var_fix) / (var_fix + var_res + var_str + var_sect)), 3), " (", round(hdi((var_fix) / (var_fix + var_res + var_str + var_sect), credMass = 0.95)[1], 3), " - ", round(hdi((var_fix) / (var_fix + var_res + var_str + var_sect), credMass = 0.95)[2], 3), ")", sep = ""), 
                  paste("Conditional = ", round(median(var_fit / (var_fit + var_res)), 3), " (", round(hdi(var_fit / (var_fit + var_res), credMass = 0.95)[1], 3), " - ", round(hdi(var_fit / (var_fit + var_res), credMass = 0.95)[2], 3), ")", sep = "")))
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(c)")
par(usr = usr)

# residuals by Groundwater
plot(apply(top_mod$BUGSoutput$sims.list$res, 2, mean) ~ dat2$gw, xlab = "Groundwater Index", ylab = "Residuals", bty = "l")
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(d)")
par(usr = usr)

dev.off()


##---------------------------------------------------------------------------------------##
## PLOTTING
##---------------------------------------------------------------------------------------##

mod.gg <- ggs(as.mcmc(top_mod))
Covs <- model.matrix(~ zdoy + I(zdoy^2) + zgw + zdoy*zgw + I(zdoy^2)*zgw, data = dat2, na.rm = FALSE)

# Parameter (beta) estimates - dot plot
jpeg("Production by Time/Figures/ProductionByTime_ParameterDotPlot.jpg", units = "in", width = 5, height = 4, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^beta", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = colnames(Covs)[-1])
dev.off()

# Parameter (log.sigma.beta: covariates on uncertainty) estimates - dot plot
jpeg("Production by Time/Figures/ProductionByTime_ParameterDotPlot_SigmaBetas.jpg", units = "in", width = 5, height = 4, res = 1500)
ggs_caterpillar(D = mod.gg, family = "log.sigma.beta", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = "zgw")
dev.off()

# alpha.stream - stream offsets to the intercept
StreamNames <- StreamNames %>% arrange(strid)
jpeg("Production by Time/Figures/ProductionByTime_ParameterDotPlot_alpha.stream.jpg", units = "in", width = 5, height = 4, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^alpha.stream.adj", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) +
  theme_bw() + ylab("") + xlab("Posterior estimate") +
  scale_y_discrete(labels = StreamNames$stream)
dev.off()

# alpha.section - section offsets to the intercept
jpeg("Production by Time/Figures/ProductionByTime_ParameterDotPlot_alpha.section.jpg", units = "in", width = 5, height = 9, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^alpha.sect.adj", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = unq_sect_str$sectid)
dev.off()

# alpha.year - year offsets to the intercept
jpeg("Production by Time/Figures/ProductionByTime_ParameterDotPlot_alpha.year.jpg", units = "in", width = 5, height = 9, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^alpha.year.adj", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = c(2021:2022))
dev.off()

#-------------------------------------------#
# Random Intercepts: Stream x Year
#-------------------------------------------#
StreamNames <- StreamNames %>% arrange(strid2)

jpeg("Production by Time/Figures/ProductionByTime_PosteriorIntercepts.jpg", units = "in", width = 3.75, height = 6, res = 1500)
par(mar = c(4,3,1.5,1), mgp = c(2.5,1,0))
# set up plot area
x <- seq(from = -3.5, to = -1.5, length.out = 100)
y <- seq(from = 1, to = 13.7, length.out = 100)
plot(y ~ x, type = "n", xlab = expression(paste("Mean ln(Production, mm m"^"-1", "day"^"-1", ")", sep = "")), ylab = "", axes = F)
axis(1)
# iterate over intercepts
ctr <- 0
for (j in nStreams:1) {
  par(xpd = NA)
  text(y = j+0.5, x = -4.0, labels = StreamNames$stream[j], pos = 4)
  par(xpd = FALSE)
  dens <- density(Mcmcdat[,"alpha.adj"] + Mcmcdat[,paste("alpha.stream.adj.",StreamNames$strid[j],".", sep = "")])
  dens$y2 <- (dens$y / max(dens$y))
  l <- min(which(dens$x >= hdi(dens, credMass = 0.95)[1]))
  h <- max(which(dens$x < hdi(dens, credMass = 0.95)[2]))
  polygon(c(dens$x[c(l, l:h, h)]), c(0+j, dens$y2[l:h]+j, 0+j), col = adjustcolor(StreamNames$cols[j], alpha.f = 0.8), lty = 0)
  lines(dens$y2+j ~ dens$x, lwd = 1.5)
  ctr <- ctr + 4
}
abline(v = median(Mcmcdat[,"alpha.adj"]), lty = 2)
par(xpd = TRUE)
dev.off()


jpeg("Production by Time/Figures/ProductionByTime_PosteriorIntercepts_withYears.jpg", units = "in", width = 3.75, height = 6, res = 1500)
par(mar = c(4,3,1.5,1), mgp = c(2.5,1,0))
# set up plot area
x <- seq(from = -3.75, to = -1.5, length.out = 100)
y <- seq(from = 1, to = 13.7, length.out = 100)
plot(y ~ x, type = "n", xlab = expression(paste("Mean ln(Production, mm m"^"-1", "day"^"-1", ")", sep = "")), ylab = "", axes = F)
axis(1)
# iterate over intercepts
ctr <- 0
for (j in nStreams:1) {
  par(xpd = NA)
  text(y = j+0.5, x = -4.25, labels = StreamNames$stream[j], pos = 4)
  for (i in 1:nYears){
    par(xpd = FALSE)
    dens <- density(Mcmcdat[,"alpha.adj"] + Mcmcdat[,paste("alpha.year.adj.",i,".", sep = "")] + Mcmcdat[,paste("alpha.stream.adj.",StreamNames$strid[j],".", sep = "")])
    dens$y2 <- (dens$y / max(dens$y))
    l <- min(which(dens$x >= hdi(dens, credMass = 0.95)[1]))
    h <- max(which(dens$x < hdi(dens, credMass = 0.95)[2]))
    polygon(c(dens$x[c(l, l:h, h)]), c(0+j, dens$y2[l:h]+j, 0+j), col = adjustcolor(StreamNames$cols[j], alpha.f = 0.5), lty = 0)
    lines(dens$y2+j ~ dens$x, lwd = 1.5, lty = i)
  }
  ctr <- ctr + 4
}
abline(v = median(Mcmcdat[,"alpha.adj"]))
par(xpd = TRUE)
legend(x = -3.5, y = 15.2, legend = c(2021, 2022), lty = 1:2, lwd = 1.5, bty = "n", horiz = TRUE, cex = 0.9)
dev.off()


#-------------------------------------------#
# Marginal Effects Plots
#-------------------------------------------#

# set up gradient legend
gradpal <- tibble(gw = seq(from = 0, to = 1, length.out = 100)) %>%
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>% 
  filter(gw >= min(dat2$gw) & gw <= max(dat2$gw))
# number of values for prediction
nvalues <- 100
# generate sequences of predictor variables and scales axes
x.gw <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = nvalues)
x.doy <- seq(from = min(dat2$zdoy), to = max(dat2$zdoy), length.out = nvalues)
# manually set up and transform axes to original scale of the data
x.axis.gw <- seq(from = 0, to = 0.8, by = 0.1)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
x.axis.doy <- c(213,244,274,305)
x.scaled.doy <- (x.axis.doy - scaletbl$means[2]) / scaletbl$sds[2]

# Model without prior length
jpeg("Production by Time/Figures/ProductionByTime_TimeGWEffect_wPoints.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# growth ~ temp + density (temp focal)
plot(log(prodlen_censor + 0.015) ~ zdoy, dat2, pch = NA, ylim = c(-6.5,0.5), xlab = "", ylab = expression(paste("Production (mm m"^"-1", "day"^"-1", ", log scale)", sep = "")), axes = F)
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2, at = log(c(-0.01, 0.01, 0.1, 1) + 0.015), labels = c(-0.01, 0.01, 0.1, 1))
box(bty = "l")
abline(h = log(0.015), lty = 2)
# interaction: doy * high gw
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*x.doy + Mcmcdat[i,"beta.2."]*(x.doy^2) + Mcmcdat[i,"beta.3."]*max(x.gw) + Mcmcdat[i,"beta.4."]*max(x.gw)*x.doy + Mcmcdat[i,"beta.5."]*max(x.gw)*(x.doy^2) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == max(coldf$gw), "cols"], 0.2), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == max(coldf$gw), "cols"]), lwd = 2)
# interaction: doy * low gw
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*x.doy + Mcmcdat[i,"beta.2."]*(x.doy^2) + Mcmcdat[i,"beta.3."]*min(x.gw) + Mcmcdat[i,"beta.4."]*min(x.gw)*x.doy + Mcmcdat[i,"beta.5."]*min(x.gw)*(x.doy^2) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == min(coldf$gw), "cols"], 0.35), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == min(coldf$gw), "cols"]), lwd = 2)
# points
points(log(prodlen_censor + 0.015) ~ zdoy, dat2, pch = 16, col = alpha(dat2$cols, 0.5))
# legend
# text(x = 1.2, y = 4.0, labels = "Groundwater\nInfluence", cex = 0.9)
gradientLegend(valRange = c(min(gradpal$gw), max(gradpal$gw)), color = gradpal$cols, length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.05, 0.03, 0.1, 0.27), dec = 2)
dev.off()


#-------------------------------------------#
# (Cumulative) Production ~ time ~ groundwater
#-------------------------------------------#

# number of groundwater values
nvalues <- 100
# numver of MCMC draws
nsims <- dim(Mcmcdat)[1]
# to project for complete range of gw influence, need to manually set min/max values
xmin <- (0 - scaletbl$means[1]) / scaletbl$sds[1]
xmax <- (1 - scaletbl$means[1]) / scaletbl$sds[1]
doys <- c(min(dat2$doy, na.rm = T):max(dat2$doy, na.rm = T))
# generate sequences of predictor variables and scales axes
x.gw <- seq(from = xmin, to = xmax, length.out = nvalues)
x.doy <- (doys - scaletbl$means[2]) / scaletbl$sds[2]
# manually set up and transform axes to original scale of the data
x.axis.gw <- seq(from = 0, to = 1, by = 0.2)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
x.axis.doy <- c(213,244,274,305)
x.scaled.doy <- (x.axis.doy - scaletbl$means[2]) / scaletbl$sds[2]
# colors
colss <- rev(hcl.colors(n = nvalues, "Viridis"))

# project production with uncertainty
pred_dist_mat <- array(data = NA, dim = c(nsims, length(x.doy), nvalues))
for (i in 1:nvalues) {
  for (j in 1:nsims) { pred_dist_mat[j,,i] <- Mcmcdat[j,"alpha.adj"] + Mcmcdat[j,"beta.1."]*x.doy + Mcmcdat[j,"beta.2."]*(x.doy^2) + Mcmcdat[j,"beta.3."]*x.gw[i] + Mcmcdat[j,"beta.4."]*x.gw[i]*x.doy + Mcmcdat[j,"beta.5."]*x.gw[i]*(x.doy^2) }
  print(i)
}
# median
pred_med_mat <- matrix(data = NA, nrow = length(x.doy), ncol = nvalues)
for (i in 1:length(x.doy)) {
  for (j in 1:nvalues) {
    pred_med_mat[i,j] <- median(pred_dist_mat[,i,j])
  }
}
# end of season cumulative production by gw
pred_eos_mat <- matrix(data = NA, nrow = nsims, ncol = nvalues)
for (i in 1:nsims) {
  for (j in 1:nvalues) {
    pred_eos_mat[i,j] <- sum(exp(pred_dist_mat[i,,j]))
  }
  print(i)
}
write_csv(as_tibble(pred_med_mat), "Production by Time/ProductionByTime_ProjectedDailyLogProduction_Median.csv")
write_csv(as_tibble(exp(pred_med_mat)), "Production by Time/ProductionByTime_ProjectedDailyProduction_Median.csv")
write_csv(as_tibble(apply(exp(pred_med_mat), 2, cumsum)), "Production by Time/ProductionByTime_ProjectedCumulativeProduction_Median.csv")

pred_med_mat <- as.matrix(read_csv("Production by Time/ProductionByTime_ProjectedDailyLogProduction_Median.csv"))


# plot cumulative production through time
jpeg("Production by Time/Figures/ProductionByTime_CumulativeProduction_ThroughTime.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(seq(from = 0, to = 22, length.out = length(doys)) ~ x.doy, dat, pch = NA, xlab = "", ylab = expression(paste("Cumulative Production (mm m"^"-1", ")", sep = "")), axes = F)
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "l")
for (i in 1:length(x.gw)) { lines(cumsum(exp(pred_med_mat[,i])) ~ x.doy, lwd = 5, col = colss[i]) }
# add lines for min/max observed gw influence
dum_vec <- matrix(data = NA, nrow = length(x.doy), ncol = 2)
for (j in 1:length(x.doy)) {
  dum_vec[j,1] <- median(Mcmcdat[,"alpha.adj"]) + median(Mcmcdat[,"beta.1."])*x.doy[j] + median(Mcmcdat[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat[,"beta.3."])*min(dat2$zgw) + median(Mcmcdat[,"beta.4."])*min(dat2$zgw)*x.doy[j] + median(Mcmcdat[,"beta.5."])*min(dat2$zgw)*(x.doy[j]^2)
  dum_vec[j,2] <- median(Mcmcdat[,"alpha.adj"]) + median(Mcmcdat[,"beta.1."])*x.doy[j] + median(Mcmcdat[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat[,"beta.3."])*max(dat2$zgw) + median(Mcmcdat[,"beta.4."])*max(dat2$zgw)*x.doy[j] + median(Mcmcdat[,"beta.5."])*max(dat2$zgw)*(x.doy[j]^2)
}
for (i in 1:2) { lines(cumsum(exp(dum_vec[,i])) ~ x.doy, lty = 1, lwd = 0.75, col = "white") }
# text(x = 0, y = 10, labels = "Groundwater \nInfluence", cex = 0.8)
gradientLegend(valRange = c(0,1), color = rev(viridis(100)), length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.07, 0.65, 0.14, 1), dec = 2)
par(usr = c(0,1,0,1))
lines(x = c(0.072, 0.138), y = c((min(dat2$gw)*(1-0.65))+0.65, (min(dat2$gw)*(1-0.65))+0.65), lty = 1, lwd = 0.75, col = "white")
lines(x = c(0.072, 0.138), y = c((max(dat2$gw)*(1-0.65))+0.65, (max(dat2$gw)*(1-0.65))+0.65), lty = 1, lwd = 0.75, col = "white")
dev.off()


# plot end of season cumulative production
ddd <- dat2 %>% mutate(prodlenraw = (lenmean - priorlenmean) * density) %>% group_by(year, stream, section, zgw) %>% summarize(n = n(), cumprod = sum(prodlenraw))

jpeg("Production by Time/Figures/ProductionByTime_CumulativeProduction_EndOfSeason_WithPoints.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(seq(from = 0, to = 80, length.out = nvalues) ~ x.gw, dat, pch = NA, xlab = "Groundwater Influence", ylab = expression(paste("Projected cumulative production (mm m"^"-1", ")", sep = "")), axes = F)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
polygon(x = c(x.gw, rev(x.gw)), y = c(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.025), rev(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.975))), col = alpha("black", 0.2), lty = 0)
lines(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.5) ~ x.gw, lwd = 2)
points(cumprod ~ zgw, ddd)
legend("topleft", legend = c("Observed", "Projected"), pch = c(1,NA), lty = c(NA,1), lwd = c(NA,2), bty = "n")
dev.off()

jpeg("Production by Time/Figures/ProductionByTime_CumulativeProduction_EndOfSeason_NoPoints.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(seq(from = 0, to = 53, length.out = nvalues) ~ x.gw, dat, pch = NA, xlab = "Groundwater Influence", ylab = expression(paste("Projected cumulative production (mm m"^"-1", ")", sep = "")), axes = F)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
polygon(x = c(x.gw, rev(x.gw)), y = c(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.025), rev(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.975))), col = alpha("black", 0.2), lty = 0)
lines(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.5) ~ x.gw, lwd = 2)
dev.off()

# plot(log(cumprod) ~ zgw, ddd)
# 
# plot(seq(from = -3, to = 5, length.out = nvalues) ~ x.gw, dat, pch = NA, xlab = "Groundwater Influence", ylab = expression(paste("Projected cumulative production (mm m"^"-1", ")", sep = "")), axes = F)
# axis(1, at = x.scaled.gw, labels = x.axis.gw)
# axis(2)
# box(bty = "l")
# polygon(x = c(x.gw, rev(x.gw)), y = log(c(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.025), rev(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.975)))), col = alpha("black", 0.2), lty = 0)
# lines(log(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.5)) ~ x.gw, lwd = 2)
# points(log(cumprod) ~ zgw, ddd)


#-------------------------------------------#
# Groundwater effect on sigma
#-------------------------------------------#

# number of values for prediction
nvalues <- 100
# generate sequences of predictor variables and scales axes
x.gw <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = nvalues)
x.axis.gw <- seq(from = 0, to = 0.8, by = 0.2)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]

jpeg("Production by Time/Figures/ProductionByTime_GWEffectOnSigma.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# growth ~ temp + density (temp focal)
plot(seq(from = 0.025, to = 0.041, length.out = length(dat2$zdoy)) ~ zgw, dat2, pch = NA, xlab = "Groundwater index", ylab = "Standard deviation of ln(Production)", axes = F)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- exp(Mcmcdat[i,"log.sigma.alpha"] + Mcmcdat[i,"log.sigma.beta"]*x.gw) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.gw, rev(x.gw)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty = 0)
lines(pred_median ~ x.gw, lwd = 2)
dev.off()


#-------------------------------------------#
# Combined Figure
#-------------------------------------------#

jpeg("Production by Time/Figures/ProductionByTime_CombinedFigure.jpg", units = "in", width = 6.75, height = 6, res = 1500)
par(mar = c(3,5,2,1), mgp = c(2.5,1,0), mfrow = c(2,2))

### Marginal effects: production ~ doy + gw
# set up gradient legend
gradpal <- tibble(gw = seq(from = 0, to = 1, length.out = 100)) %>%
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>% 
  filter(gw >= min(dat2$gw) & gw <= max(dat2$gw))
# number of values for prediction
nvalues <- 100
# generate sequences of predictor variables and scales axes
x.gw <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = nvalues)
x.doy <- seq(from = min(dat2$zdoy), to = max(dat2$zdoy), length.out = nvalues)
# manually set up and transform axes to original scale of the data
x.axis.gw <- seq(from = 0, to = 0.8, by = 0.1)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
x.axis.doy <- c(213,244,274,305)
x.scaled.doy <- (x.axis.doy - scaletbl$means[2]) / scaletbl$sds[2]
# plot
plot(log(prodlen_censor + 0.015) ~ zdoy, dat2, pch = NA, xlab = "", ylab = expression(paste("Production (mm m"^"-1", "day"^"-1", ")", sep = "")), axes = F)
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2, at = log(c(-0.01, 0.01, 0.1, 1) + 0.015), labels = c(-0.01, 0.01, 0.1, 1))
box(bty = "l")
abline(h = log(0.015), lty = 2)
# interaction: doy * high gw
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*x.doy + Mcmcdat[i,"beta.2."]*(x.doy^2) + Mcmcdat[i,"beta.3."]*max(x.gw) + Mcmcdat[i,"beta.4."]*max(x.gw)*x.doy + Mcmcdat[i,"beta.5."]*max(x.gw)*(x.doy^2) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == max(coldf$gw), "cols"], 0.2), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == max(coldf$gw), "cols"]), lwd = 2)
# interaction: doy * low gw
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*x.doy + Mcmcdat[i,"beta.2."]*(x.doy^2) + Mcmcdat[i,"beta.3."]*min(x.gw) + Mcmcdat[i,"beta.4."]*min(x.gw)*x.doy + Mcmcdat[i,"beta.5."]*min(x.gw)*(x.doy^2) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha("white", 0.5), lty=0)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == min(coldf$gw), "cols"], 0.27), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == min(coldf$gw), "cols"]), lwd = 2)
# points
points(log(prodlen_censor + 0.015) ~ zdoy, dat2, pch = 16, col = alpha(dat2$cols, 0.5))
# legend
# gradientLegend(valRange = c(min(gradpal$gw), max(gradpal$gw)), color = gradpal$cols, length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.05, 0.03, 0.1, 0.27), dec = 2)
# panel label
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(a)")
par(usr = usr)


### Groundwater effect on sigma
# number of values for prediction
nvalues <- 100
# generate sequences of predictor variables and scales axes
x.gw <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = nvalues)
x.axis.gw <- seq(from = 0, to = 0.8, by = 0.2)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
# plot
plot(seq(from = 0.025, to = 0.041, length.out = length(dat2$zdoy)) ~ zgw, dat2, pch = NA, xlab = "", ylab = "Standard deviation of ln(Production)", axes = F)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- exp(Mcmcdat[i,"log.sigma.alpha"] + Mcmcdat[i,"log.sigma.beta"]*x.gw) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.gw, rev(x.gw)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty = 0)
lines(pred_median ~ x.gw, lwd = 2)
# panel label
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(b)")
par(usr = usr)


par(mar = c(4,5,1,1))
### cumulative production through time
# number of groundwater values
nvalues <- 100
# numver of MCMC draws
nsims <- dim(Mcmcdat)[1]
# to project for complete range of gw influence, need to manually set min/max values
xmin <- (0 - scaletbl$means[1]) / scaletbl$sds[1]
xmax <- (1 - scaletbl$means[1]) / scaletbl$sds[1]
doys <- c(min(dat2$doy, na.rm = T):max(dat2$doy, na.rm = T))
# generate sequences of predictor variables and scales axes
x.gw <- seq(from = xmin, to = xmax, length.out = nvalues)
x.doy <- (doys - scaletbl$means[2]) / scaletbl$sds[2]
# manually set up and transform axes to original scale of the data
x.axis.gw <- seq(from = 0, to = 1, by = 0.2)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
x.axis.doy <- c(213,244,274,305)
x.scaled.doy <- (x.axis.doy - scaletbl$means[2]) / scaletbl$sds[2]
# colors
colss <- rev(hcl.colors(n = nvalues, "Viridis"))
# plot
plot(seq(from = 0, to = 22, length.out = length(doys)) ~ x.doy, dat, pch = NA, xlab = "Time", ylab = expression(paste("Cumulative production (mm m"^"-1", ")", sep = "")), axes = F)
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "l")
for (i in 1:length(x.gw)) { lines(cumsum(exp(pred_med_mat[,i])) ~ x.doy, lwd = 5, col = colss[i]) }
# add lines for min/max observed gw influence
dum_vec <- matrix(data = NA, nrow = length(x.doy), ncol = 2)
for (j in 1:length(x.doy)) {
  dum_vec[j,1] <- median(Mcmcdat[,"alpha.adj"]) + median(Mcmcdat[,"beta.1."])*x.doy[j] + median(Mcmcdat[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat[,"beta.3."])*min(dat2$zgw) + median(Mcmcdat[,"beta.4."])*min(dat2$zgw)*x.doy[j] + median(Mcmcdat[,"beta.5."])*min(dat2$zgw)*(x.doy[j]^2)
  dum_vec[j,2] <- median(Mcmcdat[,"alpha.adj"]) + median(Mcmcdat[,"beta.1."])*x.doy[j] + median(Mcmcdat[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat[,"beta.3."])*max(dat2$zgw) + median(Mcmcdat[,"beta.4."])*max(dat2$zgw)*x.doy[j] + median(Mcmcdat[,"beta.5."])*max(dat2$zgw)*(x.doy[j]^2)
}
for (i in 1:2) { lines(cumsum(exp(dum_vec[,i])) ~ x.doy, lty = 1, lwd = 0.75, col = "white") }
text(x = -1.1, y = 19, labels = "Groundwater\nindex", cex = 0.8)
gradientLegend(valRange = c(0,1), color = rev(viridis(100)), length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.05, 0.25, 0.12, 0.72), dec = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
lines(x = c(0.052, 0.118), y = c((min(dat2$gw)*(0.72-0.25))+0.25, (min(dat2$gw)*(0.72-0.25))+0.25), lty = 1, lwd = 0.75, col = "white")
lines(x = c(0.052, 0.118), y = c((max(dat2$gw)*(0.72-0.25))+0.25, (max(dat2$gw)*(0.72-0.25))+0.25), lty = 1, lwd = 0.75, col = "white")
text(0.07, 0.95, labels = "(c)")
par(usr = usr)


### projectd end of season cumulatve production
plot(seq(from = 0, to = 53, length.out = nvalues) ~ x.gw, dat, pch = NA, xlab = "Groundwater index", ylab = expression(paste("      Projected end-of-season \ncumulative production (mm m"^"-1", ")", sep = "")), axes = F)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
polygon(x = c(x.gw, rev(x.gw)), y = c(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.025), rev(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.975))), col = alpha("black", 0.2), lty = 0)
lines(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.5) ~ x.gw, lwd = 2)
# panel label
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(d)")
par(usr = usr)

dev.off()

# nvalues <- 100 # number of groundwater values
# nsims <- 100 # numver of MCMC draws
# xmin <- (0 - scaletbl$means[1]) / scaletbl$sds[1]
# xmax <- (1 - scaletbl$means[1]) / scaletbl$sds[1]
# doys <- c(min(dat2$doy, na.rm = T):max(dat2$doy, na.rm = T))
# x.gw <- seq(from = xmin, to = xmax, length.out = nvalues)
# x.doy <- (doys - scaletbl$means[2]) / scaletbl$sds[2]
# x.axis.gw <- seq(from = 0, to = 1, by = 0.2)
# x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
# x.axis.doy <- c(213,244,274,305)
# x.scaled.doy <- (x.axis.doy - scaletbl$means[2]) / scaletbl$sds[2]
# colss <- hcl.colors(n = nvalues, "Viridis")
# plot(seq(from = 0, to = 25, length.out = length(doys)) ~ x.doy, dat, pch = NA, xlab = "", ylab = expression(paste("Cumulative Production (mm m"^"-1", ")", sep = "")), axes = F)
# axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
# axis(2)
# box(bty = "l")
# for (i in 1:length(x.gw)) { lines(cumsum(exp(pred_med_mat[,i])) ~ x.doy, lwd = 5, col = colss[i]) }
# dum_vec <- matrix(data = NA, nrow = length(x.doy), ncol = 2)
# for (j in 1:length(x.doy)) {
#   dum_vec[j,1] <- median(Mcmcdat[,"alpha.adj"]) + median(Mcmcdat[,"beta.1."])*x.doy[j] + median(Mcmcdat[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat[,"beta.3."])*min(dat2$zgw) + median(Mcmcdat[,"beta.4."])*min(dat2$zgw)*x.doy[j] + median(Mcmcdat[,"beta.5."])*min(dat2$zgw)*(x.doy[j]^2)
#   dum_vec[j,2] <- median(Mcmcdat[,"alpha.adj"]) + median(Mcmcdat[,"beta.1."])*x.doy[j] + median(Mcmcdat[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat[,"beta.3."])*max(dat2$zgw) + median(Mcmcdat[,"beta.4."])*max(dat2$zgw)*x.doy[j] + median(Mcmcdat[,"beta.5."])*max(dat2$zgw)*(x.doy[j]^2)
# }
# for (i in 1:2) { lines(cumsum(exp(dum_vec[,i])) ~ x.doy, lty = 3, lwd = 1, col = "white") }
# # text(x = -1.2, y = 24.5, labels = "Groundwater \nInfluence", cex = 0.8)
# gradientLegend(valRange = c(0,1), color = viridis(100), length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.07, 0.5, 0.12, 0.85), dec = 2)
# usr <- par("usr")
# par(usr = c(0,1,0,1))
# text(0.07, 0.95, labels = "(b)")
# par(usr = usr)
# 
# 
# # projected end of season cumulative production by groundwater influence, with credible interval
# plot(seq(from = 0, to = 80, length.out = nvalues) ~ x.gw, dat, pch = NA, xlab = "Groundwater Influence", ylab = expression(paste("Cumulative Production (mm m"^"-1", ")", sep = "")), axes = F)
# axis(1, at = x.scaled.gw, labels = x.axis.gw)
# axis(2)
# box(bty = "l")
# polygon(x = c(x.gw, rev(x.gw)), y = c(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.025), rev(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.975))), col = alpha("black", 0.2), lty = 0)
# lines(apply(pred_eos_mat, MARGIN = 2, quantile, prob = 0.5) ~ x.gw, lwd = 2)
# points(cumprod ~ zgw, ddd)
# legend(x = -1.2, y = 77, legend = c("Observed", "Projected"), pch = c(1,NA), lty = c(NA,1), lwd = c(NA,2), bty = "n", cex = 0.8)
# usr <- par("usr")
# par(usr = c(0,1,0,1))
# text(0.07, 0.95, labels = "(c)")
# par(usr = usr)
# 
# 
# # groundwater effect on sigma
# nvalues <- 100
# x.gw <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = nvalues)
# x.axis.gw <- seq(from = 0, to = 0.8, by = 0.2)
# x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
# plot(seq(from = 0.025, to = 0.04, length.out = length(dat2$zdoy)) ~ zgw, dat2, pch = NA, xlab = "Groundwater Influence", ylab = "Standard Deviation of log(Production)", axes = F, xlim = c(min(x.scaled.gw), max(x.scaled.gw)))
# axis(1, at = x.scaled.gw, labels = x.axis.gw)
# axis(2)
# box(bty = "l")
# pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
# for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- exp(Mcmcdat[i,"log.sigma.alpha"] + Mcmcdat[i,"log.sigma.beta"]*x.gw) }
# pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
# pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
# pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
# polygon(c(x.gw, rev(x.gw)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty = 0)
# lines(pred_median ~ x.gw, lwd = 2)
# usr <- par("usr")
# par(usr = c(0,1,0,1))
# text(0.07, 0.95, labels = "(d)")
# par(usr = usr)
# 
# dev.off()






####################################################################################################
####################################################################################################
### FREQUENTIST
####################################################################################################
####################################################################################################



boxplot(dat$prodlen)



plot(log(prodlen_censor2 + 0.015) ~ doy, dat, col = dat$cols)
plot(log(prodwt_censor2 + 0.008) ~ doy, dat, col = dat$cols)


prodmod_len <- lm(log(prodlen_censor + 0.015) ~ zdoy + I(zdoy^2) + zgw + zdoy*zgw + I(zdoy^2)*zgw, data = dat2)

summary(prodmod_len)

plot(prodmod_len)
summary(prodmod_len)$adj.r.square
sqrt(mean(prodmod_len$residuals^2)) 

# predict from the fitted model
gwvec <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = 100)
gpmat <- matrix(NA, nrow = max(dat2$zdoy, na.rm = T)-min(dat2$zdoy, na.rm = T)+1, ncol = length(gwvec))
for (i in 1:length(gwvec)) { gpmat[,i] <- predict(prodmod2, newdata = list(gw = rep(gwvec[i], times = dim(gpmat)[1]), doy = c(min(dat2$doy, na.rm = T):max(dat2$doy, na.rm = T)))) }

# another way to visualize model output
doys <- seq(from = min(dat2$zdoy, na.rm = T), to = max(dat2$zdoy, na.rm = T), length.out = 100)
colsss <- viridis(2)
gwv <- c(min(dat2$zgw), max(dat2$zgw))


jpeg("Production by Time/ProductionByTime_SeasonalProduction_Length_dummyfig.jpg", units = "in", width = 7.5, height = 7, res = 1500)
par(mfrow = c(2,2), mar = c(3,4.5,0.5,0.5), mgp = c(2.5,1,0))

plot(log(prodlen_censor + 0.015) ~ zdoy, dat2, pch = 16, col = alpha(dat2$cols, 0.7), xlab = "", ylab = "Length-Based Production (log scale)", axes = F)
axis(1, at = c(213,244,274,305), labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "o")
abline(h = log(0.015), lty = 2)
for (i in 1:length(gwv)) {
  preds <- predict(prodmod_len, newdata = list(zgw = rep(gwv[i], times = length(doys)), zdoy = doys), interval = "confidence")
  polygon(c(doys, rev(doys)), c(preds[,2], rev(preds[,3])), col = scales::alpha(colsss[i], 0.2), lty=0)
  lines(preds[,1] ~ doys, col = colsss[i], lwd = 3)
}

plot(seq(from = 0, to = 0.45, length.out = length(doys)) ~ doys, dat, pch = NA, xlab = "", ylab = "Predicted Production (mm per meter per day)", axes = F)
axis(1, at = c(213,244,274,305), labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "o")
for (i in 1:length(gwv)) {
  preds <- predict(prodmod_len, newdata = list(gw = rep(gwv[i], times = length(doys)), doy = doys), interval = "confidence")
  lines(exp(preds[,1]) ~ doys, lwd = 3, col = colsss[i])
}

plot(seq(from = 0, to = 25, length.out = length(doys)) ~ doys, dat, pch = NA, xlab = "", ylab = "Cumulatative Production (mm per meter)", axes = F)
axis(1, at = c(213,244,274,305), labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "o")
for (i in 1:length(gwv)) {
  preds <- predict(prodmod_len, newdata = list(gw = rep(gwv[i], times = length(doys)), doy = doys), interval = "confidence")
  lines(cumsum(exp(preds[,1])) ~ doys, lwd = 3, col = colsss[i])
}

legend("topleft", legend = round(gwv, digits = 2), fill = colsss, bty = "n", title = "Groundwater \nInfluence")

dev.off()


