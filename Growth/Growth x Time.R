#==========================================================================================================================================#
# PROJECT NAME: Groundwater structures fish growth and production across a riverscape
# CREATOR: Jeffrey R. Baldock, WY Coop Unit, University of Wyoming, jbaldock@uwyo.edu 
# PURPOSE: Estimate effects of groundwater on temporal trends in YOY growth, 
# NOTES: 
# - follows modelling approach detailed in manuscript and in Letcher et al. (2022) Journal of Animal Ecology
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

# Censor growth data based on small sample sizes: include only if catch and priorcatch >= 5
dat <- dat %>% mutate(rgr_len_censor = ifelse(catch >= 5 & priorcatch >= 5, rgr_len, NA))
sum(!is.na(dat$rgr_len_censor))

# dat$pvalue <- ifelse(dat$pvalue > 1, NA, dat$pvalue)

# drop missing data z-score covariates 
dat2 <- dat %>% filter(!is.na(rgr_len_censor)) %>% mutate(zgw = c(scale(gw, center = T, scale = T)),
                                                          zdoy = c(scale(doy, center = T, scale = T)),
                                                          zpriorwt = c(scale(priorwtmean, center = T, scale = T)),
                                                          zpriorlen = c(scale(priorlenmean, center = T, scale = T)),
                                                          zdenslog = c(scale(density_log, center = T, scale = T)),
                                                          ztemp = c(scale(meantempper, center = T, scale = T)))

# create table of covariate means and stdevs for back calculating/plotting
scaletbl <- data.frame(var = c("gw", "doy", "priorwt", "priorlen"),
                       means = c(mean(dat2$gw), mean(dat2$doy), mean(dat2$priorwtmean), mean(dat2$priorlenmean)),
                       sds = c(sd(dat2$gw), sd(dat2$doy), sd(dat2$priorwtmean), sd(dat2$priorlenmean)))
write_csv(scaletbl, "Growth by Time/GrowthByTime_z-score_table.csv")

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

StreamNames <- dat2 %>% group_by(stream) %>% summarize(strid = unique(strid), gw = mean(gw)) %>% arrange(gw) %>% ungroup() %>% mutate(strid2 = 1:13) %>% 
  add_row(stream = c("dum1","dum2"), gw = c(0,1)) %>% 
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>%
  filter(!stream %in% c("dum1","dum2"))

# Sort 
dat2 <- dat2 %>% arrange(yrid, strid, sectid)

##---------------------------------------------------------------------------------------##
## Fit JAGS Models
##---------------------------------------------------------------------------------------##

# Given our study design (sampling shortly after/during emergence, short sampling intervals and thus overlapping size distributions among sampling events),
# much of the (bootstrapped) variation in growth rates likely stems from variability in emergence timing rather than individual variation in growth rates.
# B/c of overlapping distributions, bootstrapped measures of uncertainty are very large. Therefore, accounting for bootstrapped measures of 
# uncertainty with a state-space model formulation likely results in overly conservative estimates of the effects of covariates on growth, which may limit
# our ability to draw inferences about covariate effects AND explore residual variation in growth among streams that is not already account for in the model.
# Thus, running a model on the means (i.e., no state-space structure) is a justifiable path forward. 

# Length-based growth by time model

# specify candidate models
Covs_list <- list(
  model.matrix(~ zdoy + I(zdoy^2) + zgw + zdoy*zgw + I(zdoy^2)*zgw, data = dat2, na.rm = FALSE),
  model.matrix(~ zdoy + I(zdoy^2) + zgw + zdoy*zgw + I(zdoy^2)*zgw + zpriorlen, data = dat2, na.rm = FALSE)
)

# fit candidate JAGS models in for loop
Fit_list <- list() # empty list to store model fits
Rhat_list <- list()
for (i in 1:length(Covs_list)) {
  Covs_Sigma <- model.matrix(~ zgw, data = dat2, na.rm = FALSE) # covariates on Sigma...doesn't work if out of for loop
  jags.data <- list("nObs" = dim(Covs_list[[i]])[1], "L" = dat2$lenmean, "pL" = dat2$priorlenmean, "elap" = dat2$elapdays,                                                                       
                    "strid" = dat2$strid, "yrid" = dat2$yrid, "sectid" = dat2$sectid, "sect_per_str" = unq_sect_str$strid,                                                                        
                    "nYears" = length(unique(dat2$yrid)), "nStreams" = length(unique(dat2$strid)), "nSections" = length(unique(dat2$sectid)),                                                
                    "Covs" = Covs_list[[i]][,-1], "nCovs" = dim(Covs_list[[i]])[2]-1, 
                    "CovsMu" = apply(Covs_list[[i]][,-1], MARGIN = 2, FUN = mean, na.rm = T), "CovsSD" = apply(Covs_list[[i]][,-1], MARGIN = 2, FUN = sd, na.rm = T),
                    "Covs_Sigma" = matrix(Covs_Sigma[,-1]), "nSigmaCovs" = dim(Covs_Sigma)[2]-1)
  params <- c("alpha.adj", "alpha.year.adj", "alpha.stream.adj", "alpha.sect.adj", "beta", "gr", "Lexp", "L",
              "sigma.gr", "log.sigma.alpha", "sigma.alpha", "log.sigma.beta", "sigma.sb", 
              "sigma.yr", "sigma.str", "sigma.sect", "sigma.a", "sigma.b", 
              "loglik", "res", "var_fit", "var_fix", "var_res", "var_yr", "var_sect", "var_str", "margR2", "condR2")
  m <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "JAGS Models/GrowthModel_Letcher.txt", 
                     n.chains = 3, n.thin = 20, n.burnin = 50000, n.iter = 150000, DIC = TRUE)
  Fit_list[[i]] <- m
  MCMCtrace(m, ind = TRUE, filename = paste("Growth by Time/Traceplots/Growth_Length_m", i, "_traceplots.pdf", sep = "")) # traceplots
  Rhat_list[[i]] <- m$BUGSoutput$summary[,8][m$BUGSoutput$summary[,8] > 1.01] # any R-hat values too high?  
  print(i)
}

# any R-hat values too high?  
Rhat_list


##---------------------------------------------------------------------------------------##
## MODEL SELECTION - compare LOO and set top model
##---------------------------------------------------------------------------------------##

loo_list <- list()
for (i in 1:length(Fit_list)) {
  m <- Fit_list[[i]]
  loglik <- m$BUGSoutput$sims.list$loglik
  reff <- relative_eff(exp(loglik), chain_id = c(rep(1,5000),rep(2,5000),rep(3,5000)))
  loo_list[[i]] <- loo(loglik, r_eff = reff)
  print(i)
}
lc <- loo_compare(loo_list)
print(lc, simplify = FALSE, digits = 2)
plot(loo_list[[1]])
write.csv(data.frame(lc), "/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/Growth by Time/GrowthbyTime_LOO_output.csv")

# set the top model based on LOO
# top_mod <- Fit_list[[1]]


##---------------------------------------------------------------------------------------##
## Save MCMC samples from top model (and second best model, top + cubic temp)
##---------------------------------------------------------------------------------------##

# DOY and priorlength are strongly correlated (0.806) and we feel that growth allometry is thus accounted for by including DOY, which is the primary variable of interest
# Further, marginal effects plots of how growth changes through time are uninterpretable given estimated growth rates, likely stemming from high correlation between DOY 
# and priorlength. Finally, projected body size trajectories and estimates of growth capacity from each model are virtually identical, indicating that the two models are 
# effectively equivalent, but multicollinearity within the second model limits very important inference

# set model
top_mod <- Fit_list[[1]]

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
write_csv(as.data.frame(Mcmcdat), "Growth by Time/GrowthByTime_Length_mcmcsamps.csv")
Mcmcdat1 <- read.csv("Growth by Time/GrowthByTime_Length_mcmcsamps.csv")
write.csv(as.data.frame(param.summary), "Growth by Time/GrowthByTime_Length_ParameterSummary.csv", row.names = T)
param.summary_len <- read.csv("Growth by Time/GrowthByTime_Length_ParameterSummary.csv", row.names = 1)


# set top model
top_mod <- Fit_list[[2]]

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
write_csv(as.data.frame(Mcmcdat), "Growth by Time/GrowthByTime_LengthB_mcmcsamps.csv")
Mcmcdat2 <- read.csv("Growth by Time/GrowthByTime_LengthB_mcmcsamps.csv")

# save models as separate objects
top_mod1 <- Fit_list[[1]]
top_mod2 <- Fit_list[[2]]



# compare growth rates estimate from each model
jpeg("Growth by Time/Figures/GrowthByTime_EstGrowthModelComparison.jpg", units = "in", width = 8, height = 8, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), mfrow = c(2,2))

# temporal model with and without prior length
plot(apply(top_mod1$BUGSoutput$sims.list$gr, 2, mean) ~ apply(top_mod2$BUGSoutput$sims.list$gr, 2, mean), xlim = c(-0.1, 0.9), ylim = c(-0.1, 0.9),
     xlab = "Estimated Growth Rate - Temporal + Prior Length", ylab = "Estimated Growth Rate - Temporal only")
abline(a = 0, b = 1, lwd = 2, col = "red")
legend("topleft", legend = paste("Pearson's r = ", round(cor(apply(top_mod1$BUGSoutput$sims.list$gr, 2, mean), apply(top_mod2$BUGSoutput$sims.list$gr, 2, mean), method = "pearson"), digits = 3)), bty = "n")

plot.new() # skip panel

# temporal model without prior length and covariate model
plot(apply(top_mod1$BUGSoutput$sims.list$gr, 2, mean) ~ apply(top_mod$BUGSoutput$sims.list$gr, 2, mean), xlim = c(-0.1, 0.9), ylim = c(-0.1, 0.9),
     ylab = "Estimated Growth Rate - Temporal only", xlab = "Estimated Growth Rate - Covariate Model")
abline(a = 0, b = 1, lwd = 2, col = "red")
legend("topleft", legend = paste("Pearson's r = ", round(cor(apply(top_mod1$BUGSoutput$sims.list$gr, 2, mean), apply(top_mod$BUGSoutput$sims.list$gr, 2, mean), method = "pearson"), digits = 3)), bty = "n")

# temporal model with prior length and covariate model
plot(apply(top_mod2$BUGSoutput$sims.list$gr, 2, mean) ~ apply(top_mod$BUGSoutput$sims.list$gr, 2, mean), xlim = c(-0.1, 0.9), ylim = c(-0.1, 0.9),
     ylab = "Estimated Growth Rate - Temporal + Prior Length", xlab = "Estimated Growth Rate - Covariate Model")
abline(a = 0, b = 1, lwd = 2, col = "red")
legend("topleft", legend = paste("Pearson's r = ", round(cor(apply(top_mod2$BUGSoutput$sims.list$gr, 2, mean), apply(top_mod$BUGSoutput$sims.list$gr, 2, mean), method = "pearson"), digits = 3)), bty = "n")

dev.off()


##---------------------------------------------------------------------------------------##
## MODEL DIAGNOSTICS
##---------------------------------------------------------------------------------------##
Covs <- model.matrix(~ zdoy + I(zdoy^2) + zgw + zdoy*zgw + I(zdoy^2)*zgw, data = dat2, na.rm = FALSE)

# subset expected and observed MCMC samples
ppdat_exp <- as.matrix(Mcmcdat1[,startsWith(names(Mcmcdat1), "Lexp.")])
ppdat_obs <- as.matrix(Mcmcdat1[,startsWith(names(Mcmcdat1), "L.")])

# Bayesian p-value
sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2])

# Posterior Predictive Check: plot median posterior expected length ~ observed length with Bayesian p-value
jpeg("Growth by Time/Figures/GrowthByTime_PPCheck.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(x = seq(from = 25, to = 75, length.out = 100), y = seq(from = 25, to = 75, length.out = 100), pch = NA, xlab = "Observed Length (mm)", ylab = "Median Posterior Expected Length (mm)")
for (i in 1:dim(Covs)[1]) { points(median(Mcmcdat1[,paste("Lexp.",i,".", sep = "")]) ~ dat2$lenmean[i])}
legend("topleft", bty = "n", legend = paste("Bayesian p-value = ", round(sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2]), digits = 3), sep = ""))
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()

# histogram of residuals
jpeg("Growth by Time/Figures/GrowthByTime_ResidHist.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist((ppdat_obs - ppdat_exp), main = "", xlab = "Observed - Expected")
legend("topright", bty = "n", legend = paste("Median Resid. = ", round(median((ppdat_obs - ppdat_exp)), digits = 2), " mm", sep = ""))
abline(v = median(unlist(ppdat_obs - ppdat_exp)), col = "red", lwd = 2)
dev.off()

# RMSE
sqrt(mean((ppdat_obs - ppdat_exp)^2)) # point estimate
rmse <- vector("numeric", 15000L)
for (j in 1:dim(Mcmcdat1)[1]) {
  res <- vector("numeric", 258L)
  for (i in 1:dim(Covs_list[[1]])[1]) {
    res[i] <- (ppdat_obs[j,i] - ppdat_exp[j,i])^2
  }
  rmse[j] <- sqrt(mean(res))
  print(j)
}
jpeg("Growth by Time/Figures/GrowthByTime_RMSE.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist(rmse, xlab = "Root Mean Square Error", ylab = "Posterior Frequency", main = "")
legend("topright", bty = "n", legend = paste("Median RMSE = ", round(median(rmse), digits = 2), " mm", sep = ""))
abline(v = median(rmse), col = "red", lwd = 2)
dev.off()


# Calculate Bayesian R^2 per Gelman et al. (2018), Nakagawa and Shielzeth (2013)
# Note: this does not work (as is) for state-space model formulation...conditional and marginal both approach 1
# Note: this is "conditional" R^2...that is, variance explained by both fixed and random effects combined. 
var_fix <- Mcmcdat1[,"var_fix"]
var_fit <- Mcmcdat1[,"var_fit"]
var_str <- Mcmcdat1[,"var_str"]
var_yr <- Mcmcdat1[,"var_yr"]
var_sect <- Mcmcdat1[,"var_sect"]
var_res <- Mcmcdat1[,"var_res"]

# model derived marginal and condition R2
hist(Mcmcdat1[,"margR2"])
hist(Mcmcdat1[,"condR2"])

# conditional R2 calculated in two ways
hist((var_fix + var_yr + var_str + var_sect) / (var_fix + var_res + var_yr + var_str + var_sect))
hist(var_fit / (var_fit + var_res), add = T)
median((var_fix + var_yr + var_str + var_sect) / (var_fix + var_res + var_yr + var_str + var_sect))
median(var_fit / (var_fit + var_res))

# random components each
hist((var_str) / (var_fix + var_str + var_yr + var_sect + var_res))
hist((var_sect) / (var_fix + var_str + var_yr + var_sect + var_res))
hist((var_yr) / (var_fix + var_str + var_yr + var_sect + var_res))
### Very broad distribution of variation explained by year component likely driven by few (2) factors associated with that hierarchical component.
### Multiple references (Gelman and Hill???) note that hierarchical variance may be poorly estimated (overly broad) when number of factors is small. 
### This results in distributions of marginal and conditional R^2 that are also very broad. Thus, should rely on modes for inference, as the mode 
### (more so than the median) reduces the influence of overly broad year component. Modes align well with preliminary analysis using frequentist
### approach (lmer), which yields marginal R2 = 0.401 and conditional R2 = 0.495
# hist((var_fix) / (var_fix + var_str + var_sect + var_res)) # demonstrate how dropping the year component results in more normal distribution
# hist((var_fix) / (var_fix + var_res)) 
# 
# # function to calculate mode
# getmode <- function(v) {
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }

# Marginal/conditional plot
jpeg("Growth by Time/Figures/GrowthByTime_MargCondR2.jpg", units = "in", width = 5.5, height = 5, res = 1500)
par(mar = c(4,4,0.5,0.5), mgp = c(2.5,1,0))
hist(Mcmcdat1[,"margR2"], col = alpha("darkorange", 0.5), xlim = c(0,1), ylim = c(0,8000), breaks = seq(from = 0, to = 1, by = 0.02), main = NA, xlab = expression(Bayesian ~ R^2), ylab = "Posterior Frequency")
hist(Mcmcdat1[,"condR2"], col = alpha("forestgreen", 0.5), xlim = c(0,1), add = TRUE, breaks = seq(from = 0, to = 1, by = 0.02))
# abline(v = median(Mcmcdat1[,"margR2"]), col = "darkorange", lwd = 3)
# abline(v = median(Mcmcdat1[,"condR2"]), col = "forestgreen", lwd = 3)
legend("topleft", fill = c(alpha("darkorange", 0.5), alpha("forestgreen", 0.5)), bty = "n", cex = 0.9,
       legend = c(paste("Marginal = ", round(median(Mcmcdat1[,"margR2"]), 3), " (", round(hdi(Mcmcdat1[,"margR2"], credMass = 0.95)[1], 3), ", ", round(hdi(Mcmcdat1[,"margR2"], credMass = 0.95)[2], 3), ")", sep = ""), 
                  paste("Conditional = ", round(median(Mcmcdat1[,"condR2"]), 3), " (", round(hdi(Mcmcdat1[,"condR2"], credMass = 0.95)[1], 3), ", ", round(hdi(Mcmcdat1[,"condR2"], credMass = 0.95)[2], 3), ")", sep = "")))

dev.off()


##### Combined Figure
jpeg("Growth by Time/Figures/GrowthByTime_ModelDiagnosticsCombined.jpg", units = "in", width = 7, height = 7, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), mfrow = c(2,2))

# PP Check
plot(x = seq(from = 25, to = 80, length.out = 100), y = seq(from = 25, to = 80, length.out = 100), pch = NA, xlab = "Observed Length (mm)", ylab = "Median Posterior Expected Length (mm)", bty = "l")
for (i in 1:dim(Covs)[1]) { points(median(Mcmcdat1[,paste("Lexp.",i,".", sep = "")]) ~ dat2$lenmean[i])}
legend("top", bty = "n", legend = paste("Bayesian p-value = ", round(sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2]), digits = 3), sep = ""))
abline(a = 0, b = 1, col = "red", lwd = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(a)")
par(usr = usr)

# RMSE
hist(rmse, xlab = "Root Mean Square Error", ylab = "Posterior Frequency", main = "")
box(bty = "l")
legend("topright", bty = "n", legend = paste("Median = ", round(median(rmse), digits = 2), " mm", sep = ""))
abline(v = median(rmse), col = "red", lwd = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(b)")
par(usr = usr)

# Marg/Cond R2
hist(Mcmcdat1[,"margR2"], col = alpha("darkorange", 0.5), xlim = c(0,1), ylim = c(0,9000), breaks = seq(from = 0, to = 1, by = 0.02), main = NA, xlab = expression(Bayesian ~ R^2), ylab = "Posterior Frequency")
box(bty = "l")
hist(Mcmcdat1[,"condR2"], col = alpha("forestgreen", 0.5), xlim = c(0,1), add = TRUE, breaks = seq(from = 0, to = 1, by = 0.02))
legend("topright", fill = c(alpha("darkorange", 0.5), alpha("forestgreen", 0.5)), bty = "n", cex = 0.9,
       legend = c(paste("Marginal = ", round(median(Mcmcdat1[,"margR2"]), 3), " (", round(hdi(Mcmcdat1[,"margR2"], credMass = 0.95)[1], 3), ", ", round(hdi(Mcmcdat1[,"margR2"], credMass = 0.95)[2], 3), ")", sep = ""), 
                  paste("Conditional = ", round(median(Mcmcdat1[,"condR2"]), 3), " (", round(hdi(Mcmcdat1[,"condR2"], credMass = 0.95)[1], 3), ", ", round(hdi(Mcmcdat1[,"condR2"], credMass = 0.95)[2], 3), ")", sep = "")))
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(c)")
par(usr = usr)

# residuals by Groundwater
plot(apply(as.matrix(Mcmcdat1[,startsWith(names(Mcmcdat1), "res")]), 2, mean) ~ dat2$gw, xlab = "Groundwater Index", ylab = "Residuals", bty = "l")
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(d)")
par(usr = usr)

dev.off()



##---------------------------------------------------------------------------------------##
## PLOTTING
##---------------------------------------------------------------------------------------##

mod.gg <- ggs(as.mcmc(top_mod1))
Covs <- model.matrix(~ zdoy + I(zdoy^2) + zgw + zdoy*zgw + I(zdoy^2)*zgw, data = dat2, na.rm = FALSE)

# Parameter (beta) estimates - dot plot
jpeg("Growth by Time/Figures/GrowthByTime_ParameterDotPlot.jpg", units = "in", width = 5, height = 4, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^beta", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = colnames(Covs)[-1])
dev.off()

# Parameter (log.sigma.beta: covariates on uncertainty) estimates - dot plot
jpeg("Growth by Time/Figures/GrowthByTime_ParameterDotPlot_SigmaBetas.jpg", units = "in", width = 5, height = 4, res = 1500)
ggs_caterpillar(D = mod.gg, family = "log.sigma.beta", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = "zgw")
dev.off()

# alpha.stream - stream offsets to the intercept
StreamNames <- StreamNames %>% arrange(strid)
jpeg("Growth by Time/Figures/GrowthByTime_ParameterDotPlot_alpha.stream.jpg", units = "in", width = 5, height = 4, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^alpha.stream.adj", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) +
  theme_bw() + ylab("") + xlab("Posterior estimate") +
  scale_y_discrete(labels = StreamNames$stream)
dev.off()

# alpha.section - section offsets to the intercept
jpeg("Growth by Time/Figures/GrowthByTime_ParameterDotPlot_alpha.section.jpg", units = "in", width = 5, height = 9, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^alpha.sect.adj", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = unq_sect_str$sectid)
dev.off()


#-------------------------------------------#
# Random Intercepts: Stream x Year
#-------------------------------------------#

Mcmcdat <- Mcmcdat1
StreamNames <- StreamNames %>% arrange(strid2)


jpeg("Growth by Time/Figures/GrowthByTime_PosteriorIntercepts.jpg", units = "in", width = 3.75, height = 6, res = 1500)
par(mar = c(4,3,1.5,0.5), mgp = c(2.5,1,0))
# set up plot area
x <- seq(from = 0.3, to = 0.6, length.out = 100)
y <- seq(from = 1, to = 13.7, length.out = 100)
plot(y ~ x, type = "n", xlab = expression(paste("Mean Growth Rate (% length day"^"-1", ")", sep = "")), ylab = "", axes = F)
axis(1)
# iterate over intercepts
ctr <- 0
for (j in nStreams:1) {
  par(xpd = NA)
  text(y = j+0.5, x = 0.25, labels = StreamNames$stream[j], pos = 4)
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
legend(x = 0.35, y = 15.2, legend = c(2021, 2022), lty = 1:2, lwd = 1.5, bty = "n", horiz = TRUE, cex = 0.9)
dev.off()


jpeg("Growth by Time/Figures/GrowthByTime_PosteriorIntercepts_withYears.jpg", units = "in", width = 3.75, height = 6, res = 1500)
par(mar = c(4,3,1.5,0.5), mgp = c(2.5,1,0))
# set up plot area
x <- seq(from = 0.3, to = 0.6, length.out = 100)
y <- seq(from = 1, to = 13.7, length.out = 100)
plot(y ~ x, type = "n", xlab = expression(paste("Mean Growth Rate (% length day"^"-1", ")", sep = "")), ylab = "", axes = F)
axis(1)
# iterate over intercepts
ctr <- 0
for (j in nStreams:1) {
  par(xpd = NA)
  text(y = j+0.5, x = 0.25, labels = StreamNames$stream[j], pos = 4)
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
legend(x = 0.35, y = 15.2, legend = c(2021, 2022), lty = 1:2, lwd = 1.5, bty = "n", horiz = TRUE, cex = 0.9)
dev.off()


#-------------------------------------------#
# Marginal Effects Plots
#-------------------------------------------#

# set up gradient legend
gradpal <- tibble(gw = seq(from = 0, to = 1, length.out = 100)) %>%
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>% 
  filter(gw >= min(dat2$gw) & gw <= max(dat2$gw))
# gradientLegend(valRange = c(min(gradpal$gw), max(gradpal$gw)), color = gradpal$cols, length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.05, 0.05, 0.12, 0.4), dec = 2)
# number of values for prediction
nvalues <- 100
# generate sequences of predictor variables and scales axes
x.gw <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = nvalues)
x.doy <- seq(from = min(dat2$zdoy), to = max(dat2$zdoy), length.out = nvalues)
x.plen <- seq(from = min(dat2$zpriorlen), to = max(dat2$zpriorlen), length.out = nvalues)
# manually set up and transform axes to original scale of the data
x.axis.gw <- seq(from = 0, to = 0.8, by = 0.1)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
x.axis.doy <- c(213,244,274,305)
x.scaled.doy <- (x.axis.doy - scaletbl$means[2]) / scaletbl$sds[2]

# Model without prior length
jpeg("Growth by Time/Figures/GrowthByTime_TimeGWEffect_wPoints.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# growth ~ temp + density (temp focal)
plot(seq(from = 0, to = 0.7, length.out = length(dat2$zdoy)) ~ zdoy, dat2, pch = NA, xlab = "", ylab = expression(paste("Growth Rate (mm day"^"-1", ")", sep = "")), axes = F)
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "l")
# interaction: doy * high gw
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat1), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat1[i,"alpha.adj"] + Mcmcdat1[i,"beta.1."]*x.doy + Mcmcdat1[i,"beta.2."]*(x.doy^2) + Mcmcdat1[i,"beta.3."]*max(x.gw) + Mcmcdat1[i,"beta.4."]*max(x.gw)*x.doy + Mcmcdat1[i,"beta.5."]*max(x.gw)*(x.doy^2) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == max(coldf$gw), "cols"], 0.2), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == max(coldf$gw), "cols"]), lwd = 2)
# interaction: doy * low gw
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat1), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat1[i,"alpha.adj"] + Mcmcdat1[i,"beta.1."]*x.doy + Mcmcdat1[i,"beta.2."]*(x.doy^2) + Mcmcdat1[i,"beta.3."]*min(x.gw) + Mcmcdat1[i,"beta.4."]*min(x.gw)*x.doy + Mcmcdat1[i,"beta.5."]*min(x.gw)*(x.doy^2) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == min(coldf$gw), "cols"], 0.2), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == min(coldf$gw), "cols"]), lwd = 2)
# points
for (i in 1:dim(dat2)[1]) { points(median(Mcmcdat1[,paste("gr.",i,".", sep = "")]) ~ dat2$zdoy[i], pch = 16, col = alpha(dat2$cols, 0.5)[i]) }
# legend
text(x = 1.2, y = 4.0, labels = "Groundwater\nInfluence", cex = 0.9)
gradientLegend(valRange = c(min(gradpal$gw), max(gradpal$gw)), color = gradpal$cols, length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.82, 0.65, 0.87, 1), dec = 2)
dev.off()


# Model with prior length
jpeg("Growth by Time/Figures/GrowthByTime_TimeGWEffect_wPoints_WithPriorLen.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# growth ~ temp + density (temp focal)
plot(seq(from = -0.1, to = 1, length.out = length(dat2$zdoy)) ~ zdoy, dat2, pch = NA, xlab = "", ylab = expression(paste("Growth Rate (mm day"^"-1", ")", sep = "")), axes = F)
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "l")
# interaction: doy * high gw
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat2), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat2[i,"alpha.adj"] + Mcmcdat2[i,"beta.1."]*x.doy + Mcmcdat2[i,"beta.2."]*(x.doy^2) + Mcmcdat2[i,"beta.3."]*max(x.gw) + Mcmcdat2[i,"beta.5."]*max(x.gw)*x.doy + Mcmcdat2[i,"beta.6."]*max(x.gw)*(x.doy^2) + Mcmcdat2[i,"beta.4."]*((lenvecB_med_mat[j,86] - scaletbl$means[4]) / scaletbl$sds[4]) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == max(coldf$gw), "cols"], 0.2), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == max(coldf$gw), "cols"]), lwd = 2)
# interaction: doy * low gw
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat2), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat2[i,"alpha.adj"] + Mcmcdat2[i,"beta.1."]*x.doy + Mcmcdat2[i,"beta.2."]*(x.doy^2) + Mcmcdat2[i,"beta.3."]*min(x.gw) + Mcmcdat2[i,"beta.5."]*min(x.gw)*x.doy + Mcmcdat2[i,"beta.6."]*min(x.gw)*(x.doy^2) + Mcmcdat2[i,"beta.4."]*((lenvecB_med_mat[j,5] - scaletbl$means[4]) / scaletbl$sds[4]) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == min(coldf$gw), "cols"], 0.2), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == min(coldf$gw), "cols"]), lwd = 2)
# points
for (i in 1:dim(dat2)[1]) { points(median(Mcmcdat2[,paste("gr.",i,".", sep = "")]) ~ dat2$zdoy[i], pch = 16, col = alpha(dat2$cols, 0.5)[i]) }
# legend
text(x = 1.2, y = 4.0, labels = "Groundwater\nInfluence", cex = 0.9)
gradientLegend(valRange = c(min(gradpal$gw), max(gradpal$gw)), color = gradpal$cols, length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.82, 0.65, 0.87, 1), dec = 2)
dev.off()


####### For Presentation
png("Growth by Time/Figures/Presentation Figs/GrowthByTime_TimeGWEffect_Blank.jpg", units = "in", width = 4, height = 3.75, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), bg = NA)
plot(seq(from = 0, to = 0.7, length.out = length(dat2$zdoy)) ~ zdoy, dat2, pch = NA, xlab = "", ylab = expression(paste("Growth Rate (mm day"^"-1", ")", sep = "")), axes = F)
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "l")
dev.off()


png("Growth by Time/Figures/Presentation Figs/GrowthByTime_TimeGWEffect_HighGW.jpg", units = "in", width = 4, height = 3.75, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), bg = NA)
plot(seq(from = 0, to = 0.7, length.out = length(dat2$zdoy)) ~ zdoy, dat2, pch = NA, xlab = "", ylab = "", axes = F)
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat1), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat1[i,"alpha.adj"] + Mcmcdat1[i,"beta.1."]*x.doy + Mcmcdat1[i,"beta.2."]*(x.doy^2) + Mcmcdat1[i,"beta.3."]*max(x.gw) + Mcmcdat1[i,"beta.4."]*max(x.gw)*x.doy + Mcmcdat1[i,"beta.5."]*max(x.gw)*(x.doy^2) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == max(coldf$gw), "cols"], 0.2), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == max(coldf$gw), "cols"]), lwd = 2)
dev.off()


png("Growth by Time/Figures/Presentation Figs/GrowthByTime_TimeGWEffect_LowGW.jpg", units = "in", width = 4, height = 3.75, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), bg = NA)
plot(seq(from = 0, to = 0.7, length.out = length(dat2$zdoy)) ~ zdoy, dat2, pch = NA, xlab = "", ylab = "", axes = F)
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat1), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat1[i,"alpha.adj"] + Mcmcdat1[i,"beta.1."]*x.doy + Mcmcdat1[i,"beta.2."]*(x.doy^2) + Mcmcdat1[i,"beta.3."]*min(x.gw) + Mcmcdat1[i,"beta.4."]*min(x.gw)*x.doy + Mcmcdat1[i,"beta.5."]*min(x.gw)*(x.doy^2) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == min(coldf$gw), "cols"], 0.2), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == min(coldf$gw), "cols"]), lwd = 2)
dev.off()


#-------------------------------------------#
# Growth capacity ~ time ~ groundwater
#-------------------------------------------#

# number of groundwater values
nvalues <- 100

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

# empty matrices to store projected body size
lenvec_med_mat <- matrix(data = NA, nrow = length(x.doy)+1, ncol = nvalues)
lenvecB_med_mat <- matrix(data = NA, nrow = length(x.doy)+1, ncol = nvalues)
gr <- matrix(data = NA, nrow = length(x.doy), ncol = nvalues)

# starting body size in mm -- approx. average size in mid-August
lenvec_med_mat[1,] <- 35 
lenvecB_med_mat[1,] <- 35 

# Iteratively project body size through time
for (i in 1:nvalues) {
  for (j in 1:length(x.doy)) {
    # lenvec_med_mat[j+1,i] <- lenvec_med_mat[j,i] + (median(Mcmcdat1[,"alpha.adj"]) + median(Mcmcdat1[,"beta.1."])*x.doy[j] + median(Mcmcdat1[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat1[,"beta.3."])*x.gw[i] + median(Mcmcdat1[,"beta.4."])*x.gw[i]*x.doy[j] + median(Mcmcdat1[,"beta.5."])*x.gw[i]*(x.doy[j]^2))
    gr[j,i] <- (median(Mcmcdat2[,"alpha.adj"]) + median(Mcmcdat2[,"beta.1."])*x.doy[j] + median(Mcmcdat2[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat2[,"beta.3."])*x.gw[i] + median(Mcmcdat2[,"beta.4."])*((lenvecB_med_mat[j,i] - scaletbl$means[4]) / scaletbl$sds[4]) + median(Mcmcdat2[,"beta.5."])*x.gw[i]*x.doy[j] + median(Mcmcdat2[,"beta.6."])*x.gw[i]*(x.doy[j]^2))
    lenvecB_med_mat[j+1,i] <- lenvecB_med_mat[j,i] + gr[j,i]
    }
  print(i)
}


# projected growth capcity without prior length
jpeg("Growth by Time/Figures/GrowthByTime_GrowthCapacity_Length.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
colss <- rev(hcl.colors(n = nvalues, "Viridis"))
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(seq(from = lenvec_med_mat[1,1], to = max(lenvec_med_mat), length.out = length(doys)) ~ x.doy, dat2, pch = NA, axes = F, xlab = "", ylab = "Length (mm)")
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "l")
for (i in 1:dim(lenvec_med_mat)[2]) { lines(lenvec_med_mat[-79,i] ~ x.doy, col = colss[i], lwd = 1.5) }
# add lines for min/max observed gw influence
dum_vec <- matrix(data = NA, nrow = length(x.doy)+1, ncol = 2)
dum_vec[1,] <- 35 
for (j in 1:length(x.doy)) {
  dum_vec[j+1,1] <- dum_vec[j,1] + (median(Mcmcdat1[,"alpha.adj"]) + median(Mcmcdat1[,"beta.1."])*x.doy[j] + median(Mcmcdat1[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat1[,"beta.3."])*min(dat2$zgw) + median(Mcmcdat1[,"beta.4."])*min(dat2$zgw)*x.doy[j] + median(Mcmcdat1[,"beta.5."])*min(dat2$zgw)*(x.doy[j]^2))
  dum_vec[j+1,2] <- dum_vec[j,2] + (median(Mcmcdat1[,"alpha.adj"]) + median(Mcmcdat1[,"beta.1."])*x.doy[j] + median(Mcmcdat1[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat1[,"beta.3."])*max(dat2$zgw) + median(Mcmcdat1[,"beta.4."])*max(dat2$zgw)*x.doy[j] + median(Mcmcdat1[,"beta.5."])*max(dat2$zgw)*(x.doy[j]^2))
}
for (i in 1:2) { lines(dum_vec[-79,i] ~ x.doy, lty = 1, lwd = 0.75, col = "white") }
# legend
text(x = 1.25, y = 50, labels = "Groundwater \nindex", cex = 0.8)
gradientLegend(valRange = c(0,1), color = rev(viridis(100)), length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.82, 0.05, 0.87, 0.4), dec = 2)
par(usr = c(0,1,0,1))
lines(x = c(0.822, 0.868), y = c((min(dat2$gw)*(0.4-0.05))+0.05, (min(dat2$gw)*(0.4-0.05))+0.05), lty = 1, lwd = 0.75, col = "white")
lines(x = c(0.822, 0.868), y = c((max(dat2$gw)*(0.4-0.05))+0.05, (max(dat2$gw)*(0.4-0.05))+0.05), lty = 1, lwd = 0.75, col = "white")
# # inset figure showing relative growth advantage
# par(fig = c(0.11,0.7,0.4,1), mar = c(4,4,1,1), mgp = c(1.2,0.4,0), new = T)
# plot(((lenvec_med_mat[79,] - lenvec_med_mat[79,1]) / lenvec_med_mat[79,1])*100 ~ x.gw, type = "l", lwd = 2, ylab = "Growth Advantage (%)", xlab = "Groundwater Influence", axes = F, cex.lab = 0.8, ylim = c(0,86))
# axis(1, cex.axis = 0.8, tck = -0.04, at = x.scaled.gw, labels = x.axis.gw)
# axis(2, cex.axis = 0.8, tck = -0.04)
# box(bty = "l")
dev.off()


# Projected growth capacity with prior length
jpeg("Growth by Time/Figures/GrowthByTime_GrowthCapacity_Length_WithPriorLen.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(seq(from = lenvecB_med_mat[1,1], to = max(lenvecB_med_mat), length.out = length(doys)) ~ x.doy, dat2, pch = NA, axes = F, xlab = "", ylab = "Length (mm)")
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "l")
for (i in 1:dim(lenvecB_med_mat)[2]) { lines(lenvecB_med_mat[-79,i] ~ x.doy, col = colss[i], lwd = 1.5) }
# add lines for min/max observed gw influence
dum_vec <- matrix(data = NA, nrow = length(x.doy)+1, ncol = 2)
dum_vec[1,] <- 35 
for (j in 1:length(x.doy)) {
  dum_vec[j+1,1] <- dum_vec[j,1] + (median(Mcmcdat2[,"alpha.adj"]) + median(Mcmcdat2[,"beta.1."])*x.doy[j] + median(Mcmcdat2[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat2[,"beta.3."])*min(dat2$zgw) + median(Mcmcdat2[,"beta.4."])*((dum_vec[j,1] - scaletbl$means[4]) / scaletbl$sds[4]) + median(Mcmcdat2[,"beta.5."])*min(dat2$zgw)*x.doy[j] + median(Mcmcdat2[,"beta.6."])*min(dat2$zgw)*(x.doy[j]^2))
  dum_vec[j+1,2] <- dum_vec[j,2] + (median(Mcmcdat2[,"alpha.adj"]) + median(Mcmcdat2[,"beta.1."])*x.doy[j] + median(Mcmcdat2[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat2[,"beta.3."])*max(dat2$zgw) + median(Mcmcdat2[,"beta.4."])*((dum_vec[j,2] - scaletbl$means[4]) / scaletbl$sds[4]) + median(Mcmcdat2[,"beta.5."])*max(dat2$zgw)*x.doy[j] + median(Mcmcdat2[,"beta.6."])*max(dat2$zgw)*(x.doy[j]^2))
}
for (i in 1:2) { lines(dum_vec[-79,i] ~ x.doy, lty = 1, lwd = 0.75, col = "white") }
# legend
text(x = 1.25, y = 45, labels = "Groundwater \nindex", cex = 0.8)
gradientLegend(valRange = c(0,1), color = rev(viridis(100)), length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.82, 0.05, 0.87, 0.4), dec = 2)
par(usr = c(0,1,0,1))
lines(x = c(0.822, 0.868), y = c((min(dat2$gw)*(0.4-0.05))+0.05, (min(dat2$gw)*(0.4-0.05))+0.05), lty = 1, lwd = 0.75, col = "white")
lines(x = c(0.822, 0.868), y = c((max(dat2$gw)*(0.4-0.05))+0.05, (max(dat2$gw)*(0.4-0.05))+0.05), lty = 1, lwd = 0.75, col = "white")
# # inset figure showing relative growth advantage
# par(fig = c(0.11,0.7,0.4,1), mar = c(4,4,1,1), mgp = c(1.2,0.4,0), new = T)
# plot(((lenvec_med_mat[79,] - lenvec_med_mat[79,1]) / lenvec_med_mat[79,1])*100 ~ x.gw, type = "l", lwd = 2, ylab = "Growth Advantage (%)", xlab = "Groundwater Influence", axes = F, cex.lab = 0.8, ylim = c(0,86))
# axis(1, cex.axis = 0.8, tck = -0.04, at = x.scaled.gw, labels = x.axis.gw)
# axis(2, cex.axis = 0.8, tck = -0.04)
# box(bty = "l")
dev.off()


plot(seq(from = -0.1, to = 1, length.out = length(doys)) ~ x.doy, dat2, pch = NA, axes = F, xlab = "", ylab = "Growth")
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "l")
for (i in 1:dim(gr)[2]) { lines(gr[,i] ~ x.doy, col = colss[i], lwd = 1.5) }
for (i in 1:dim(dat2)[1]) { points(median(Mcmcdat2[,paste("gr.",i,".", sep = "")]) ~ dat2$zdoy[i], pch = 16, col = alpha(dat2$cols, 0.5)[i]) }


# project body size similar to above but with credible intervals
pred_dist_mat <- array(data = NA, dim = c(nrow(Mcmcdat), length(x.doy), nvalues))
pred_med_mat <- matrix(data = NA, nrow = length(x.doy), ncol = nvalues)
pred_low_mat <- matrix(data = NA, nrow = length(x.doy), ncol = nvalues)
pred_upp_mat <- matrix(data = NA, nrow = length(x.doy), ncol = nvalues)
lenvec_med_mat <- matrix(data = NA, nrow = length(x.doy)+1, ncol = nvalues)
lenvec_low_mat <- matrix(data = NA, nrow = length(x.doy)+1, ncol = nvalues)
lenvec_upp_mat <- matrix(data = NA, nrow = length(x.doy)+1, ncol = nvalues)
lenvec_med_mat[1,] <- 35 # starting body size in mm (25) of grams (0.10) -- Approx. size at emergence
lenvec_low_mat[1,] <- 35 # starting body size in mm (25) of grams (0.10) -- Approx. size at emergence
lenvec_upp_mat[1,] <- 35 # starting body size in mm (25) of grams (0.10) -- Approx. size at emergence
for (i in 1:nvalues) {
  for (j in 1:nrow(Mcmcdat)) { pred_dist_mat[j,,i] <- Mcmcdat[j,"alpha.adj"] + Mcmcdat[j,"beta.1."]*x.doy + Mcmcdat[j,"beta.2."]*(x.doy^2) + Mcmcdat[j,"beta.3."]*x.gw[i] + Mcmcdat[j,"beta.4."]*x.gw[i]*x.doy + Mcmcdat[j,"beta.5."]*x.gw[i]*(x.doy^2) }
  pred_med_mat[,i] <- apply(pred_dist_mat[,,i], MARGIN = 2, quantile, prob = 0.5)
  pred_low_mat[,i] <- apply(pred_dist_mat[,,i], MARGIN = 2, quantile, prob = 0.025)
  pred_upp_mat[,i] <- apply(pred_dist_mat[,,i], MARGIN = 2, quantile, prob = 0.975)
  for (j in 1:length(x.doy)) {
    lenvec_med_mat[j+1,i] <- lenvec_med_mat[j,i] + pred_med_mat[j,i]
    lenvec_low_mat[j+1,i] <- lenvec_low_mat[j,i] + pred_low_mat[j,i]
    lenvec_upp_mat[j+1,i] <- lenvec_upp_mat[j,i] + pred_upp_mat[j,i]
  }
  print(i)
}


# # projected growth capcity without prior length
# jpeg("Growth by Time/Figures/GrowthByTime_GrowthCapacity_Length.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
# colss <- hcl.colors(n = nvalues, "Viridis")
# par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# plot(seq(from = lenvec_med_mat[1,1], to = max(lenvec_med_mat), length.out = length(doys)) ~ x.doy, dat2, pch = NA, axes = F, xlab = "", ylab = "Length (mm)")
# axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
# axis(2)
# box(bty = "l")
# for (i in 1:dim(lenvec_med_mat)[2]) { lines(lenvec_med_mat[-79,i] ~ x.doy, col = colss[i], lwd = 1.5) }
# # add lines for min/max observed gw influence
# dum_vec <- matrix(data = NA, nrow = length(x.doy)+1, ncol = 2)
# dum_vec[1,] <- 35 
# for (j in 1:length(x.doy)) {
#   dum_vec[j+1,1] <- dum_vec[j,1] + (median(Mcmcdat1[,"alpha.adj"]) + median(Mcmcdat1[,"beta.1."])*x.doy[j] + median(Mcmcdat1[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat1[,"beta.3."])*min(dat2$zgw) + median(Mcmcdat1[,"beta.4."])*min(dat2$zgw)*x.doy[j] + median(Mcmcdat1[,"beta.5."])*min(dat2$zgw)*(x.doy[j]^2))
#   dum_vec[j+1,2] <- dum_vec[j,2] + (median(Mcmcdat1[,"alpha.adj"]) + median(Mcmcdat1[,"beta.1."])*x.doy[j] + median(Mcmcdat1[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat1[,"beta.3."])*max(dat2$zgw) + median(Mcmcdat1[,"beta.4."])*max(dat2$zgw)*x.doy[j] + median(Mcmcdat1[,"beta.5."])*max(dat2$zgw)*(x.doy[j]^2))
# }
# for (i in 1:2) { lines(dum_vec[-79,i] ~ x.doy, lty = 3, lwd = 1, col = "white") }
# # legend
# text(x = 1.25, y = 50, labels = "Groundwater \nInfluence", cex = 0.8)
# gradientLegend(valRange = c(0,1), color = viridis(100), length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.82, 0.05, 0.87, 0.4), dec = 2)
# dev.off()

# projected end of season body size by groundwater influence, with credible interval
jpeg("Growth by Time/Figures/GrowthByTime_GrowthCapacity_EndOfSeasonLengthByGroundwater.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(lenvec_med_mat[79,] ~ x.gw, pch = NA, ylab = "End of Season Length (mm)", xlab = "Groundwater Influence", ylim = c(55,75), bty = "l", axes = F)
polygon(x = c(x.gw, rev(x.gw)), y = c(lenvec_low_mat[79,], rev(lenvec_upp_mat[79,])), col = alpha("black", 0.2), lty = 0)
lines(lenvec_med_mat[79,] ~ x.gw, lwd = 2)
axis(1, cex.axis = 0.8, tck = -0.04, at = x.scaled.gw, labels = x.axis.gw)
axis(2, cex.axis = 0.8, tck = -0.04)
box(bty = "l")
dev.off()

####
# Write out/store projected length matrices
write_csv(as_tibble(lenvec_med_mat), "Growth by Time/GrowthByTime_ProjectedLength_50perc.csv")
write_csv(as_tibble(lenvec_low_mat), "Growth by Time/GrowthByTime_ProjectedLength_025perc.csv")
write_csv(as_tibble(lenvec_upp_mat), "Growth by Time/GrowthByTime_ProjectedLength_975perc.csv")

lenvec_med_mat <- as.matrix(read_csv("Growth by Time/GrowthByTime_ProjectedLength_50perc.csv"))
lenvec_low_mat <- as.matrix(read_csv("Growth by Time/GrowthByTime_ProjectedLength_025perc.csv"))
lenvec_upp_mat <- as.matrix(read_csv("Growth by Time/GrowthByTime_ProjectedLength_975perc.csv"))


#-------------------------------------------#
# Groundwater effect on sigma
#-------------------------------------------#

# number of values for prediction
nvalues <- 100
# generate sequences of predictor variables and scales axes
x.gw <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = nvalues)
x.axis.gw <- seq(from = 0, to = 0.8, by = 0.2)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]

# Model without prior length
jpeg("Growth by Time/Figures/GrowthByTime_GWEffectOnSigma.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# growth ~ temp + density (temp focal)
plot(seq(from = 0.12, to = 0.28, length.out = length(dat2$zdoy)) ~ zgw, dat2, pch = NA, xlab = "Groundwater index", ylab = "Standard Deviation of Growth Rate", axes = F)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat1), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- exp(Mcmcdat1[i,"log.sigma.alpha"] + Mcmcdat1[i,"log.sigma.beta"]*x.gw) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.gw, rev(x.gw)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty = 0)
lines(pred_median ~ x.gw, lwd = 2)
dev.off()


#-------------------------------------------#
# Combined Figure
#-------------------------------------------#

# Model without prior length
jpeg("Growth by Time/Figures/GrowthByTime_CombinedFigure.jpg", units = "in", width = 6.5, height = 6, res = 1500)
par(mar = c(3,4,2,1), mgp = c(2.5,1,0), mfrow = c(2,2))

# growth ~ time + groundwater
nvalues <- 100
x.gw <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = nvalues)
x.doy <- seq(from = min(dat2$zdoy), to = max(dat2$zdoy), length.out = nvalues)
x.axis.gw <- seq(from = 0, to = 0.8, by = 0.1)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
x.axis.doy <- c(213,244,274,305)
x.scaled.doy <- (x.axis.doy - scaletbl$means[2]) / scaletbl$sds[2]
# set up plot
plot(seq(from = 0, to = 0.7, length.out = length(dat2$zdoy)) ~ zdoy, dat2, pch = NA, xlab = "", ylab = expression(paste("Growth rate (mm day"^"-1", ")", sep = "")), axes = F)
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "l")
# interaction: doy * high gw
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat1), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat1[i,"alpha.adj"] + Mcmcdat1[i,"beta.1."]*x.doy + Mcmcdat1[i,"beta.2."]*(x.doy^2) + Mcmcdat1[i,"beta.3."]*max(x.gw) + Mcmcdat1[i,"beta.4."]*max(x.gw)*x.doy + Mcmcdat1[i,"beta.5."]*max(x.gw)*(x.doy^2) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == max(coldf$gw), "cols"], 0.2), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == max(coldf$gw), "cols"]), lwd = 2)
# interaction: doy * low gw
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat1), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat1[i,"alpha.adj"] + Mcmcdat1[i,"beta.1."]*x.doy + Mcmcdat1[i,"beta.2."]*(x.doy^2) + Mcmcdat1[i,"beta.3."]*min(x.gw) + Mcmcdat1[i,"beta.4."]*min(x.gw)*x.doy + Mcmcdat1[i,"beta.5."]*min(x.gw)*(x.doy^2) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha("white", 0.5), lty=0)
polygon(c(x.doy, rev(x.doy)), c(pred_upper, rev(pred_lower)), col = scales::alpha(coldf[coldf$gw == min(coldf$gw), "cols"], 0.27), lty=0)
lines(pred_median ~ x.doy, col = as.character(coldf[coldf$gw == min(coldf$gw), "cols"]), lwd = 2)
# points
for (i in 1:dim(dat2)[1]) { points(median(Mcmcdat1[,paste("gr.",i,".", sep = "")]) ~ dat2$zdoy[i], pch = 16, col = alpha(dat2$cols, 0.5)[i]) }
# legend
gradpal <- tibble(gw = seq(from = 0, to = 1, length.out = 100)) %>%
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>% 
  filter(gw >= min(dat2$gw) & gw <= max(dat2$gw))
# text(x = 1.2, y = 4.0, labels = "Groundwater\nindex", cex = 0.9)
# gradientLegend(valRange = c(min(gradpal$gw), max(gradpal$gw)), color = gradpal$cols, length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.82, 0.65, 0.87, 1), dec = 2)
# panel label
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(a)")
par(usr = usr)


# groundwater effect on sigma
nvalues <- 100
x.gw <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = nvalues)
x.axis.gw <- seq(from = 0.2, to = 0.8, by = 0.2)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
plot(seq(from = 0.1, to = 0.3, length.out = length(dat2$zdoy)) ~ zgw, dat2, pch = NA, xlab = "", ylab = "Standard deviation of length", axes = F)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat1), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- exp(Mcmcdat1[i,"log.sigma.alpha"] + Mcmcdat1[i,"log.sigma.beta"]*x.gw) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.gw, rev(x.gw)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty = 0)
lines(pred_median ~ x.gw, lwd = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(b)")
par(usr = usr)

par(mar = c(4,4,1,1))

# projected growth capcity without prior length
gradpal <- tibble(gw = seq(from = 0, to = 1, length.out = 100)) %>%
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>% 
  filter(gw >= min(dat2$gw) & gw <= max(dat2$gw))
nvalues <- 100
xmin <- (0 - scaletbl$means[1]) / scaletbl$sds[1]
xmax <- (1 - scaletbl$means[1]) / scaletbl$sds[1]
doys <- c(min(dat2$doy, na.rm = T):max(dat2$doy, na.rm = T))
x.gw <- seq(from = xmin, to = xmax, length.out = nvalues)
x.doy <- (doys - scaletbl$means[2]) / scaletbl$sds[2]
x.axis.gw <- seq(from = 0, to = 1, by = 0.2)
x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
x.axis.doy <- c(213,244,274,305)
x.scaled.doy <- (x.axis.doy - scaletbl$means[2]) / scaletbl$sds[2]
colss <- rev(hcl.colors(n = nvalues, "Viridis"))
plot(seq(from = lenvec_med_mat[1,1], to = max(lenvec_med_mat), length.out = length(doys)) ~ x.doy, dat2, pch = NA, axes = F, xlab = "Time", ylab = "Projected length (mm)")
axis(1, at = x.scaled.doy, labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "l")
for (i in 1:dim(lenvec_med_mat)[2]) { lines(lenvec_med_mat[-79,i] ~ x.doy, col = colss[i], lwd = 1.5) }
# add lines for min/max observed gw influence
dum_vec <- matrix(data = NA, nrow = length(x.doy)+1, ncol = 2)
dum_vec[1,] <- 35 
for (j in 1:length(x.doy)) {
  dum_vec[j+1,1] <- dum_vec[j,1] + (median(Mcmcdat1[,"alpha.adj"]) + median(Mcmcdat1[,"beta.1."])*x.doy[j] + median(Mcmcdat1[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat1[,"beta.3."])*min(dat2$zgw) + median(Mcmcdat1[,"beta.4."])*min(dat2$zgw)*x.doy[j] + median(Mcmcdat1[,"beta.5."])*min(dat2$zgw)*(x.doy[j]^2))
  dum_vec[j+1,2] <- dum_vec[j,2] + (median(Mcmcdat1[,"alpha.adj"]) + median(Mcmcdat1[,"beta.1."])*x.doy[j] + median(Mcmcdat1[,"beta.2."])*(x.doy[j]^2) + median(Mcmcdat1[,"beta.3."])*max(dat2$zgw) + median(Mcmcdat1[,"beta.4."])*max(dat2$zgw)*x.doy[j] + median(Mcmcdat1[,"beta.5."])*max(dat2$zgw)*(x.doy[j]^2))
}
for (i in 1:2) { lines(dum_vec[-79,i] ~ x.doy, lty = 1, lwd = 0.75, col = "white") }
# legend
text(x = 1.2, y = 54, labels = "Groundwater \nindex", cex = 0.8)
gradientLegend(valRange = c(0,1), color = rev(viridis(100)), length = 0.5, inside = TRUE, side = 2, n.seg = 3, pos = c(0.80, 0.05, 0.87, 0.5), dec = 2)
# panel label
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(c)")
lines(x = c(0.802, 0.868), y = c((min(dat2$gw)*(0.5-0.05))+0.05, (min(dat2$gw)*(0.5-0.05))+0.05), lty = 1, lwd = 0.75, col = "white")
lines(x = c(0.802, 0.868), y = c((max(dat2$gw)*(0.5-0.05))+0.05, (max(dat2$gw)*(0.5-0.05))+0.05), lty = 1, lwd = 0.75, col = "white")
par(usr = usr)


# projected end of season body size by groundwater influence, with credible interval
plot(lenvec_med_mat[79,] ~ x.gw, pch = NA, ylab = "Projected end of season length (mm)", xlab = "Groundwater index", ylim = c(55,75), bty = "l", axes = F)
polygon(x = c(x.gw, rev(x.gw)), y = c(lenvec_low_mat[79,], rev(lenvec_upp_mat[79,])), col = alpha("black", 0.2), lty = 0)
lines(lenvec_med_mat[79,] ~ x.gw, lwd = 2)
axis(1, at = x.scaled.gw, labels = x.axis.gw)
axis(2)
box(bty = "l")
# panel label
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(d)")
par(usr = usr)

dev.off()





##---------------------------------------------------------------------------------------##
## How does the modeled temporal trend in growth rate affect end of season body size across 
## a gradient of GW input? Ideally, we could just compare measured/observed end of season 
## body size across streams. However, this doesn't account for differences in the time of 
## sampling (which differs by up to 7 days), nor does it account for differences in spawn/
## emergence timing and thus size at first sampling. Weight at first sampling differs by 
## almost an order of magnitude (0.078-0.651), indicating considerable differences in 
## emergence timing. As a result, using empirical data, there is no relationship between 
## groundwater index and end of season body size. There is also an apparent trend towards 
## larger body size at first sampling in snowmelt streams (earlier emergence due to warmer 
## incubation temps), which would inflat end of season size in snowmelt stream relative to 
## GW streams, leading to our finding of no relationship using the empirical data. See 
## exploratory analysis in DataTable.R script. To account for these issues, we instead can 
## simulate body size/weight trajectories for a fish of a given size given modeled trends 
## in growth rates across a gradient of groundwater input and then compare simulated end 
## of season weights to infer relative growth benefits of spring vs snowmelt fed streams.
## Really, we're comparing scope for growth across these systems...per discussion below.
## 
## There are some interesting ideas and discussion to pursue here with respect to water temp, 
## spawn timing, developmental rates, emergence timing and early growth conditions that provide 
## very stimulating fodder for a couple Discussion paragraphs. For example, while fish may be 
## able to enter GW streams early to spawn b/c no flooding, temps are SO cold that fry still 
## emerge later than in runoff streams that are warmer. I.e., they aren't entirely able to  
## compensate for cold temps by spawning early. However, the extended growing season in GW 
## streams appears to allow fish to catch up to counterparts in runoff streams such that fish
## are more or less the same size by the onset of winter, which would explain lack of 
## relationship between empirical end of season body size and GW index. 
##---------------------------------------------------------------------------------------##
