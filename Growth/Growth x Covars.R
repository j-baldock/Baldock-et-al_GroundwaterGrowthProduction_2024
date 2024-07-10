#==========================================================================================================================================#
# PROJECT NAME: Groundwater structures fish growth and production across a riverscape
# CREATOR: Jeffrey R. Baldock, WY Coop Unit, University of Wyoming, jbaldock@uwyo.edu 
# PURPOSE: Estimate effects of temperature and density on variation in YOY growth, 
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
dat2_gr <- dat %>% filter(!is.na(rgr_len_censor)) %>% mutate(zgw = c(scale(gw, center = T, scale = T)),
                                                             zdenslog = c(scale(density_log, center = T, scale = T)),
                                                             zcpuelog = c(scale(cpue_log, center = T, scale = T)),
                                                             zeffdlog = c(scale(effden_permeter_log, center = T, scale = T)),
                                                             ztemp = c(scale(meantempper, center = T, scale = T)),
                                                             zpriorwt = c(scale(priorwtmean, center = T, scale = T)),
                                                             zpriorlen = c(scale(priorlenmean, center = T, scale = T)),
                                                             zpvalue = c(scale(pvalue, center = T, scale = T)))

# create table of covariate means and stdevs for back calculating/plotting
scaletbl_gr <- data.frame(var = c("gw", "denslog", "cpuelog", "effdlog", "temp", "priorwt", "priorlen", "pvalue"),
                          means = c(mean(dat2_gr$gw), mean(dat2_gr$density_log), mean(dat2_gr$cpue_log), mean(dat2_gr$effden_permeter_log), mean(dat2_gr$meantempper), mean(dat2_gr$priorwtmean), mean(dat2_gr$priorlenmean), mean(dat2_gr$pvalue)),
                          sds = c(sd(dat2_gr$gw), sd(dat2_gr$density_log), sd(dat2_gr$cpue_log), sd(dat2_gr$effden_permeter_log), sd(dat2_gr$meantempper), sd(dat2_gr$priorwtmean), sd(dat2_gr$priorlenmean), sd(dat2_gr$pvalue)))
write_csv(scaletbl_gr, "Growth by Covars/GrowthByCovars_z-score_table.csv")

# preliminary models to select most appropriate density metric
# summary(lm(rgr_len_censor ~ zdenslog + ztemp + I(ztemp^2) + zdenslog*ztemp, dat2))
# summary(lm(rgr_len_censor ~ zcpuelog + ztemp + I(ztemp^2) + zcpuelog*ztemp, dat2))
# summary(lm(rgr_len_censor ~ zeffdlog + ztemp + I(ztemp^2) + zeffdlog*ztemp, dat2))
# 
# # pairs plots
# boxplot(dat2 %>% dplyr::select(rgr_len_censor, zgw, zdenslog, zcpuelog, zeffdlog, ztemp))
# ggpairs(dat2 %>% dplyr::select(rgr_len_censor, zgw, zdenslog, zcpuelog, zeffdlog, ztemp))
# ggpairs(dat2 %>% dplyr::select(density, cpue, effden_permeter, density_log, cpue_log, effden_permeter_log))
# All density metrics, raw and logged, are highly correlated (raw r > 0.7 and logged r > 0.82)
# Adjusted R2 values are greater for density and cpue than for effective density
# Density (logged) is probably the most intuitive metric 
# Early objective will be to describe how groundwater induces characteristic patterns of water temp and density
# So it makes sense to use density (logged) for modeling growth, to maintain consistency and interpretability


# correlation between relative growth rates based on length vs weight
# growth rates based on mean size
# plot(rgr_wt_censor ~ rgr_len_censor, dat2)
# cor(dat2$rgr_wt_censor, dat2$rgr_len_censor)


##------------------------------------------------##
## Misc JAGS fields
##------------------------------------------------##

# grouping variables as factors for JAGS model
dat2_gr$stream <- as.factor(dat2_gr$stream)
dat2_gr$sectid <- as.factor(dat2_gr$sectid)
dat2_gr$strid <- as.factor(dat2_gr$strid)
dat2_gr$yrid <- as.factor(dat2_gr$yrid)

# tibble linking each section to each stream for JAGS model with section nested within stream
unq_sect_str_gr <- dat2_gr %>% group_by(strid) %>% reframe(sectid = unique(sectid)) %>% ungroup()

# some additional data for JAGS
nYears_gr <- length(unique(dat2_gr$yrid))
nSections_gr <- length(unique(dat2_gr$sectid))
nStreams_gr <- length(unique(dat2_gr$strid))

#
StreamNames_gr <- dat2_gr %>% group_by(stream) %>% summarize(strid = unique(strid), gw = mean(gw)) %>% arrange(gw) %>% ungroup() %>% 
  add_row(stream = c("dum1","dum2"), gw = c(0,1)) %>% # add dummy rows to ensure right color scheme
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>%
  filter(!stream %in% c("dum1","dum2")) %>%
  mutate(strid2 = 1:13)

# Sort 
dat2_gr <- dat2_gr %>% arrange(yrid, strid, sectid)


##---------------------------------------------------------------------------------------##
## Fit Candidate JAGS Models
##---------------------------------------------------------------------------------------##

# Given our study design (sampling shortly after/during emergence, short sampling intervals and thus overlapping size distributions among sampling events),
# much of the (bootstrapped) variation in growth rates likely stems from variability in emergence timing rather than individual variation in growth rates.
# B/c of overlapping distributions, bootstrapped measures of uncertainty are very large. Therefore, accounting for bootstrapped measures of 
# uncertainty with a state-space model formulation likely results in overly conservative estimates of the effects of covariates on growth, which may limit
# our ability to draw inferences about covariate effects AND explore residual variation in growth among streams that is not already account for in the model.
# Thus, running a model on the means (i.e., no state-space structure) is a justifiable path forward. 

##### Length-based growth

# specify candidate models, covariates on sigma are identical
Covs_list <- list(
  model.matrix(~ zdenslog + ztemp + I(ztemp^2) + zdenslog*ztemp + zpriorlen, data = dat2_gr, na.rm = FALSE),
  model.matrix(~ zdenslog + ztemp + zdenslog*ztemp + zpriorlen, data = dat2_gr, na.rm = FALSE),
  model.matrix(~ zdenslog + ztemp + zpriorlen, data = dat2_gr, na.rm = FALSE),
  model.matrix(~ zdenslog + zpriorlen, data = dat2_gr, na.rm = FALSE),
  model.matrix(~ ztemp + I(ztemp^2) + zpriorlen, data = dat2_gr, na.rm = FALSE),
  model.matrix(~ ztemp + zpriorlen, data = dat2_gr, na.rm = FALSE),
  model.matrix(~ zpriorlen, data = dat2_gr, na.rm = FALSE)
)

# fit candidate JAGS models in for loop
Fit_list_gr <- list() # empty list to store model fits
Rhat_list <- list()
for (i in 1:length(Covs_list)) {
  Covs_Sigma <- model.matrix(~ zgw, data = dat2_gr, na.rm = FALSE) # covariates on Sigma
  if (i %in% c(1:6)) {
    jags.data <- list("nObs" = dim(Covs_list[[i]])[1], "L" = dat2_gr$lenmean, "pL" = dat2_gr$priorlenmean, "elap" = dat2_gr$elapdays,                                                                       
                      "strid" = dat2_gr$strid, "yrid" = dat2_gr$yrid, "sectid" = dat2_gr$sectid, "sect_per_str" = unq_sect_str_gr$strid,                                                                        
                      "nYears" = length(unique(dat2_gr$yrid)), "nStreams" = length(unique(dat2_gr$strid)), "nSections" = length(unique(dat2_gr$sectid)),                                                
                      "Covs" = Covs_list[[i]][,-1], "nCovs" = dim(Covs_list[[i]])[2]-1, 
                      "CovsMu" = apply(Covs_list[[i]][,-1], MARGIN = 2, FUN = mean, na.rm = T), "CovsSD" = apply(Covs_list[[i]][,-1], MARGIN = 2, FUN = sd, na.rm = T),
                      "Covs_Sigma" = matrix(Covs_Sigma[,-1]), "nSigmaCovs" = dim(Covs_Sigma)[2]-1)
  } else {
    jags.data <- list("nObs" = dim(Covs_list[[i]])[1], "L" = dat2_gr$lenmean, "pL" = dat2_gr$priorlenmean, "elap" = dat2_gr$elapdays,                                                                       
                      "strid" = dat2_gr$strid, "yrid" = dat2_gr$yrid, "sectid" = dat2_gr$sectid, "sect_per_str" = unq_sect_str_gr$strid,                                                                        
                      "nYears" = length(unique(dat2_gr$yrid)), "nStreams" = length(unique(dat2_gr$strid)), "nSections" = length(unique(dat2_gr$sectid)),                                                
                      "Covs" = matrix(Covs_list[[i]][,-1]), "nCovs" = dim(Covs_list[[i]])[2]-1, 
                      "CovsMu" = apply(matrix(Covs_list[[i]][,-1]), MARGIN = 2, FUN = mean, na.rm = T), "CovsSD" = apply(matrix(Covs_list[[i]][,-1]), MARGIN = 2, FUN = sd, na.rm = T),
                      "Covs_Sigma" = matrix(Covs_Sigma[,-1]), "nSigmaCovs" = dim(Covs_Sigma)[2]-1)
  }
  params <- c("alpha.adj", "alpha.year.adj", "alpha.stream.adj", "alpha.sect.adj", "beta", "gr", "Lexp", "L",
              "sigma.gr", "log.sigma.alpha", "sigma.alpha", "log.sigma.beta", "sigma.sb", 
              "sigma.yr", "sigma.str", "sigma.sect", "sigma.a", "sigma.b", 
              "loglik", "res", "var_fit", "var_fix", "var_res", "var_yr", "var_sect", "var_str", "margR2", "condR2")
  m <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "JAGS Models/GrowthModel_Letcher.txt", 
                     n.chains = 3, n.thin = 20, n.burnin = 50000, n.iter = 150000, DIC = TRUE)
  Fit_list_gr[[i]] <- m
  MCMCtrace(m, ind = TRUE, filename = paste("Growth by Covars/Traceplots/Growth_Length_m", i, "_traceplots.pdf", sep = "")) # traceplots
  Rhat_list[[i]] <- m$BUGSoutput$summary[,8][m$BUGSoutput$summary[,8] > 1.01] # any R-hat values too high?  
  print(i)
}

# any R-hat values too high?  
Rhat_list


##---------------------------------------------------------------------------------------##
## MODEL SELECTION - compare LOO and set top model
##---------------------------------------------------------------------------------------##

loo_list <- list()
for (i in 1:length(Fit_list_gr)) {
  m <- Fit_list_gr[[i]]
  loglik <- m$BUGSoutput$sims.list$loglik
  reff <- relative_eff(exp(loglik), chain_id = c(rep(1,5000),rep(2,5000),rep(3,5000)))
  loo_list[[i]] <- loo(loglik, r_eff = reff)
  print(i)
}
lc <- loo_compare(loo_list)
print(lc, simplify = FALSE, digits = 2)
plot(loo_list[[2]])
write.csv(data.frame(lc), "/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/Growth by Covars/GrowthbyCovars_LOO_output.csv")

# set the top model based on LOO
top_mod <- Fit_list_gr[[2]]
Covs <- Covs_list[[2]]


##---------------------------------------------------------------------------------------##
## Save MCMC samples from top model (and second best model, top + cubic temp)
##---------------------------------------------------------------------------------------##

# Save MCMC samples from top model
top_mod.mcmc <- as.data.frame(as.matrix(as.mcmc(top_mod)))
modelout <- top_mod$BUGSoutput
# generate MCMC samples and store as an array
McmcList <- vector("list", length = dim(modelout$sims.array)[2])
for(i in 1:length(McmcList)) { McmcList[[i]] = as.mcmc(modelout$sims.array[,i,]) }
# rbind MCMC samples from 3 chains and save as object
Mcmcdat <- rbind(McmcList[[1]], McmcList[[2]], McmcList[[3]])
param.summary <- modelout$summary

# save model output
write_csv(as.data.frame(Mcmcdat), "Growth by Covars/GrowthByCovars_TopModel_mcmcsamps.csv")
Mcmcdat <- read.csv("Growth by Covars/GrowthByCovars_TopModel_mcmcsamps.csv")
write.csv(as.data.frame(param.summary), "Growth by Covars/GrowthByCovars_TopModel_ParameterSummary.csv", row.names = T)
param.summary <- read.csv("Growth by Covars/GrowthByCovars_TopModel_ParameterSummary.csv", row.names = 1)

# Mcmcdat_full <- read.csv("Growth by Covars/GrowthByCovars_FullModel_mcmcsamps.csv")


##---------------------------------------------------------------------------------------##
## MODEL DIAGNOSTICS
##---------------------------------------------------------------------------------------##

# subset expected and observed MCMC samples
ppdat_exp <- as.matrix(Mcmcdat[,startsWith(names(Mcmcdat), "Lexp.")])
ppdat_obs <- as.matrix(Mcmcdat[,startsWith(names(Mcmcdat), "L.")])

# Bayesian p-value
sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2])

# Posterior Predictive Check: plot median posterior expected length ~ observed length with Bayesian p-value
jpeg("Growth by Covars/Figures/GrowthByCovars_PPCheck.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(x = seq(from = 25, to = 75, length.out = 100), y = seq(from = 25, to = 75, length.out = 100), pch = NA, xlab = "Observed Length (mm)", ylab = "Median Posterior Expected Length (mm)")
for (i in 1:dim(Covs)[1]) { points(median(Mcmcdat[,paste("Lexp.",i,".", sep = "")]) ~ dat2_gr$lenmean[i])}
legend("topleft", bty = "n", legend = paste("Bayesian p-value = ", round(sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2]), digits = 3), sep = ""))
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()

# histogram of residuals
jpeg("Growth by Covars/Figures/GrowthByCovars_ResidHist.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist((ppdat_obs - ppdat_exp), main = "", xlab = "Observed - Expected")
legend("topright", bty = "n", legend = paste("Median Resid. = ", round(median((ppdat_obs - ppdat_exp)), digits = 2), " mm", sep = ""))
abline(v = median(unlist(ppdat_obs - ppdat_exp)), col = "red", lwd = 2)
dev.off()

# RMSE
sqrt(mean((ppdat_obs - ppdat_exp)^2)) # point estimate
rmse <- vector("numeric", 15000L)
for (j in 1:dim(Mcmcdat)[1]) {
  res <- vector("numeric", 258L)
  for (i in 1:dim(Covs_list[[1]])[1]) {
    res[i] <- (ppdat_obs[j,i] - ppdat_exp[j,i])^2
  }
  rmse[j] <- sqrt(mean(res))
  print(j)
}
jpeg("Growth by Covars/Figures/GrowthByCovars_RMSE.jpg", units = "in", width = 5, height = 4.5, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist(rmse, xlab = "Root Mean Square Error (mm)", ylab = "Posterior Frequency", main = "")
legend("topright", bty = "n", legend = paste("Median RMSE = ", round(median(rmse), digits = 2), " mm", sep = ""))
abline(v = median(rmse), col = "red", lwd = 2)
dev.off()


# Calculate Bayesian R^2 per Gelman et al. (2018), Nakagawa and Shielzeth (2013)
# Note: this does not work (as is) for state-space model formulation...conditional and marginal both approach 1
# Note: this is "conditional" R^2...that is, variance explained by both fixed and random effects combined. 
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
hist((var_fix + var_yr + var_str + var_sect) / (var_fix + var_res + var_yr + var_str + var_sect))
hist(var_fit / (var_fit + var_res), add = T)
median((var_fix + var_yr + var_str + var_sect) / (var_fix + var_res + var_yr + var_str + var_sect))
median(var_fit / (var_fit + var_res))

# random components each
hist((var_str) / (var_fix + var_str + var_yr + var_sect + var_res))
hist((var_sect) / (var_fix + var_str + var_yr + var_sect + var_res))
hist((var_yr) / (var_fix + var_str + var_yr + var_sect + var_res))
# ### Very broad distribution of variation explained by year component likely driven by few (2) factors associated with that hierarchical component.
# ### Multiple references (Gelman and Hill???) note that hierarchical variance may be poorly estimated (overly broad) when number of factors is small. 
# ### This results in distributions of marginal and conditional R^2 that are also very broad. Thus, should rely on modes for inference, as the mode 
# ### (more so than the median) reduces the influence of overly broad year component. Modes align well with preliminary analysis using frequentist
# ### approach (lmer), which yields marginal R2 = 0.401 and conditional R2 = 0.495
# hist((var_fix) / (var_fix + var_str + var_yr + var_sect + var_res)) # marginal R2 with year component
# hist(top_mod$BUGSoutput$sims.list$margR2)
# hist((var_fix) / (var_fix + var_str + var_sect + var_res)) # demonstrate how dropping the year component results in more normal distribution

# function to calculate mode
# getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Marginal/conditional plot
jpeg("Growth by Covars/Figures/GrowthByCovars_MargCondR2.jpg", units = "in", width = 5.5, height = 5, res = 1500)
par(mar = c(4,4,0.5,0.5), mgp = c(2.5,1,0))
hist(Mcmcdat[,"margR2"], col = alpha("darkorange", 0.5), xlim = c(0,1), ylim = c(0,7000), breaks = seq(from = 0, to = 1, by = 0.02), main = NA, xlab = expression(Bayesian ~ R^2), ylab = "Posterior Frequency")
hist(Mcmcdat[,"condR2"], col = alpha("forestgreen", 0.5), xlim = c(0,1), add = TRUE, breaks = seq(from = 0, to = 1, by = 0.02))
# abline(v = median(Mcmcdat[,"margR2"]), col = "darkorange", lwd = 3)
# abline(v = median(Mcmcdat[,"condR2"]), col = "forestgreen", lwd = 3)
legend("topleft", fill = c(alpha("darkorange", 0.5), alpha("forestgreen", 0.5)), bty = "n", cex = 0.9,
       legend = c(paste("Marginal = ", round(median(Mcmcdat[,"margR2"]), 3), " (", round(hdi(Mcmcdat[,"margR2"], credMass = 0.95)[1], 3), ", ", round(hdi(Mcmcdat[,"margR2"], credMass = 0.95)[2], 3), ")", sep = ""), 
                  paste("Conditional = ", round(median(Mcmcdat[,"condR2"]), 3), " (", round(hdi(Mcmcdat[,"condR2"], credMass = 0.95)[1], 3), ", ", round(hdi(Mcmcdat[,"condR2"], credMass = 0.95)[2], 3), ")", sep = "")))
dev.off()


##### Combined Figure
jpeg("Growth by Covars/Figures/GrowthByCovars_ModelDiagnosticsCombined.jpg", units = "in", width = 7, height = 7, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), mfrow = c(2,2))

# PP Check
plot(x = seq(from = 25, to = 75, length.out = 100), y = seq(from = 25, to = 75, length.out = 100), pch = NA, xlab = "Observed Length (mm)", ylab = "Median Posterior Expected Length (mm)", bty = "l")
for (i in 1:dim(Covs)[1]) { points(median(Mcmcdat[,paste("Lexp.",i,".", sep = "")]) ~ dat2_gr$lenmean[i])}
legend("top", bty = "n", legend = paste("Bayesian p-value = ", round(sum(ppdat_exp > ppdat_obs) / (dim(ppdat_obs)[1]*dim(ppdat_obs)[2]), digits = 3), sep = ""))
abline(a = 0, b = 1, col = "red", lwd = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(a)")
par(usr = usr)

# RMSE
hist(rmse, xlab = "Root Mean Square Error (mm)", ylab = "Posterior Frequency", main = "", bty = "l")
box(bty = "l")
legend("topright", bty = "n", legend = paste("Median = ", round(median(rmse), digits = 2), " mm", sep = ""))
abline(v = median(rmse), col = "red", lwd = 2)
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(b)")
par(usr = usr)

# Marg/Cond R2
hist(Mcmcdat[,"margR2"], col = alpha("darkorange", 0.5), xlim = c(0,1), ylim = c(0,7500), breaks = seq(from = 0, to = 1, by = 0.02), main = NA, xlab = expression(Bayesian ~ R^2), ylab = "Posterior Frequency", bty = "l")
box(bty = "l")
hist(Mcmcdat[,"condR2"], col = alpha("forestgreen", 0.5), xlim = c(0,1), add = TRUE, breaks = seq(from = 0, to = 1, by = 0.02))
legend("topright", fill = c(alpha("darkorange", 0.5), alpha("forestgreen", 0.5)), bty = "n", cex = 0.9,
       legend = c(paste("Marginal = ", round(median(Mcmcdat[,"margR2"]), 3), " (", round(hdi(Mcmcdat[,"margR2"], credMass = 0.95)[1], 3), ", ", round(hdi(Mcmcdat[,"margR2"], credMass = 0.95)[2], 3), ")", sep = ""), 
                  paste("Conditional = ", round(median(Mcmcdat[,"condR2"]), 3), " (", round(hdi(Mcmcdat[,"condR2"], credMass = 0.95)[1], 3), ", ", round(hdi(Mcmcdat[,"condR2"], credMass = 0.95)[2], 3), ")", sep = "")))
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(c)")
par(usr = usr)

# residuals by Groundwater
plot(apply(as.matrix(Mcmcdat[,startsWith(names(Mcmcdat), "res")]), 2, mean) ~ dat2_gr$gw, xlab = "Groundwater Index", ylab = "Residuals", bty = "l")
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.07, 0.95, labels = "(d)")
par(usr = usr)

dev.off()

##---------------------------------------------------------------------------------------##
## PLOTTING 
##---------------------------------------------------------------------------------------##

# Parameter (beta) estimates - dot plot
mod.gg <- ggs(as.mcmc(top_mod))
Covs <- Covs

# beta - primary covariates
jpeg("Growth by Covars/Figures/GrowthByCovars_ParameterDotPlot_beta.jpg", units = "in", width = 5, height = 4, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^beta", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = colnames(Covs)[-1])
dev.off()

# log.sigma.beta - covariates on sigma
jpeg("Growth by Covars/Figures/GrowthByCovars_ParameterDotPlot_log.sigma.beta.jpg", units = "in", width = 5, height = 4, res = 1500)
ggs_caterpillar(D = mod.gg, family = "log.sigma.beta", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) +
  theme_bw() + ylab("") + xlab("Posterior estimate") +
  scale_y_discrete(labels = "zgw")
dev.off()

# alpha.stream - stream offsets to the intercept
StreamNames_gr <- StreamNames_gr %>% arrange(strid)
jpeg("Growth by Covars/Figures/GrowthByCovars_ParameterDotPlot_alpha.stream.jpg", units = "in", width = 5, height = 4, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^alpha.stream.adj", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) +
  theme_bw() + ylab("") + xlab("Posterior estimate") +
  scale_y_discrete(labels = StreamNames_gr$stream)
dev.off()

# alpha.section - section offsets to the intercept
jpeg("Growth by Covars/Figures/GrowthByCovars_ParameterDotPlot_alpha.section.jpg", units = "in", width = 5, height = 9, res = 1500)
ggs_caterpillar(D = mod.gg, family = "^alpha.sect.adj", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = unq_sect_str_gr$sectid)
dev.off()

ggs_caterpillar(D = mod.gg, family = "^alpha.year.adj", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + 
  scale_y_discrete(labels = c(2021:2022))

#-------------------------------------------#
# Random Intercepts: Stream x Year
#-------------------------------------------#

StreamNames_gr <- StreamNames_gr %>% arrange(strid2)

jpeg("Growth by Covars/Figures/GrowthByCovars_PosteriorIntercepts.jpg", units = "in", width = 3.75, height = 6, res = 1500)
par(mar = c(4,3,1.5,0.5), mgp = c(2.5,1,0))
# set up plot area
x <- seq(from = 0.25, to = 0.48, length.out = 100)
y <- seq(from = 1, to = 13.7, length.out = 100)
plot(y ~ x, type = "n", xlab = expression(paste("Mean Growth Rate (mm day"^"-1", ")", sep = "")), ylab = "", axes = F)
axis(1)
# iterate over intercepts
ctr <- 0
for (j in 1:nStreams_gr) {
  par(xpd = NA)
  text(y = j+0.5, x = 0.2, labels = StreamNames_gr$stream[j], pos = 4)
  par(xpd = FALSE)
  dens <- density(Mcmcdat[,"alpha.adj"] + Mcmcdat[,paste("alpha.stream.adj.",StreamNames_gr$strid[j],".", sep = "")])
  dens$y2 <- (dens$y / max(dens$y))
  l <- min(which(dens$x >= hdi(dens, credMass = 0.95)[1]))
  h <- max(which(dens$x < hdi(dens, credMass = 0.95)[2]))
  polygon(c(dens$x[c(l, l:h, h)]), c(0+j, dens$y2[l:h]+j, 0+j), col = adjustcolor(StreamNames_gr$cols[j], alpha.f = 0.8), lty = 0)
  lines(dens$y2+j ~ dens$x, lwd = 1.5)
  ctr <- ctr + 4
}
abline(v = median(Mcmcdat[,"alpha.adj"]), lty = 2)
par(xpd = TRUE)
# legend(x = 0.25, y = 15.2, legend = c(2021, 2022), lty = 1:2, lwd = 1.5, bty = "n", horiz = TRUE, cex = 0.9)
dev.off()


jpeg("Growth by Covars/Figures/GrowthByCovars_PosteriorIntercepts_withYears.jpg", units = "in", width = 3.75, height = 6, res = 1500)
par(mar = c(4,3,1.5,0.5), mgp = c(2.5,1,0))
# set up plot area
x <- seq(from = 0.2, to = 0.5, length.out = 100)
y <- seq(from = 1, to = 13.7, length.out = 100)
plot(y ~ x, type = "n", xlab = expression(paste("Mean Growth Rate (mm day"^"-1", ")", sep = "")), ylab = "", axes = F)
axis(1)
# iterate over intercepts
ctr <- 0
for (j in 1:nStreams_gr) {
  par(xpd = NA)
  text(y = j+0.5, x = 0.15, labels = StreamNames_gr$stream[j], pos = 4)
  for (i in 1:nYears_gr){
    par(xpd = FALSE)
    dens <- density(Mcmcdat[,"alpha.adj"] + Mcmcdat[,paste("alpha.year.adj.",i,".", sep = "")] + Mcmcdat[,paste("alpha.stream.adj.",StreamNames_gr$strid[j],".", sep = "")])
    dens$y2 <- (dens$y / max(dens$y))
    l <- min(which(dens$x >= hdi(dens, credMass = 0.95)[1]))
    h <- max(which(dens$x < hdi(dens, credMass = 0.95)[2]))
    polygon(c(dens$x[c(l, l:h, h)]), c(0+j, dens$y2[l:h]+j, 0+j), col = adjustcolor(StreamNames_gr$cols[j], alpha.f = 0.5), lty = 0)
    lines(dens$y2+j ~ dens$x, lwd = 1.5, lty = i)
  }
  ctr <- ctr + 4
}
abline(v = median(Mcmcdat[,"alpha.adj"]))
par(xpd = TRUE)
legend(x = 0.25, y = 15.2, legend = c(2021, 2022), lty = 1:2, lwd = 1.5, bty = "n", horiz = TRUE, cex = 0.9)
dev.off()


jpeg("Growth by Covars/Figures/GrowthByCovars_PosteriorIntercepts_withSections.jpg", units = "in", width = 3.75, height = 6, res = 1500)
par(mar = c(4,3,1.5,0.5), mgp = c(2.5,1,0))
# set up plot area
x <- seq(from = 0.2, to = 0.5, length.out = 100)
y <- seq(from = 1, to = 13.7, length.out = 100)
plot(y ~ x, type = "n", xlab = expression(paste("Mean Growth Rate (mm day"^"-1", ")", sep = "")), ylab = "", axes = F)
axis(1)
# iterate over intercepts
ctr <- 0
dd <- unq_sect_str_gr %>% mutate(sectid2 = c(1:dim(unq_sect_str_gr)[1]))
for (j in 1:nStreams_gr) {
  par(xpd = NA)
  text(y = j+0.5, x = 0.15, labels = StreamNames_gr$stream[j], pos = 4)
  ddd <- dd %>% filter(strid == j)
  for (i in 1:nYears_gr){
    for (k in 1:dim(ddd)[1]) {
      par(xpd = FALSE)
      dens <- density(Mcmcdat[,"alpha.adj"] + Mcmcdat[,paste("alpha.year.adj.",i,".", sep = "")] + Mcmcdat[,paste("alpha.sect.adj.",as.numeric(ddd$sectid2[k]),".", sep = "")])
      dens$y2 <- (dens$y / max(dens$y))
      l <- min(which(dens$x >= hdi(dens, credMass = 0.95)[1]))
      h <- max(which(dens$x < hdi(dens, credMass = 0.95)[2]))
      polygon(c(dens$x[c(l, l:h, h)]), c(0+j, dens$y2[l:h]+j, 0+j), col = adjustcolor(StreamNames_gr$cols[j], alpha.f = 0.3), lty = 0)
      lines(dens$y2+j ~ dens$x, lwd = 1.5, lty = i)
    }
  }
  ctr <- ctr + 4
}
abline(v = median(Mcmcdat[,"alpha.adj"]))
par(xpd = TRUE)
legend(x = 0.25, y = 15.2, legend = c(2021, 2022), lty = 1:2, lwd = 1.5, bty = "n", horiz = TRUE, cex = 0.9)
dev.off()


#-------------------------------------------#
# Marginal Effects Plots
#-------------------------------------------#

range(dat2_gr$ztemp)
dd <- dat2_gr[dat2_gr$zdenslog >= quantile(dat2_gr$zdenslog, 0.75),]
range(dd$ztemp) # -0.823, 0.891
dd <- dat2_gr[dat2_gr$zdenslog <= quantile(dat2_gr$zdenslog, 0.25),]
range(dd$ztemp)


nvalues <- 100
# zsCols <- rev(hcl.colors(2, "Zissou 1"))
zsPal <- colorRampPalette(hcl.colors(10, "Zissou 1"))

# generate sequences of predictor variables and scales axes
# x.gw <- seq(from = min(dat2$zgw), to = max(dat2$zgw), length.out = nvalues)
x.temp <- seq(from = min(dat2_gr$ztemp), to = max(dat2_gr$ztemp), length.out = nvalues)
x.temp2 <- seq(from = -0.823, to = 0.891, length.out = nvalues)
x.logdens <- seq(from = min(dat2_gr$zdenslog), to = max(dat2_gr$zdenslog), length.out = nvalues)
# x.priorlen <- seq(from = min(dat2$zpriorlen), to = max(dat2$zpriorlen), length.out = nvalues)
# manually set up and transform axes to original scale of the data
# x.axis.gw <- seq(from = 0, to = 0.8, by = 0.1)
# x.scaled.gw <- (x.axis.gw - scaletbl$means[1]) / scaletbl$sds[1]
x.axis.logdens <- seq(from = -3, to = 1, by = 1)
x.scaled.logdens <- (x.axis.logdens - scaletbl_gr$means[2]) / scaletbl_gr$sds[2]
x.axis.temp <- seq(from = 4, to = 16, by = 2)
x.scaled.temp <- (x.axis.temp - scaletbl_gr$means[5]) / scaletbl_gr$sds[5]
# x.axis.priorlen <- seq(from = 20, to = 70, by = 10)
# x.scaled.priorlen <- (x.axis.priorlen - scaletbl$means[7]) / scaletbl$sds[7]

jpeg("Growth by Covars/Figures/GrowthByCovars_MarginalEffects_TempFocus.jpg", units = "in", width = 3.75, height = 3.75, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# growth ~ temp + density (temp focal)
plot(seq(from = -0.2, to = 1, length.out = length(dat2_gr$ztemp)) ~ ztemp, dat2_gr, pch = NA, axes = F, 
     #col = alpha(zsPal(10)[as.numeric(cut(dat2$zdenslog, breaks = 10))], 0.7),
     xlab = expression(paste("Mean Temperature ("^"o", "C)", sep = "")), ylab = expression(paste("Growth Rate (mm day"^"-1", ")", sep = "")))
axis(1, at = x.scaled.temp, labels = x.axis.temp)
axis(2)
box(bty = "l")
# interaction: temp * high density
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*max(x.logdens) + Mcmcdat[i,"beta.2."]*x.temp2 + Mcmcdat[i,"beta.4."]*max(x.logdens)*x.temp2 }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.temp2, rev(x.temp2)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty=0)
lines(pred_median ~ x.temp2, lty = 2, lwd = 2)
# interaction: temp * low density
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*min(x.logdens) + Mcmcdat[i,"beta.2."]*x.temp + Mcmcdat[i,"beta.4."]*min(x.logdens)*x.temp }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.temp, rev(x.temp)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty=0)
lines(pred_median ~ x.temp, lty = 1, lwd = 2)
# points
for (i in 1:dim(dat2_gr)[1]) {
  points(median(Mcmcdat[,paste("gr.",i,".", sep = "")]) ~ dat2_gr$ztemp[i], pch = 16, col = alpha(zsPal(10)[as.numeric(cut(dat2_gr$zdenslog, breaks = 10))], 0.5)[i])
}
# legends
usr <- par("usr")
par(usr = c(0,1,0,1))
legend(x = 0, y = 0.675, legend = c("Min.", "Max."), lwd = 2, lty = 1:2, bty = "n", cex = 0.8)
# point legend
text(x = 0.15, y = 0.97, labels = "log(Density)", cex = 0.8)
gradientLegend(valRange = c(min(dat2_gr$zdenslog), max(dat2_gr$zdenslog)), color = alphaPalette(zsPal(100), f.seq = rep(0.9,100)), length = 0.5, inside = TRUE, side = 2, n.seg = 1, pos = c(0.05, 0.7, 0.1, 0.9), dec = 2)
par(usr = usr)
dev.off()


jpeg("Growth by Covars/Figures/GrowthByCovars_MarginalEffects_DensityFocus.jpg", units = "in", width = 3.75, height = 3.75, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# growth ~ temp + density (density focal)
plot(seq(from = -0.2, to = 1, length.out = length(dat2_gr$zdenslog)) ~ zdenslog, dat2_gr, pch = NA, axes = F, 
     #col = alpha(zsPal(10)[as.numeric(cut(dat2$zdenslog, breaks = 10))], 0.7),
     xlab = "log(Density)", ylab = expression(paste("Growth Rate (mm day"^"-1", ")", sep = "")))
axis(1, at = x.scaled.logdens, labels = x.axis.logdens)
axis(2)
box(bty = "l")
# interaction: temp * high density
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*(x.logdens) + Mcmcdat[i,"beta.2."]*max(x.temp) + Mcmcdat[i,"beta.4."]*(x.logdens)*max(x.temp) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.logdens, rev(x.logdens)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty=0)
lines(pred_median ~ x.logdens, lty = 2, lwd = 2)
# interaction: temp * low density
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*(x.logdens) + Mcmcdat[i,"beta.2."]*min(x.temp) + Mcmcdat[i,"beta.4."]*(x.logdens)*min(x.temp) }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.logdens, rev(x.logdens)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty=0)
lines(pred_median ~ x.logdens, lty = 1, lwd = 2)
# points
for (i in 1:dim(dat2_gr)[1]) {
  points(median(Mcmcdat[,paste("gr.",i,".", sep = "")]) ~ dat2_gr$zdenslog[i], pch = 16, col = alpha(zsPal(10)[as.numeric(cut(dat2_gr$ztemp, breaks = 10))], 0.5)[i])
}

# legend("bottom", legend = c("Min. Temp", "Max. Temp"), lwd = 2, lty = 1:2, cex = 0.9, bty = "n")
# legend(x = -1.4, y = -0.10, legend = c("", ""), pch = 16, col = zsPal(2), bty = "n", cex = 0.9)
legend("top", legend = c("Low temp", "High temp"), lwd = 2, lty = 1:2, bty = "n", cex = 0.8)

dev.off()


############ For talk
png("Growth by Covars/Figures/Presentation Figs/GrowthByCovars_MarginalEffects_PointsOnly.png", units = "in", width = 3.75, height = 3.75, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), bg = NA)
plot(seq(from = -0.2, to = 1, length.out = length(dat2_gr$ztemp)) ~ ztemp, dat2_gr, pch = NA, axes = F, xlab = expression(paste("Mean Temperature ("^"o", "C)", sep = "")), ylab = expression(paste("Growth Rate (mm day"^"-1", ")", sep = "")))
axis(1, at = x.scaled.temp, labels = x.axis.temp)
axis(2)
box(bty = "l")
# points
for (i in 1:dim(dat2_gr)[1]) { points(median(Mcmcdat[,paste("gr.",i,".", sep = "")]) ~ dat2_gr$ztemp[i], pch = 16, col = alpha(zsPal(10)[as.numeric(cut(dat2_gr$zdenslog, breaks = 10))], 0.5)[i]) }
# legends
usr <- par("usr")
par(usr = c(0,1,0,1))
legend(x = 0, y = 0.675, legend = c("Min.", "Max."), lwd = 2, lty = 1:2, bty = "n", cex = 0.8)
# point legend
text(x = 0.15, y = 0.97, labels = "log(Density)", cex = 0.8)
gradientLegend(valRange = c(min(dat2_gr$zdenslog), max(dat2_gr$zdenslog)), color = alphaPalette(zsPal(100), f.seq = rep(0.9,100)), length = 0.5, inside = TRUE, side = 2, n.seg = 1, pos = c(0.05, 0.7, 0.1, 0.9), dec = 2)
par(usr = usr)
dev.off()


png("Growth by Covars/Figures/Presentation Figs/GrowthByCovars_MarginalEffects_MaxReg.png", units = "in", width = 3.75, height = 3.75, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), bg = NA)
plot(seq(from = -0.2, to = 1, length.out = length(dat2_gr$ztemp)) ~ ztemp, dat2_gr, pch = NA, axes = F, xlab = "", ylab = "")
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*max(x.logdens) + Mcmcdat[i,"beta.2."]*x.temp2 + Mcmcdat[i,"beta.4."]*max(x.logdens)*x.temp2 }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.temp2, rev(x.temp2)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty=0)
lines(pred_median ~ x.temp2, lty = 2, lwd = 2)
dev.off()


png("Growth by Covars/Figures/Presentation Figs/GrowthByCovars_MarginalEffects_MinReg.png", units = "in", width = 3.75, height = 3.75, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0), bg = NA)
plot(seq(from = -0.2, to = 1, length.out = length(dat2_gr$ztemp)) ~ ztemp, dat2_gr, pch = NA, axes = F, xlab = "", ylab = "")
pred_dist <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nvalues)
for (i in 1:nrow(pred_dist)) { pred_dist[i,] <- Mcmcdat[i,"alpha.adj"] + Mcmcdat[i,"beta.1."]*min(x.logdens) + Mcmcdat[i,"beta.2."]*x.temp + Mcmcdat[i,"beta.4."]*min(x.logdens)*x.temp }
pred_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.5)
polygon(c(x.temp, rev(x.temp)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.2), lty=0)
lines(pred_median ~ x.temp, lty = 1, lwd = 2)
dev.off()



#-------------------------------------------#
# Groundwater effect on sigma
#-------------------------------------------#

# number of values for prediction
nvalues <- 100
# generate sequences of predictor variables and scales axes
x.gw <- seq(from = min(dat2_gr$zgw), to = max(dat2_gr$zgw), length.out = nvalues)
x.axis.gw <- seq(from = 0, to = 0.8, by = 0.2)
x.scaled.gw <- (x.axis.gw - scaletbl_gr$means[1]) / scaletbl_gr$sds[1]

# Model without prior length
jpeg("Growth by Covars/Figures/GrowthByCovars_GWEffectOnSigma.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# growth ~ temp + density (temp focal)
plot(seq(from = 0.14, to = 0.24, length.out = length(dat2_gr$zgw)) ~ zgw, dat2_gr, pch = NA, xlab = "Groundwater Influence", ylab = "Standard Deviation of Growth Rate", axes = F, xlim = c(min(x.scaled.gw), max(x.scaled.gw)))
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
# for (i in 1:dim(dat2)[1]) { points(median(Mcmcdat1[,paste("sigma.gr.",i,".", sep = "")]) ~ dat2$zgw[i]) }

dev.off()



#-------------------------------------------#
# plot residuals and calculated growth rates against groundwater influence
#-------------------------------------------#

# residuals
jpeg("Growth by Covars/Figures/GrowthByCovars_ResidualsByGroundwater.jpg", units = "in", width = 4, height = 4, res = 1500)
par(mar = c(3.5,3.5,1,1), mgp = c(2,0.8,0))
plot(apply(top_mod$BUGSoutput$sims.list$res, 2, mean) ~ dat2_gr$gw, xlab = "Groundwater Influence", ylab = "Residuals", bty = "l")
dev.off()

# model-estimated growth rates
jpeg("Growth by Covars/Figures/GrowthByCovars_GrowthByGroundwater.jpg", units = "in", width = 4, height = 4, res = 1500)
par(mar = c(3.5,3.5,1,1), mgp = c(2,0.8,0))
plot(apply(top_mod$BUGSoutput$sims.list$gr, 2, mean) ~ dat2_gr$gw, xlab = "Groundwater Influence", ylab = expression(paste("Growth Rate (mm day"^"-1", ")", sep = "")), bty = "l")
dev.off()

# CV of growth by groundwater
grtib <- dat2_gr %>% select(stream, year, section, doy, gw) %>% mutate(growth = apply(top_mod$BUGSoutput$sims.list$gr, 2, mean))
cvtib <- grtib %>% group_by(stream, section, gw) %>% summarize(cvgr = sd(growth)/ mean(growth)) %>% ungroup()
cvmod <- lm(log(cvgr) ~ gw, cvtib)
summary(cvmod)
plot(cvmod)
x.gw <- seq(from = min(cvtib$gw), to = max(cvtib$gw), length.out = 100)
preds <- predict(cvmod, newdata = list(gw = x.gw), interval = "confidence")

jpeg("Growth by Covars/Figures/GrowthByCovars_CoefVarGrowthByGroundwater.jpg", units = "in", width = 4, height = 4, res = 1500)
par(mar = c(3.5,3.5,1,1), mgp = c(2,0.8,0))
plot(cvgr ~ gw, cvtib, xlab = "Groundwater Influence", ylab = "CV of Growth Rates", bty = "l", pch = NA)
polygon(x = c(x.gw, rev(x.gw)), y = exp(c(preds[,2], rev(preds[,3]))), col = alpha("black", 0.2), lty = 0)
lines(exp(preds[,1]) ~ x.gw, lwd = 2)
points(cvgr ~ gw, cvtib)
dev.off()









##---------------------------------------------------------------------------------------##
## Frequentist Approach
##---------------------------------------------------------------------------------------##
# dat2$zlogdens <- dat2$zeffdlog
# 
# # preliminary model
# grmod <- lm(specgrowthlencensor ~ zlogdens + ztemp + I(ztemp^2) + zlogdens*ztemp + zpriorlen, dat2)
# summary(grmod)
# # plot(grmod)
# summary(grmod)$adj.r.squared # model explains 33.7% percent of variation in growth rates...OK (no random effects here)
# sqrt(mean(grmod$residuals^2)) # RMSE of 1.9% body mass per day...not too bad
# car::vif(grmod, type = "predictor")
# 
# # PLOT
# # predict from the fitted model
# png("Figures/GrowthCovars.png", units = "in", width = 6, height = 6, res = 1000)
# par(mfrow = c(2,2), mar = c(4,4,1,1), mgp = c(2,0.8,0))
# newcols <- hcl.colors(2, "Geyser")
# nvals <- 100
# temps <- seq(from = min(dat2$ztemp, na.rm = T), to = max(dat2$ztemp, na.rm = T), length.out = nvals)
# denss <- seq(from = min(dat2$zlogdens, na.rm = T), to = max(dat2$zlogdens, na.rm = T), length.out = nvals)
# lens <- seq(from = min(dat2$zpriorlen, na.rm = T), to = max(dat2$zpriorlen, na.rm = T), length.out = nvals)
# 
# # growth ~ temp + density: gw = mean (temp focal)
# plot(specgrowthlen ~ ztemp, dat2, pch = NA, xlab = "Mean Temperature", ylab = "Growth Rate (% mass per day)")
# 
# nd <- list(zlogdens = rep(min(denss), times = nvals), ztemp = temps, zpriorlen = rep(0, times = nvals))
# preds <- predict(grmod, newdata = nd, interval = "confidence")
# polygon(c(temps, rev(temps)), c(preds[,2], rev(preds[,3])), col = scales::alpha(newcols[1], 0.2), lty=0)
# lines(preds[,1] ~ temps, col = newcols[1], lwd = 3)
# 
# nd <- list(zlogdens = rep(max(denss), times = nvals), ztemp = temps, zpriorlen = rep(0, times = nvals))
# preds <- predict(grmod, newdata = nd, interval = "confidence")
# polygon(c(temps, rev(temps)), c(preds[,2], rev(preds[,3])), col = scales::alpha(newcols[2], 0.2), lty=0)
# lines(preds[,1] ~ temps, col = newcols[2], lwd = 3)
# 
# legend("topleft", legend = c("Min Densiy", "Max Density"), fill = newcols, bty = "n")
# 
# # growth ~ temp + density: gw = mean (density focal)
# plot(specgrowthlen ~ zlogdens, dat2, pch = NA, xlab = "Conspecific YOY Density", ylab = "Growth Rate (% mass per day)")
# 
# nd <- list(zlogdens = denss, ztemp = rep(min(temps), times = nvals), zpriorlen = rep(0, times = nvals))
# preds <- predict(grmod, newdata = nd, interval = "confidence")
# polygon(c(denss, rev(denss)), c(preds[,2], rev(preds[,3])), col = scales::alpha(newcols[1], 0.2), lty=0)
# lines(preds[,1] ~ denss, col = newcols[1], lwd = 3)
# 
# nd <- list(zlogdens = denss, ztemp = rep(max(temps), times = nvals), zpriorlen = lens)
# preds <- predict(grmod, newdata = nd, interval = "confidence")
# polygon(c(denss, rev(denss)), c(preds[,2], rev(preds[,3])), col = scales::alpha(newcols[2], 0.2), lty=0)
# lines(preds[,1] ~ denss, col = newcols[2], lwd = 3)
# 
# legend("topright", legend = c("Min Temp", "Max Temp"), fill = newcols, bty = "n")
# 
# # growth ~ temp alone
# plot(specgrowthlen ~ ztemp, dat2, pch = NA, xlab = "Mean Temperature", ylab = "Growth Rate (% mass per day)")
# nd <- list(zlogdens = rep(0, times = nvals), ztemp = temps, zpriorlen = rep(0, times = nvals))
# preds <- predict(grmod, newdata = nd, interval = "confidence")
# polygon(c(temps, rev(temps)), c(preds[,2], rev(preds[,3])), col = scales::alpha("black", 0.2), lty=0)
# lines(preds[,1] ~ temps, col = "black", lwd = 3)
# 
# plot(specgrowthlen ~ zpriorlen, dat2, pch = NA, xlab = "Prior Length (mm)", ylab = "Growth Rate (% mass per day)")
# nd <- list(zlogdens = rep(0, times = nvals), ztemp = rep(0, times = nvals), zpriorlen = lens)
# preds <- predict(grmod, newdata = nd, interval = "confidence")
# polygon(c(lens, rev(lens)), c(preds[,2], rev(preds[,3])), col = scales::alpha("black", 0.2), lty=0)
# lines(preds[,1] ~ lens, col = "black", lwd = 3)
# 
# # any patterns in residuals related to groundwater input? 
# plot(grmod$residuals ~ dat2$gw)
# summary(lm(grmod$residuals ~ dat2$gw))
# # No. Indicates that variation in growth among streams and over time is driven by temperature and density. And that other factors, such as food quantity/quality may not vary with groundwater input.
# 
# # # growth ~ gw
# # plot(growthwtmeancensor ~ zgw, dat2, pch = NA, xlab = "Groundwater Influence", ylab = "Growth Rate (% mass per day)")
# # nd <- list(zlogdens = rep(0, times = nvals), ztemp = rep(0, times = nvals), zgw = gws)
# # preds <- predict(grmod, newdata = nd, interval = "confidence")
# # polygon(c(gws, rev(gws)), c(preds[,2], rev(preds[,3])), col = scales::alpha("black", 0.2), lty=0)
# # lines(preds[,1] ~ gws, col = "black", lwd = 3)
# 
# dev.off()
