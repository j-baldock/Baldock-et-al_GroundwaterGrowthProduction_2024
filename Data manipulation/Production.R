# PRODUCTION

library(tidyverse)
library(lubridate)



size <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/GrowthData_working_withAges_YOYonly.csv") %>%
  mutate(section = factor(section)) %>%
  mutate(stream = factor(stream)) %>%
  mutate(streamnum = as.numeric(stream)) %>%
  mutate(year = year(date)) %>%
  mutate(doy = yday(date)) %>%
  filter(!is.na(lengthmm) & !is.na(weightg))

density <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/YOY_CatchDensity_2021-2022.csv") # density
density$density100 <- density$density * 100

#-------------------------------------#
# Seasonal Production
#-------------------------------------#

streams <- unique(size$stream)
years <- unique(size$year)
nsims <- 1000 # number of bootstrap simulations

# matrices for production in grams and production in mm, respectively
prodmatwt <- array(data = NA, dim = c(nsims, 5, 4, length(streams), 2))
prodmatlen <- array(data = NA, dim = c(nsims, 5, 4, length(streams), 2))

for(m in 1:length(years)) {
  for(l in 1:length(streams)) {
    s <- size %>% filter(stream == streams[l] & year == years[m])
    dates <- unique(s$date)
    if(dim(s)[1] == 0) next
    # bootstrap appraoch to estimate daily production (w/ uncertainty) between sampling events, sensu Kaylor et al (2021)
    for(i in 2:length(dates)) {
      for(k in 1:4) {
        d0 <- s %>% filter(date == dates[i-1] & section == k)   # data from former date, section k
        if(dim(d0)[1] == 0) next
        d1 <- s %>% filter(date == dates[i] & section == k)     # data from latter date, section k
        if(dim(d1)[1] == 0) next
        for(j in 1:nsims) {
          w0 <- sample(d0$weightg, 1)                               # randomly sample weight from former date
          w1 <- sample(d1$weightg, 1)                               # randomly sample weight from latter date
          l0 <- sample(d0$lengthmm, 1)
          l1 <- sample(d1$lengthmm, 1)
          t0 <- unique(d0$doy)                                      # day of year of former sampling date
          t1 <- unique(d1$doy)                                      # day of year of latter sampling date
          dens <- density %>% filter(stream == streams[l] & date == dates[i] & section == k)
          dens1 <- unique(dens$density100)
          prodmatwt[j,i,k,l,m] <- ((w1 - w0) * (dens1)) / (t1 - t0)
          prodmatlen[j,i,k,l,m] <- ((l1 - l0) * (dens1)) / (t1 - t0)
        }
      }
    }
    # plot output
    par(mfrow = c(2,2), mar = c(3,3,1,1), oma = c(0,0,2,0))
    for(i in 1:4) {
      # boxplots of bootstrapped growth rates
      boxplot(prodmatwt[,,i,l,m], ylim = range(prodmatwt, na.rm = T))
      abline(h = 0, lty = 2)
      # points of bootstrapped growth rates
      for(j in 1:length(dates)) {
        points(prodmatwt[,j,i,l,m] ~ rep(j,nsims), pch = 1, col = scales::alpha("dodgerblue", 0.5))
      }
      # points of mean growth rates and standard deviations derived from bootstraps
      points(apply(prodmatwt[,,i,l,m], 2, mean) ~ c(1:5), pch = 16, col = "orange", cex = 1.5)
      arrows(c(1:5), apply(prodmatwt[,,i,l,m], 2, mean), c(1:5), apply(prodmatwt[,,i,l,m], 2, mean) + apply(prodmatwt[,,i,l,m], 2, sd), col = "orange", lwd = 2, length = 0)
      arrows(c(1:5), apply(prodmatwt[,,i,l,m], 2, mean), c(1:5), apply(prodmatwt[,,i,l,m], 2, mean) - apply(prodmatwt[,,i,l,m], 2, sd), col = "orange", lwd = 2, length = 0)
      # plots mean growth rates derived from changes in mean weight between sampling events
      #points(spgrome[i,] ~ c(1:5), pch = 16, col = "darkred", cex = 1.5)
      legend("topleft", legend = paste("S", i, sep = ""), bty = "n")
    }
    mtext(paste(streams[l], years[m], sep = " "), side = 3, outer = TRUE)
  }
}

# visualize output
par(mfrow = c(2,2), mar = c(3,3,1,1), oma = c(0,0,2,0))
for(m in 1:length(years)) {
  for(l in 1:length(streams)) {
    if(sum(is.na(prodmatwt[,,,l,m])) == 20000) next
    for(i in 1:4) {
      # boxplots of bootstrapped growth rates
      boxplot(prodmatwt[,,i,l,m], ylim = range(prodmatwt[,,,l,m], na.rm = T))
      abline(h = 0, lty = 2)
      # points of bootstrapped growth rates
      for(j in 1:length(dates)) {
        points(prodmatwt[,j,i,l,m] ~ rep(j,nsims), pch = 1, col = scales::alpha("dodgerblue", 0.5))
      }
      # points of mean growth rates and standard deviations derived from bootstraps
      points(apply(prodmatwt[,,i,l,m], 2, mean) ~ c(1:5), pch = 16, col = "orange", cex = 1.5)
      arrows(c(1:5), apply(prodmatwt[,,i,l,m], 2, mean), c(1:5), apply(prodmatwt[,,i,l,m], 2, mean) + apply(prodmatwt[,,i,l,m], 2, sd), col = "orange", lwd = 2, length = 0)
      arrows(c(1:5), apply(prodmatwt[,,i,l,m], 2, mean), c(1:5), apply(prodmatwt[,,i,l,m], 2, mean) - apply(prodmatwt[,,i,l,m], 2, sd), col = "orange", lwd = 2, length = 0)
      # plots mean growth rates derived from changes in mean weight between sampling events
      #points(spgrome[i,] ~ c(1:5), pch = 16, col = "darkred", cex = 1.5)
      legend("topleft", legend = paste("S", i, sep = ""), bty = "n")
    }
    mtext(paste(streams[l], years[m], sep = " "), side = 3, outer = TRUE)
  }
}

# summarize production (mean and sd) and save to file
prodsum <- data.frame(stream = character(), 
                       date = character(),
                       priordate = character(), 
                       reach = numeric(), 
                       year = numeric(), 
                       prodwtmean = double(), 
                       prodwtsd = double(),
                       prodwtq5 = double(),
                       prodwtq25 = double(),
                       prodwtq50 = double(),
                       prodwtq75 = double(),
                       prodwtq95 = double(),
                       prodlenmean = double(),
                       prodlensd = double(),
                       prodlenq5 = double(),
                       prodlenq25 = double(),
                       prodlenq50 = double(),
                       prodlenq75 = double(),
                       prodlenq95 = double(),
                       stringsAsFactors = FALSE)

for(m in 1:dim(prodmatwt)[5]) {       # years
  for(l in 1:dim(prodmatwt)[4]) {     # streams
    s <- size %>% filter(stream == streams[l] & year == years[m])
    dates <- unique(s$date)
    for(i in 1:dim(prodmatwt)[2]) {   # sample periods/dates
      for(j in 1:dim(prodmatwt)[3]) { # reaches
        priordate <- ifelse(length(dates[i-1]) == 0, NA, as.character(dates[i-1]))
        new <- data.frame(as.character(streams[l]), 
                          as.character(dates[i]), 
                          as.character(priordate), 
                          j, 
                          years[m], 
                          mean(prodmatwt[,i,j,l,m]), 
                          sd(prodmatwt[,i,j,l,m]),
                          quantile(prodmatwt[,i,j,l,m], probs = 0.05, na.rm = TRUE),
                          quantile(prodmatwt[,i,j,l,m], probs = 0.25, na.rm = TRUE),
                          quantile(prodmatwt[,i,j,l,m], probs = 0.50, na.rm = TRUE),
                          quantile(prodmatwt[,i,j,l,m], probs = 0.75, na.rm = TRUE),
                          quantile(prodmatwt[,i,j,l,m], probs = 0.95, na.rm = TRUE),
                          mean(prodmatlen[,i,j,l,m]), 
                          sd(prodmatlen[,i,j,l,m]),
                          quantile(prodmatlen[,i,j,l,m], probs = 0.05, na.rm = TRUE),
                          quantile(prodmatlen[,i,j,l,m], probs = 0.25, na.rm = TRUE),
                          quantile(prodmatlen[,i,j,l,m], probs = 0.50, na.rm = TRUE),
                          quantile(prodmatlen[,i,j,l,m], probs = 0.75, na.rm = TRUE),
                          quantile(prodmatlen[,i,j,l,m], probs = 0.95, na.rm = TRUE))
        prodsum[nrow(prodsum) + 1,] <- new
      }
    }
  }
}

write_csv(prodsum, "YOY_Production_Summary.csv")




###############################################
# code to look at temp, gw, and time as drivers of production
# Not sure how relevant/realistic this is given that we cannot account for fish outmigrating from tribs


prodmod <- lm(prodwtmean ~ ztemp + I(ztemp^2) + zgw, dat2)

newcols <- hcl.colors(2, "Geyser")
nvals <- 100
temps <- seq(from = min(dat2$ztemp, na.rm = T), to = max(dat2$ztemp, na.rm = T), length.out = nvals)
gws <- seq(from = min(dat2$zgw, na.rm = T), to = max(dat2$zgw, na.rm = T), length.out = nvals)
# production ~ temp + temp^2
plot(prodwtmean ~ ztemp, dat2, xlab = "Mean Temperature", ylab = "Production (g/100m/day)")
nd <- list(ztemp = temps, zgw =  rep(0, times = nvals))
preds <- predict(prodmod, newdata = nd, interval = "confidence")
polygon(c(temps, rev(temps)), c(preds[,2], rev(preds[,3])), col = scales::alpha(newcols[1], 0.2), lty=0)
lines(preds[,1] ~ temps, col = newcols[1], lwd = 3)
# growth ~ gw
plot(prodwtmean ~ zgw, dat2, pch = NA, ylim = c(0,2), xlab = "Groundwater Influence", ylab = "Production (g/100m/day)")
nd <- list(ztemp = rep(0, times = nvals), zgw = gws)
preds <- predict(prodmod, newdata = nd, interval = "confidence")
polygon(c(gws, rev(gws)), c(preds[,2], rev(preds[,3])), col = scales::alpha("black", 0.2), lty=0)
lines(preds[,1] ~ gws, col = "black", lwd = 3)


########################
# fit the linear modelpr
prodmod2 <- lm(log(prodwtmean) ~ doy + I(doy^2) + gw + doy*gw + I(doy^2)*gw, data = dat2)
summary(prodmod2)
plot(prodmod2)
summary(prodmod2)$adj.r.squared # model explains 47.7% percent of variation in growth rates...pretty good
sqrt(mean(prodmod2$residuals^2)) # RMSE of 1.7% body mass per day...honestly not too bad

# predict from the fitted model
gwvec <- seq(from = min(dat2$gw), to = max(dat2$gw), length.out = 100)
gpmat <- matrix(NA, nrow = max(dat2$doy, na.rm = T)-min(dat2$doy, na.rm = T)+1, ncol = length(gwvec))
for (i in 1:length(gwvec)) { gpmat[,i] <- predict(prodmod2, newdata = list(gw = rep(gwvec[i], times = dim(gpmat)[1]), doy = c(min(dat2$doy, na.rm = T):max(dat2$doy, na.rm = T)))) }

# another way to visualize model output
plot(log(prodwtmean) ~ doy, dat2, axes = F, xlab = "", ylab = "Production")
axis(1, at = c(213,244,274,305), labels = c("Aug","Sept","Oct","Nov"))
axis(2)
box(bty = "o")
doys <- c(min(dat2$doy, na.rm = T):max(dat2$doy, na.rm = T))
colsss <- c("blue", "purple", "red")
gwv <- c(min(gwvec), median(gwvec), max(gwvec))
for (i in 1:length(gwv)) {
  preds <- predict(prodmod2, newdata = list(gw = rep(gwv[i], times = length(doys)), doy = doys), interval = "confidence")
  polygon(c(doys, rev(doys)), c(preds[,2], rev(preds[,3])), col = scales::alpha(colsss[i], 0.2), lty=0)
  lines(preds[,1] ~ doys, col = colsss[i], lwd = 3)
}
legend("topright", legend = round(gwv, digits = 2), fill = colsss, bty = "n", title = "Groundwater \nInfluence")







