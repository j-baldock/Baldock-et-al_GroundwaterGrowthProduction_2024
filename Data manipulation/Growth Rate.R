

library(tidyverse)
library(lubridate)
library(RColorBrewer)
library(ggpubr)

setwd("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth")
getwd()

dat <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/GrowthData_working_withAges_YOYonly.csv") %>%
  mutate(section = factor(section)) %>%
  mutate(stream = factor(stream)) %>%
  mutate(streamnum = as.numeric(stream)) %>%
  mutate(year = year(date)) %>%
  mutate(doy = yday(date))


# # Violin plots by weight
# par(mfrow = c(2,2))
# for (i in 1:4) {
#   d <- subset(dat2, dat2$section == i)
#   boxplot(weightg ~ date, d)
#   stripchart(weightg ~ date, d, 
#              method = "jitter",
#              pch = 16, col = "dodgerblue", vertical = TRUE, add = TRUE)
# }

#-------------------------#
# End of season body size 
#-------------------------#
# could play around with end date (final vs 4th sampling event), given low/lack of catches in some streams in Oct/Nov

# weight
streams <- unique(dat$stream)
par(mfrow = c(4,4), mar = c(2.5,2.5,0.5,0.5))
for(l in 1:length(streams)) {
  dat2 <- dat %>% filter(stream == streams[l] & date > "2022-01-01" & !is.na(lengthmm) & !is.na(weightg))
  dates <- unique(dat2$date)
  dat2 <- dat2 %>% filter(date == dates[5])
  boxplot(weightg ~ section, dat2, ylim = c(0,6))
  abline(h = 1, col = "darkred")
  legend("topleft", legend = streams[l], bty = "n")
}

# length
par(mfrow = c(4,4), mar = c(2.5,2.5,0.5,0.5))
for(l in 1:length(streams)) {
  dat2 <- dat %>% filter(stream == streams[l] & date > "2022-01-01" & !is.na(lengthmm) & !is.na(weightg))
  dates <- unique(dat2$date)
  dat2 <- dat2 %>% filter(date == dates[5])
  boxplot(lengthmm ~ section, dat2, ylim = c(0,100))
  abline(h = 50, col = "darkred")
  legend("topleft", legend = streams[l], bty = "n")
}

#-------------------------#
# Total Summer Growth Rate
#-------------------------#
# late emergence in some streams/asynchronous emergence among streams complicates things
# no fish caught during first sampling event in some streams --> can't calculate total growth rate
# could calc growth between second and final sampling events, but that overlooks early august when the majority of growth may be achieved in some streams (spread, cliff, willow)

streams <- unique(dat$stream)
years <- unique(dat$year)
nsims <- 1000 # number of bootstrap simulations
spgrotot <- array(data = NA, dim = c(nsims, 4, length(streams), 2))

for(m in 1:length(years)) {
  for(l in 1:length(streams)) {
    dat2 <- dat %>% filter(stream == streams[l] & year == years[m] & !is.na(lengthmm) & !is.na(weightg))
    dates <- unique(dat2$date)
    
    # bootstrap appraoch to estimate specific growth rates (w/ uncertainty) across critical period, sensu Kaylor et al (2021)
    for(k in 1:4) {
      d0 <- dat2 %>% filter(date == min(dates), section == k)     # data from first (August) sampling event, section k
      if(dim(d0)[1] == 0) next
      d1 <- dat2 %>% filter(date == max(dates), section == k)     # data from last (Oct/Nov) samping event, section k
      if(dim(d1)[1] == 0) next
      for(j in 1:nsims) {
        x0 <- sample(d0$weightg, 1)                               # randomly sample weight from former date
        x1 <- sample(d1$weightg, 1)                               # randomly sample weight from latter date
        t0 <- unique(d0$doy)                                      # day of year of former sampling date
        t1 <- unique(d1$doy)                                      # day of year of latter sampling date
        spgrotot[j,k,l,m] <- ((log(x1) - log(x0)) / (t1 - t0)) * 100   # specific growth rate
      }
    }
  }
}

# summarize (mean and sd) and input to matrix
spgrototsum <- data.frame(stream = character(), 
                       date = character(),
                       priordate = character(), 
                       reach = numeric(), 
                       year = numeric(), 
                       growthwtmean = double(), 
                       growthwtsd = double(),
                       growthwtq5 = double(),
                       growthwtq25 = double(),
                       growthwtq50 = double(),
                       growthwtq75 = double(),
                       growthwtq95 = double(),
                       stringsAsFactors = FALSE)

for(m in 1:dim(spgrotot)[4]) {       # years
  for(l in 1:dim(spgrotot)[3]) {     # streams
    dat2 <- dat %>% filter(stream == streams[l] & year == years[m] & !is.na(lengthmm) & !is.na(weightg))
    dates <- unique(dat2$date)
    for(i in 1:dim(spgrotot)[2]) {   # reaches
      new <- data.frame(as.character(streams[l]), 
                        as.character(max(dates)), 
                        as.character(min(dates)), 
                        i, 
                        years[m], 
                        mean(spgrotot[,i,l,m]), 
                        sd(spgrotot[,i,l,m]),
                        quantile(spgrotot[,i,l,m], probs = 0.05, na.rm = TRUE),
                        quantile(spgrotot[,i,l,m], probs = 0.25, na.rm = TRUE),
                        quantile(spgrotot[,i,l,m], probs = 0.50, na.rm = TRUE),
                        quantile(spgrotot[,i,l,m], probs = 0.75, na.rm = TRUE),
                        quantile(spgrotot[,i,l,m], probs = 0.95, na.rm = TRUE))
      spgrototsum[nrow(spgrototsum) + 1,] <- new
    }
  }
}

# write growth summary file
write_csv(spgrototsum, "YOY_SpecificGrowthRate_SummerTotal_Summary.csv")
spgrototsum <- read_csv("YOY_SpecificGrowthRate_SummerTotal_Summary.csv")

# plot
par(mfrow = c(4,4), mar = c(2.5,2.5,0.5,0.5))
for(l in 1:length(streams)) {
  # boxplots of bootstrapped growth rates
  boxplot(spgrotot[,,l], ylim = range(spgrotot, na.rm = T))
  abline(h = 0, lty = 2) 
  abline(h = mean(spgrotot, na.rm = T), col = "darkred")
  # points of bootstrapped growth rates
  for(j in 1:4) {
    points(spgrotot[,j,l] ~ rep(j,nsims), pch = 1, col = scales::alpha("dodgerblue", 0.5))
  }
  # points of mean growth rates and standard deviations derived from bootstraps
  points(apply(spgrotot[,,l], 2, mean) ~ c(1:4), pch = 16, col = "orange", cex = 1.5)
  arrows(c(1:4), apply(spgrotot[,,l], 2, mean), c(1:4), apply(spgrotot[,,l], 2, mean) + apply(spgrotot[,,l], 2, sd), col = "orange", lwd = 2, length = 0)
  arrows(c(1:4), apply(spgrotot[,,l], 2, mean), c(1:4), apply(spgrotot[,,l], 2, mean) - apply(spgrotot[,,l], 2, sd), col = "orange", lwd = 2, length = 0)
  
  legend("topleft", legend = streams[l], bty = "n")
}


#---------------------------------------------#
# Seasonal growth rate
#---------------------------------------------#

streams <- unique(dat$stream)
years <- unique(dat$year)
nsims <- 1000 # number of bootstrap simulations

# WEIGHT
spgroweight <- array(data = NA, dim = c(nsims, 5, 4, length(streams), 2))
for(m in 1:length(years)) {
  for(l in 1:length(streams)) {
    dat2 <- dat %>% filter(stream == streams[l] & year == years[m] & !is.na(lengthmm) & !is.na(weightg))
    dates <- unique(dat2$date)
    if(dim(dat2)[1] == 0) next
    # bootstrap appraoch to estimate specific growth rates (w/ uncertainty) between sampling events, sensu Kaylor et al (2021)
    for(i in 2:length(dates)) {
      for(k in 1:4) {
        d0 <- dat2 %>% filter(date == dates[i-1] & section == k)   # data from former date, section k
        if(dim(d0)[1] == 0) next
        d1 <- dat2 %>% filter(date == dates[i] & section == k)     # data from latter date, section k
        if(dim(d1)[1] == 0) next
        for(j in 1:nsims) {
          x0 <- sample(d0$weightg, 1)                               # randomly sample weight from former date
          x1 <- sample(d1$weightg, 1)                               # randomly sample weight from latter date
          t0 <- unique(d0$doy)                                      # day of year of former sampling date
          t1 <- unique(d1$doy)                                      # day of year of latter sampling date
          spgroweight[j,i,k,l,m] <- ((log(x1) - log(x0)) / (t1 - t0)) * 100   # specific growth rate
        }
      }
    }
    # plot output
    par(mfrow = c(2,2), mar = c(3,3,1,1), oma = c(0,0,2,0))
    for(i in 1:4) {
      # boxplots of bootstrapped growth rates
      boxplot(spgroweight[,,i,l,m], ylim = range(spgroweight, na.rm = T))
      abline(h = 0, lty = 2)
      # points of bootstrapped growth rates
      for(j in 1:length(dates)) {
        points(spgroweight[,j,i,l,m] ~ rep(j,nsims), pch = 1, col = scales::alpha("dodgerblue", 0.5))
      }
      # points of mean growth rates and standard deviations derived from bootstraps
      points(apply(spgroweight[,,i,l,m], 2, mean) ~ c(1:5), pch = 16, col = "orange", cex = 1.5)
      arrows(c(1:5), apply(spgroweight[,,i,l,m], 2, mean), c(1:5), apply(spgroweight[,,i,l,m], 2, mean) + apply(spgroweight[,,i,l,m], 2, sd), col = "orange", lwd = 2, length = 0)
      arrows(c(1:5), apply(spgroweight[,,i,l,m], 2, mean), c(1:5), apply(spgroweight[,,i,l,m], 2, mean) - apply(spgroweight[,,i,l,m], 2, sd), col = "orange", lwd = 2, length = 0)
      # plots mean growth rates derived from changes in mean weight between sampling events
      #points(spgrome[i,] ~ c(1:5), pch = 16, col = "darkred", cex = 1.5)
      legend("topleft", legend = paste("S", i, sep = ""), bty = "n")
    }
    mtext(paste(streams[l], years[m], sep = " "), side = 3, outer = TRUE)
  }
}

# LENGTH
spgrolength <- array(data = NA, dim = c(nsims, 5, 4, length(streams), 2))
for(m in 1:length(years)) {
  for(l in 1:length(streams)) {
    dat2 <- dat %>% filter(stream == streams[l] & year == years[m] & !is.na(lengthmm) & !is.na(weightg))
    dates <- unique(dat2$date)
    if(dim(dat2)[1] == 0) next
    # bootstrap appraoch to estimate specific growth rates (w/ uncertainty) between sampling events, sensu Kaylor et al (2021)
    for(i in 2:length(dates)) {
      for(k in 1:4) {
        d0 <- dat2 %>% filter(date == dates[i-1], section == k)   # data from former date, section k
        if(dim(d0)[1] == 0) next
        d1 <- dat2 %>% filter(date == dates[i], section == k)     # data from latter date, section k
        if(dim(d1)[1] == 0) next
        for(j in 1:nsims) {
          x0 <- sample(d0$lengthmm, 1)                               # randomly sample weight from former date
          x1 <- sample(d1$lengthmm, 1)                               # randomly sample weight from latter date
          t0 <- unique(d0$doy)                                      # day of year of former sampling date
          t1 <- unique(d1$doy)                                      # day of year of latter sampling date
          spgrolength[j,i,k,l,m] <- ((log(x1) - log(x0)) / (t1 - t0)) * 100   # specific growth rate
        }
      }
    }
    # plot output
    par(mfrow = c(2,2), mar = c(3,3,1,1), oma = c(0,0,2,0))
    for(i in 1:4) {
      # boxplots of bootstrapped growth rates
      boxplot(spgrolength[,,i,l,m], ylim = range(spgrolength, na.rm = T))
      abline(h = 0, lty = 2)
      # points of bootstrapped growth rates
      for(j in 1:length(dates)) {
        points(spgrolength[,j,i,l,m] ~ rep(j,nsims), pch = 1, col = scales::alpha("dodgerblue", 0.5))
      }
      # points of mean growth rates and standard deviations derived from bootstraps
      points(apply(spgrolength[,,i,l,m], 2, mean) ~ c(1:5), pch = 16, col = "orange", cex = 1.5)
      arrows(c(1:5), apply(spgrolength[,,i,l,m], 2, mean), c(1:5), apply(spgrolength[,,i,l,m], 2, mean) + apply(spgrolength[,,i,l,m], 2, sd), col = "orange", lwd = 2, length = 0)
      arrows(c(1:5), apply(spgrolength[,,i,l,m], 2, mean), c(1:5), apply(spgrolength[,,i,l,m], 2, mean) - apply(spgrolength[,,i,l,m], 2, sd), col = "orange", lwd = 2, length = 0)
      # plots mean growth rates derived from changes in mean weight between sampling events
      legend("topleft", legend = paste("S", i, sep = ""), bty = "n")
    }
    mtext(paste(streams[l], years[m], sep = " "), side = 3, outer = TRUE)
  }
}

# summarize (mean and sd) and input to matrix
spgrosum <- data.frame(stream = character(), 
                       date = character(),
                       priordate = character(), 
                       reach = numeric(), 
                       year = numeric(), 
                       growthwtmean = double(), 
                       growthwtsd = double(),
                       growthwtq5 = double(),
                       growthwtq25 = double(),
                       growthwtq50 = double(),
                       growthwtq75 = double(),
                       growthwtq95 = double(),
                       growthlenmean = double(),
                       growthlensd = double(),
                       growthlenq5 = double(),
                       growthlenq25 = double(),
                       growthlenq50 = double(),
                       growthlenq75 = double(),
                       growthlenq95 = double(),
                       stringsAsFactors = FALSE)

for(m in 1:dim(spgroweight)[5]) {       # years
  for(l in 1:dim(spgroweight)[4]) {     # streams
    dat2 <- dat %>% filter(stream == streams[l] & year == years[m] & !is.na(lengthmm) & !is.na(weightg))
    dates <- unique(dat2$date)
    for(i in 1:dim(spgroweight)[2]) {   # sample periods/dates
      for(j in 1:dim(spgroweight)[3]) { # reaches
        priordate <- ifelse(length(dates[i-1]) == 0, NA, as.character(dates[i-1]))
        new <- data.frame(as.character(streams[l]), 
                          as.character(dates[i]), 
                          as.character(priordate), 
                          j, 
                          years[m], 
                          mean(spgroweight[,i,j,l,m]), 
                          sd(spgroweight[,i,j,l,m]),
                          quantile(spgroweight[,i,j,l,m], probs = 0.05, na.rm = TRUE),
                          quantile(spgroweight[,i,j,l,m], probs = 0.25, na.rm = TRUE),
                          quantile(spgroweight[,i,j,l,m], probs = 0.50, na.rm = TRUE),
                          quantile(spgroweight[,i,j,l,m], probs = 0.75, na.rm = TRUE),
                          quantile(spgroweight[,i,j,l,m], probs = 0.95, na.rm = TRUE),
                          mean(spgrolength[,i,j,l,m]), 
                          sd(spgrolength[,i,j,l,m]),
                          quantile(spgrolength[,i,j,l,m], probs = 0.05, na.rm = TRUE),
                          quantile(spgrolength[,i,j,l,m], probs = 0.25, na.rm = TRUE),
                          quantile(spgrolength[,i,j,l,m], probs = 0.50, na.rm = TRUE),
                          quantile(spgrolength[,i,j,l,m], probs = 0.75, na.rm = TRUE),
                          quantile(spgrolength[,i,j,l,m], probs = 0.95, na.rm = TRUE))
        spgrosum[nrow(spgrosum) + 1,] <- new
      }
    }
  }
}

# write growth summary file
write_csv(spgrosum, "YOY_SpecificGrowthRate_Summary.csv")
spgrosum <- read_csv("YOY_SpecificGrowthRate_Summary.csv")







#-------------------------#
# Proportional seasonal growth rate
#-------------------------#

#### fully bootstrap approach
#### does not appear to work well given overlapping size distributions between sampling events
#### means do not match medians, means do not sum to 1
streams <- unique(dat$stream)
nsims <- 100 # number of bootstrap simulations
prgro <- array(data = NA, dim = c(nsims, 5, 4, length(streams)))
# length/weight at emergence from fry trapping
dEmer <- read_csv("/Volumes/GoogleDrive/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Databases/Redd cap data 2019-2021.csv") %>%
  filter(date > "2021-01-01") %>% select(length, weight) %>% drop_na() %>% rename(lengthmm = length, weightg = weight)

for(l in 1:length(streams)) {
  dat2 <- dat %>% filter(stream == streams[l] & date > "2022-01-01" & !is.na(lengthmm) & !is.na(weightg))
  dates <- unique(dat2$date)
  
  # bootstrap approach to estimate proportion of total growth (weight) attained between each sampling event, bootstrapped to account for uncertainty
  for(i in 1:length(dates)) {
    for(k in 1:4) {
      dFin <- dat2 %>% filter(date == max(dates), section == k) # data from final sampling date
      if(dim(dFin)[1] == 0) next
      if(i == 1) d0 <- dEmer else d0 <- dat2 %>% filter(date == dates[i-1], section == k)   # data from former date, section k
      if(dim(d0)[1] == 0 & i == 2) d00 <- dEmer else d00 <- d0
      if(dim(d00)[1] == 0) next
      d1 <- dat2 %>% filter(date == dates[i], section == k)     # data from latter date, section k
      if(dim(d1)[1] == 0) next
      for(j in 1:nsims) {
        xEmer <- sample(dEmer$weightg, 1)                              # randomly sample weight from emergence data
        xFin <- sample(dFin$weightg, 1)                               # randomly sample weight from final date
        x0 <- sample(d00$weightg, 1)                                   # randomly sample weight from former sampling event
        x1 <- sample(d1$weightg, 1)                                   # randomly sample weight from latter sampling event
        prgro[j,i,k,l] <- (x1 - x0) / (xFin - xEmer)                  # calculate proportional growth 
      }
    }
  }
  # plot output
  par(mfrow = c(2,2), mar = c(3,3,1,1), oma = c(0,0,2,0))
  for(i in 1:4) {
    # boxplots of bootstrapped growth rates
    boxplot(prgro[,,i,l], ylim = c(-1,1))
    abline(h = 0, lty = 2)
    abline(h = 1, lty = 1, col = "red")
    # points of bootstrapped growth rates
    for(j in 1:length(dates)) {
      points(prgro[,j,i,l] ~ rep(j,nsims), pch = 1, col = scales::alpha("dodgerblue", 0.5))
    }
    # points of mean growth rates and standard deviations derived from bootstraps
    points(apply(prgro[,,i,l], 2, mean) ~ c(1:5), pch = 16, col = "orange", cex = 1.5)
    arrows(c(1:5), apply(prgro[,,i,l], 2, mean), c(1:5), apply(prgro[,,i,l], 2, mean) + apply(prgro[,,i,l], 2, sd), col = "orange", lwd = 2, length = 0)
    arrows(c(1:5), apply(prgro[,,i,l], 2, mean), c(1:5), apply(prgro[,,i,l], 2, mean) - apply(prgro[,,i,l], 2, sd), col = "orange", lwd = 2, length = 0)
    # plots mean growth rates derived from changes in mean weight between sampling events
    #points(spgrome[i,] ~ c(1:5), pch = 16, col = "darkred", cex = 1.5)
    legend("topleft", legend = paste("S", i, sep = ""), bty = "n")
    legend("topright", legend = sum(apply(prgro[,,i,l], 2, mean), na.rm = T), bty = "n")
  }
  mtext(streams[l], side = 3, outer = TRUE)
}

# means only: seasonal proportional growth rate based on changes in mean weight between sampling events
# works well (sum to 1), but no uncertainty
prgro <- array(data = NA, dim = c(4, 5, length(streams)))
for(l in 1:length(streams)) {
  dat2 <- dat %>% filter(stream == streams[l] & date > "2022-01-01" & !is.na(lengthmm) & !is.na(weightg))
  dates <- unique(dat2$date)
  
  # bootstrap approach to estimate proportion of total growth (weight) attained between each sampling event, bootstrapped to account for uncertainty
  for(i in 1:length(dates)) {
    for(k in 1:4) {
      dFin <- dat2 %>% filter(date == max(dates), section == k) # data from final sampling date
      if(dim(dFin)[1] == 0) next
      if(i == 1) d0 <- dEmer else d0 <- dat2 %>% filter(date == dates[i-1], section == k)   # data from former date, section k
      if(dim(d0)[1] == 0 & i == 2) d00 <- dEmer else d00 <- d0
      if(dim(d00)[1] == 0) next
      d1 <- dat2 %>% filter(date == dates[i], section == k)     # data from latter date, section k
      if(dim(d1)[1] == 0) next
      
      xEmer <- mean(dEmer$weightg, na.rm = T)                              # mean weight from emergence data
      xFin <- mean(dFin$weightg, na.rm = T)                               # mean weight from final date
      x0 <- mean(d00$weightg, na.rm = T)                                   # mean weight from former sampling event
      x1 <- mean(d1$weightg, na.rm = T)                                   # mean weight from latter sampling event
      prgro[k,i,l] <- (x1 - x0) / (xFin - xEmer)
    }
  }
  # plot output
  par(mfrow = c(2,2), mar = c(3,3,1,1), oma = c(0,0,2,0))
  for(k in 1:4) {
    # boxplots of bootstrapped growth rates
    plot(prgro[k,,l] ~ c(1:5), ylim = c(-1,1))
    abline(h = 0, lty = 2)
    abline(h = 1, lty = 1, col = "red")
    legend("topleft", legend = paste("S", k, sep = ""), bty = "n")
    legend("topright", legend = sum(prgro[k,,l], na.rm = T), bty = "n")
  }
  mtext(streams[l], side = 3, outer = TRUE)
}




# mean growth rates, sensu Kaylor et al (2021)
# this seems to underestimate growth rates for some sampling periods, relative to the mean of bootstrapped growth rates
# spgrome <- matrix(data = NA, nrow = 4, ncol = length(dates))
# for(i in 2:length(dates)) {
#   for(k in 1:4) {
#       d0 <- dat2 %>% filter(date == dates[i-1], section == k)   # data from former date, section k
#       d1 <- dat2 %>% filter(date == dates[i], section == k)     # data from latter date, section k
#       spgrome[k,i] <- ((log(mean(d1$weightg)) - log(mean(d0$weightg))) / (unique(d1$doy) - unique(d0$doy))) * 100   # specific growth rate
#   }
# }





t.test(spgro[,2,1], spgro[,3,1], alternative = "two.sided", var.equal = FALSE)


