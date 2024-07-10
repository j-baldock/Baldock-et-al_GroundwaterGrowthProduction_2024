#==========================================================================================================================================#
# PROJECT NAME: Groundwater structures fish growth and production across a riverscape
# CREATOR: Jeffrey R. Baldock, WY Coop Unit, University of Wyoming, jbaldock@uwyo.edu 
# PURPOSE: Plots exploring the relationship between redd density and yoy density (supplemental figures)
# NOTES:
#==========================================================================================================================================#

library(tidyverse)
library(lubridate)

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

# join colors to data
dat <- dat %>% left_join(coldf)

# fill in blank section lengths
dat$seclen[is.na(dat$seclen)] <- 150

# calculate redd density...add very small value to allow log of 0 redd density
dat <- dat %>% mutate(redddens = (avgredds/seclen), logredddens = log(redddens+0.001))

# summarize by maximum yoy density...2022 only as redd counts were not conducted in 2021
dat2 <- dat %>% filter(year %in% c(2022)) %>% 
  group_by(stream, section, gw, cols, year, yrid, redddens, logredddens) %>% 
  summarize(maxden = max(density, na.rm = T)) %>% 
  mutate(maxden_log = log(maxden + 0.01), gw_log = log(gw))

# dat2 <- dat %>% group_by(stream, year, section, logredddens, cols) %>% summarise(maxyoydens = max(density, na.rm = T))


##---------------------------------------------------------------------------------------##
## Plots
##---------------------------------------------------------------------------------------##

png("Figures/YOYDensityxReddDensity.png", units = "in", width = 9, height = 9, res = 1000)
par(mfrow = c(2,2))

# PLOT: YOY density ~ log redd density, all sampling events
plot(density*100 ~ logredddens, dat, pch = 21, bg = dat$cols, xlab = "Redd Density (per 100m), log scale", ylab = "YOY Density (per 100m)", main = "All Sampling Events", axes = F, ylim = c(0,400))
axis(1, at = log(c(1,5,10,25,50,100)), labels = c(0,5,10,25,50,100))
axis(2)
box(bty = "o")
legend("topleft", legend = coldf$stream, pch = 21, pt.bg  = coldf$cols, cex = 0.8, bty = "n")

# PLOT: log YOY density ~ log redd density, all sampling events
plot(log((density*100)+1) ~ logredddens, dat, pch = 21, bg = dat$cols, xlab = "Redd Density (per 100m), log scale", ylab = "YOY Density (per 100m), log scale", main = "All Sampling Events", axes = F, ylim = c(0,6))
axis(1, at = log(c(1,5,10,25,50,100)), labels = c(0,5,10,25,50,100))
axis(2, at = log(c(1,2,6,11,26,51,101,201,401)), labels = c(1,2,6,11,26,51,101,201,401)-1)
box(bty = "o")

# PLOT: YOY density ~ log redd density, Event at maximum yoy density
plot(maxyoydens*100 ~ logredddens, dat2, pch = 21, bg = dat2$cols, xlab = "Redd Density (per 100m), log scale", ylab = "Maximum YOY Density (per 100m)", main = "Event at Max YOY Density", axes = F, ylim = c(0,400))
axis(1, at = log(c(1,5,10,25,50,100)), labels = c(0,5,10,25,50,100))
axis(2)
box(bty = "o")

# PLOT: log YOY density ~ log redd density, Event at maximum yoy density
plot(log((maxyoydens*100)+1) ~ logredddens, dat2, pch = 21, bg = dat2$cols, xlab = "Redd Density (per 100m), log scale", ylab = "Maximum YOY Density (per 100m), log scale", main = "Event at Max YOY Density", axes = F, ylim = c(0,6))
axis(1, at = log(c(1,5,10,25,50,100)), labels = c(0,5,10,25,50,100))
axis(2, at = log(c(1,2,6,11,26,51,101,201,401)), labels = c(1,2,6,11,26,51,101,201,401)-1)
box(bty = "o")

dev.off()


### PLOT: Figure S10: redd density (log) ~ groundwater & peak YOY density (log) ~ redd density (log)
cor.test(dat2$gw, dat2$logredddens)
cor.test(dat2$maxden_log, dat2$logredddens)
# set up gradient legend
gradpal <- tibble(gw = seq(from = 0, to = 1, length.out = 100)) %>%
  mutate(cols = vrPal(100)[as.numeric(cut(gw, breaks = 100))]) %>% 
  filter(gw >= min(dat2$gw) & gw <= max(dat2$gw))


jpeg("Figures/FigureS10_YoyReddDensityxGW.jpg", units = "in", width = 4.25, height = 7.75, res = 1000)
par(mfrow = c(2,1), mar = c(4,4,1,1), mgp = c(2.5,1,0))

plot(logredddens ~ gw, dat2, pch = 21, bg = dat2$cols, bty = "l", xlab = "Groundwater index", ylab = "ln(Redd density, redds/meter)")
# legend("topleft", legend = c(paste("(a) Pearson's r = ", round(cor(dat2$gw, dat2$logredddens), digits = 3))), bty = "n")
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.02, 0.97, labels = "(a)", adj = 0, cex = 1.2)
text(0.02, 0.88, labels = paste("Pearson's r = ", round(cor(dat2$gw, dat2$logredddens), digits = 3)), adj = 0)
par(usr = usr)

plot(maxden_log ~ logredddens, dat2, pch = 21, bg = dat2$cols, bty = "l", xlab = "ln(Redd density, redds/meter)", ylab = "ln(Peak YOY density, fish/meter)")
# legend("topleft", legend = c(paste("(b) Pearson's r = ", round(cor(dat2$maxden_log, dat2$logredddens), digits = 3))), bty = "n")
usr <- par("usr")
par(usr = c(0,1,0,1))
text(0.02, 0.97, labels = "(b)", adj = 0, cex = 1.2)
text(0.02, 0.88, labels = paste("Pearson's r = ", round(cor(dat2$maxden_log, dat2$logredddens), digits = 3)), adj = 0)
text(0.87, 0.40, labels = "Groundwater \nindex", cex = 0.85)
par(usr = usr)
gradientLegend(valRange = c(min(gradpal$gw), max(gradpal$gw)), color = gradpal$cols, length = 0.5, inside = TRUE, 
               side = 2, n.seg = 3, pos = c(0.80, 0.05, 0.85, 0.30), dec = 2)

dev.off()