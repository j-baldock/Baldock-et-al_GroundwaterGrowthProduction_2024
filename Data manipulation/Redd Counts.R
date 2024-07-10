# REDD COUNTS

library(tidyverse)
library(lubridate)

dat <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Redd Counts/Redd counts_YOY growth sections_2022.csv")
dat["min"][is.na(dat["min"])] <- 0
dat["max"][is.na(dat["max"])] <- 0
dat$avg <- rowMeans(dat[,8:9], na.rm = T)
dat <- dat %>% group_by(stream, section, date) %>% summarize(min = sum(min), max = sum(max), avg = sum(avg))

write_csv(dat, "YOY_ReddCounts_Summary.csv")
dat <- read_csv("YOY_ReddCounts_Summary.csv")

# load section lengths and temp data...convert Amedian to groundwater input
seclen <- read_csv("YOY_CatchDensity_2021-2022.csv") %>% group_by(stream, section) %>% summarise(seclen = max(seclen)) %>% ungroup()
temp <- read_csv("YOY_Stream_TempMetrics.csv") %>% group_by(stream) %>% summarize(gw = 1-unique(Amedian)) %>% ungroup() %>% arrange(gw)
temp$cols <- (hcl.colors(dim(temp)[1], "Viridis"))
coldf <- temp %>% dplyr::select(stream, cols)
write_csv(coldf, "StreamColors.csv")

# simple plot showing groundwater input among streams
png("Figures/GroundwaterInput.png", units = "in", width = 6, height = 5.5, res = 1000)
par(mar = c(4,7,2,2))
barplot(temp$gw, names.arg = temp$stream, col = temp$cols, las = 1, xlab = "Groundwater Input", horiz = TRUE)
dev.off()

# Join data, prelimineary plotting
dat <- dat %>% left_join(temp) %>% left_join(seclen) %>% mutate(redddens = (avg/seclen)*100)
plot(avg ~ gw, dat, pch = 21, bg = dat$cols, xlab = "Groundwater Input", ylab = "Redd Count")
plot(redddens ~ gw, dat, pch = 21, bg = dat$cols, xlab = "Groundwater Input", ylab = "Redd Density (redds per 100m)")

# PLOT: logged redd density ~ groundwater input: higher spawning densities in spring-fed streams
png("Figures/ReddDensityxGW.png", units = "in", width = 6, height = 5.5, res = 1000)
plot(log(redddens+1) ~ jitter(gw, factor = 2), dat, pch = 21, bg = dat$cols, xlab = "Groundwater Input", ylab = "log(Redd Density), jittered", main = "2022 Redd Density",axes = F)
axis(1)
axis(2, at = log(c(1,5,10,25,50,100)), labels = c(0,5,10,25,50,100))
box(bty = "o")
legend("topleft", legend = temp$stream, pch = 21, pt.bg = temp$cols, cex = 0.8, bty = "n")
dev.off()
