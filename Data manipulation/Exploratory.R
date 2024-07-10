

install.packages("stringi")
library(stringi)

library(tidyverse)
library(lubridate)
library(RColorBrewer)
library(ggpubr)

getwd()

dat <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/GrowthData_working_withAges_YOYonly.csv") %>%
  mutate(section = factor(section)) %>%
  mutate(stream = factor(stream)) %>%
  mutate(streamnum = as.numeric(stream)) %>%
  mutate(year = year(date))

dd <- dat %>% mutate(doy = yday(date))

# summarize data into stream/date/reach means and standard deviations
lwsumm <- dat %>% group_by(stream, date, section) %>% 
  summarise(lenmean = mean(lengthmm, na.rm = T), lensd = sd(lengthmm, na.rm = T), lenmed = median(lengthmm, na.rm = T),
            wtmean = mean(weightg, na.rm = T), wtsd = sd(weightg, na.rm = T), wtmed = median(weightg, na.rm = T)) %>% ungroup()

dat %>% filter(date > "2022-01-01") %>%
  ggplot() + geom_boxplot(aes(group = date, x = date, y = weightg)) + facet_wrap(~ stream)

write_csv(lwsumm, "YOY_LengthWeight_Summary.csv")



# define color palettes
mycols12 <- brewer.pal(12, "Set3")
mycols4 <- brewer.pal(4, "Set2")

streams <- unique(dat$stream)

###### LENGTH
# Violin plots by length - by section
for (i in 1:length(streams)) {
  d <- subset(dat, dat$stream == streams[i])
  print(d %>% ggplot(aes(x = section, y = lengthmm)) + 
          geom_violin() + 
          geom_jitter(shape = 16, position = position_jitter(0.2)) +
          labs(title = streams[i]) +
          facet_wrap(~ date))
}
# Violin plots by length - sections pooled
for (i in 1:length(streams)) {
  d <- subset(dat, dat$stream == streams[i])
  print(d %>% ggplot(aes(x = date, y = lengthmm, group = date)) + 
          geom_violin() + 
          geom_jitter(shape = 16, position = position_jitter(1)) +
          labs(title = streams[i]) +
          facet_wrap(~ year, scales = "free_x"))
}

###### WEIGHT
# Violin plots by weight
for (i in 1:length(streams)) {
  d <- subset(dat, dat$stream == streams[i])
  print(d %>% ggplot(aes(x = section, y = weightg)) + 
          geom_violin() + 
          geom_jitter(shape = 16, position = position_jitter(0.2)) +
          labs(title = streams[i]) +
          facet_wrap(~ date))
}
# Violin plots by weight - sections pooled
for (i in 1:length(streams)) {
  d <- subset(dat, dat$stream == streams[i])
  print(d %>% ggplot(aes(x = date, y = weightg, group = date)) + 
          geom_violin() + 
          geom_jitter(shape = 16, position = position_jitter(1)) +
          labs(title = streams[i]) +
          facet_wrap(~ year, scales = "free_x"))
}

# boxplot
for (i in 1:length(streams)) {
  d <- subset(dat, dat$stream == streams[i])
  print(d %>% ggplot(aes(x = section, y = lengthmm)) + 
          geom_boxplot() + 
          labs(title = streams[i]) +
          facet_wrap(~ date))
}

###### ABUNDANCE
abundance <- dat %>% group_by(stream, date, year, section) %>% summarize(abun = n()) 

# for (i in 1:length(streams)) {
#   d <- subset(abundance, abundance$stream == streams[i])
#   print(d %>% ggplot(aes(x = date, y = abun, color = section)) +
#           geom_line() + geom_point() +
#           labs(title = streams[i]) +
#           facet_wrap(~ year, scales = "free_x"))
# }

abundance %>% ggplot(aes(x = date, y = abun, color = section)) +
  geom_line() + geom_point() +
  labs(title = streams[i]) +
  facet_wrap(~ year + stream, scales = "free")


###########################
# weight ~ length

# facet by stream
p <- ggscatter(dat, x = "lengthmm", y = "weightg", 
               xlab = "Length (mm)", ylab = "Weight (g)",
               shape = 21, color = "black", fill = "section", facet.by = "stream")
ggpar(p, legend = "none")

# color by stream
p <- ggscatter(dat, x = "lengthmm", y = "weightg", 
               xlab = "Length (mm)", ylab = "Weight (g)",
               shape = 21, color = "black", fill = "stream")
ggpar(p, legend = "right")

# condition factor by stream
datlw <- subset(dat, !is.na(lengthmm) & !is.na(weightg))
datlw$doy <- yday(datlw$date)
lwmod <- lm(log(weightg) ~ log(lengthmm), datlw)
summary(lwmod)
plot(lwmod)
datlw$resid <- lwmod$residuals
boxplot(resid ~ stream, datlw)
mm <- lme4::lmer(resid ~ (1|stream), datlw)
lme4::ranef(mm)

m1 <- lm(log(weightg) ~ log(lengthmm), datlw)
m2 <- lm(log(weightg) ~ log(lengthmm) + doy, datlw)
m3 <- lme4::lmer(log(weightg) ~ log(lengthmm) + (1|stream), datlw)
m4 <- lme4::lmer(log(weightg) ~ log(lengthmm) + (log(lengthmm)|stream), datlw)
m5 <- lme4::lmer(log(weightg) ~ log(lengthmm) + doy + (1|stream), datlw)
MuMIn::AICc(m1, m2, m3, m4, m5)

###########################
# length ~ time

p <- ggline(filter(dat, !is.na(lengthmm)), x = "date", y = "lengthmm", 
            add = "mean_sd", color = "section",
            xlab = "", ylab = "Length (mm)", facet.by = "stream", )
ggpar(p, legend = "right")


###########################
# weight ~ time

p <- ggline(filter(dat, !is.na(lengthmm)), x = "date", y = "weightg", 
            add = "mean_sd", color = "section", facet.by = "stream",
            xlab = "", ylab = "Weight (g)")
ggpar(p, legend = "right") 
set_palette(p, mycols4)

subAge0 <- subset(dat, dat$stream == "3 Channel" | dat$stream == "Fish" |
                    dat$stream == "Spread" | dat$stream == "Fall")
p <- ggline(filter(subAge0, !is.na(lengthmm)), x = "date", y = "weightg", 
            add = "mean_sd", color = "section", facet.by = "stream", nrow = 1,
            xlab = "", ylab = "Weight (g)")
ggpar(p, legend = "right") 
set_palette(p, mycols4)


###########################
# abundance ~ time

datAbunSum <- datAbun %>% group_by(stream, date, section, method, elap) %>% summarise(count = n())
datAbunSum$cpue <- datAbunSum$count / datAbunSum$elap

p <- ggline(datAbunSum, x = "date", y = "cpue", 
            color = "section", facet.by = "stream",
            xlab = "", ylab = "CPUE")
ggpar(p, legend = "right") 


###########################
# condition factor ~ time
