

library(tidyverse)
library(lubridate)
library(RColorBrewer)

# setwd("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth")
getwd()

dat <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Growth/GrowthData_working_withAges_YOYonly.csv") %>%
  mutate(section = factor(section)) %>%
  mutate(stream = factor(stream)) %>%
  mutate(streamnum = as.numeric(stream)) %>%
  mutate(year = year(date)) %>%
  mutate(doy = yday(date))

streams <- unique(dat$stream)
years <- unique(dat$year)
nsims <- 10 # number of bootstrap simulations

tempfiledf <- data.frame(file = list.files("Bioenergetics/Temperature files"),
                               year = rep(years, times = 13),
                               stream = c("Blackrock", "Blackrock",
                                          "Blacktail", "Blacktail",
                                          "Cliff", "Cliff",
                                          "Crystal", "Crystal",
                                          "Fall", "Fall",
                                          "Fish", "Fish",
                                          "Flat", "Flat",
                                          "Lower Bar BC", "Lower Bar BC",
                                          "Mosquito", "Mosquito",
                                          "Spread", "Spread",
                                          "3 Channel", "3 Channel",
                                          "Upper Bar BC", "Upper Bar BC",
                                          "Willow", "Willow"))

# initialize data frame
bioedat <- data.frame(stream = character(),
                      date = Date(),
                      priordate = Date(), 
                      section = numeric(), 
                      year = numeric(), 
                      NRun = numeric(), 
                      Run_Name = character(),
                      Species_Num = numeric(),
                      Species_txt = character(),
                      Initial_W = double(), 
                      fit.to = character(),
                      fit.to.value = double(),
                      First_day = double(),
                      Last_day = double(),
                      N_indiv = numeric(),
                      Oxycal = numeric(),
                      calc.pop_mort = logical(),
                      calc.spawn	= logical(),
                      calc.nut	= logical(),
                      calc.contaminant	= logical(),
                      Init_pred_conc	= numeric(),
                      contam_eq = numeric(),
                      T_File	= character(),
                      Diet_File = character(),
                      Prey_E_File = character(),
                      Indigest_Prey_File	= character(),
                      Pred_E_File = numeric(),
                      Mort_File = numeric(),
                      Reprod_File = numeric(),
                      Contam_Prey_C_File = numeric(),
                      Contam_assim_File = numeric(),
                      Contam_TE_File = numeric(),
                      Phos_Ae_File = numeric(),
                      Phos_Co_Pred_File = numeric(),
                      Phos_Co_Prey_File = numeric(),
                      Nit_Ae_File = numeric(),
                      Nit_Co_Pred_File = numeric(),
                      Nit_Co_Prey_File = numeric(),
                      FB4_Log_File = character(),
                      Save_Daily_Output = logical(),
                      Daily_Output_File = numeric(),
                      stringsAsFactors = FALSE)

for(m in 1:length(years)) {
  for(l in 1:length(streams)) {
    dat2 <- dat %>% filter(stream == streams[l] & year == years[m] & !is.na(lengthmm) & !is.na(weightg))
    dates <- unique(dat2$date)
    if(dim(dat2)[1] == 0) next
    for(i in 2:length(dates)) {
      for(k in 1:4) {
        d0 <- dat2 %>% filter(date == dates[i-1] & section == k)   # data from former date, section k
        if(dim(d0)[1] == 0) next
        d1 <- dat2 %>% filter(date == dates[i] & section == k)     # data from latter date, section k
        if(dim(d1)[1] == 0) next
        for(j in 1:nsims) {
          new <-      data.frame(stream = as.character(streams[l]),
                                 date = dates[i],
                                 priordate = dates[i-1], 
                                 section = k, 
                                 year = years[m], 
                                 NRun = NA, 
                                 Run_Name = paste(streams[l], dates[i], k, years[m], j, sep = "_"),
                                 Species_Num = 23,
                                 Species_txt = "Cutthroat trout",
                                 Initial_W = sample(d0$weightg, 1), 
                                 fit.to = "Weight",
                                 fit.to.value = sample(d1$weightg, 1),
                                 First_day = unique(d0$doy),
                                 Last_day = unique(d1$doy),
                                 N_indiv = 1,
                                 Oxycal = 13560,
                                 calc.pop_mort = FALSE,
                                 calc.spawn	= FALSE,
                                 calc.nut	= FALSE,
                                 calc.contaminant	= FALSE,
                                 Init_pred_conc	= NA,
                                 contam_eq	= NA,
                                 T_File	= paste("Main Inputs/BaldockBioE/Temperature Files/", tempfiledf$file[which(tempfiledf$year == years[m] & tempfiledf$stream == streams[l])], sep = ""),
                                 Diet_File = "Main Inputs/BaldockBioE/Diet_prop_inverts.csv",
                                 Prey_E_File = "Main Inputs/BaldockBioE/Prey_E_inverts_4000.csv",
                                 Indigest_Prey_File	= "Main Inputs/BaldockBioE/Indigestible_Prey_inverts.csv",
                                 Pred_E_File = NA,
                                 Mort_File = NA,
                                 Reprod_File = NA,
                                 Contam_Prey_C_File = NA,
                                 Contam_assim_File = NA,
                                 Contam_TE_File = NA,
                                 Phos_Ae_File = NA,
                                 Phos_Co_Pred_File = NA,
                                 Phos_Co_Prey_File = NA,
                                 Nit_Ae_File = NA,
                                 Nit_Co_Pred_File = NA,
                                 Nit_Co_Prey_File = NA,
                                 FB4_Log_File = "FB4_Log_File_BaldockYOYGrowth.csv",
                                 Save_Daily_Output = FALSE,
                                 Daily_Output_File = NA,
                                 stringsAsFactors = FALSE)
          bioedat[nrow(bioedat) + 1,] <- new
        }
      }
    }
  }
}
bioedat$NRun <- 1:nrow(bioedat)
view(bioedat)

write_csv(bioedat, "Bioenergetics/YOYGrowth_BioE_DesignFile.csv")


#---------------------------------------#
# Using mean weight
#---------------------------------------#

# Using the bootstrap approach per Kaylor et al (2021) and file produced above. There is so much overlap in size distributions
# that many errors are thrown by the bioenergetics software when initial weight is much larger than final weight. Thus, use 
# change in mean weights per sampling event to estimate mean consumption.

# initialize data frame
bioedat2 <- data.frame(stream = character(),
                      date = Date(),
                      priordate = Date(), 
                      section = numeric(), 
                      year = numeric(), 
                      catch = numeric(),
                      priorcatch = numeric(),
                      NRun = numeric(), 
                      Run_Name = character(),
                      Species_Num = numeric(),
                      Species_txt = character(),
                      Initial_W = double(), 
                      fit.to = character(),
                      fit.to.value = double(),
                      First_day = double(),
                      Last_day = double(),
                      N_indiv = numeric(),
                      Oxycal = numeric(),
                      calc.pop_mort = logical(),
                      calc.spawn	= logical(),
                      calc.nut	= logical(),
                      calc.contaminant	= logical(),
                      Init_pred_conc	= numeric(),
                      contam_eq = numeric(),
                      T_File	= character(),
                      Diet_File = character(),
                      Prey_E_File = character(),
                      Indigest_Prey_File	= character(),
                      Pred_E_File = numeric(),
                      Mort_File = numeric(),
                      Reprod_File = numeric(),
                      Contam_Prey_C_File = numeric(),
                      Contam_assim_File = numeric(),
                      Contam_TE_File = numeric(),
                      Phos_Ae_File = numeric(),
                      Phos_Co_Pred_File = numeric(),
                      Phos_Co_Prey_File = numeric(),
                      Nit_Ae_File = numeric(),
                      Nit_Co_Pred_File = numeric(),
                      Nit_Co_Prey_File = numeric(),
                      FB4_Log_File = character(),
                      Save_Daily_Output = logical(),
                      Daily_Output_File = numeric(),
                      stringsAsFactors = FALSE)

for(m in 1:length(years)) {
  for(l in 1:length(streams)) {
    dat2 <- dat %>% filter(stream == streams[l] & year == years[m] & !is.na(lengthmm) & !is.na(weightg))
    dates <- unique(dat2$date)
    if(dim(dat2)[1] == 0) next
    for(i in 2:length(dates)) {
      for(k in 1:4) {
        d0 <- dat2 %>% filter(date == dates[i-1] & section == k)   # data from former date, section k
        if(dim(d0)[1] == 0) next
        d1 <- dat2 %>% filter(date == dates[i] & section == k)     # data from latter date, section k
        if(dim(d1)[1] == 0) next
        new <-      data.frame(stream = as.character(streams[l]),
                               date = dates[i],
                               priordate = dates[i-1], 
                               section = k, 
                               year = years[m], 
                               catch = nrow(d0),
                               priorcatch = nrow(d1),
                               NRun = NA, 
                               Run_Name = paste(streams[l], dates[i], k, years[m], sep = "_"),
                               Species_Num = 23,
                               Species_txt = "Cutthroat trout",
                               Initial_W = mean(d0$weightg, na.rm = T), 
                               fit.to = "Weight",
                               fit.to.value = mean(d1$weightg, na.rm = T),
                               First_day = unique(d0$doy),
                               Last_day = unique(d1$doy),
                               N_indiv = 1,
                               Oxycal = 13560,
                               calc.pop_mort = FALSE,
                               calc.spawn	= FALSE,
                               calc.nut	= FALSE,
                               calc.contaminant	= FALSE,
                               Init_pred_conc	= NA,
                               contam_eq	= NA,
                               T_File	= paste("Main Inputs/BaldockBioE/Temperature Files/", tempfiledf$file[which(tempfiledf$year == years[m] & tempfiledf$stream == streams[l])], sep = ""),
                               Diet_File = "Main Inputs/BaldockBioE/Diet_prop_inverts.csv",
                               Prey_E_File = "Main Inputs/BaldockBioE/Prey_E_inverts_4000.csv",
                               Indigest_Prey_File	= "Main Inputs/BaldockBioE/Indigestible_Prey_inverts.csv",
                               Pred_E_File = NA,
                               Mort_File = NA,
                               Reprod_File = NA,
                               Contam_Prey_C_File = NA,
                               Contam_assim_File = NA,
                               Contam_TE_File = NA,
                               Phos_Ae_File = NA,
                               Phos_Co_Pred_File = NA,
                               Phos_Co_Prey_File = NA,
                               Nit_Ae_File = NA,
                               Nit_Co_Pred_File = NA,
                               Nit_Co_Prey_File = NA,
                               FB4_Log_File = "FB4_Log_File_BaldockYOYGrowthMeans.csv",
                               Save_Daily_Output = FALSE,
                               Daily_Output_File = NA,
                               stringsAsFactors = FALSE)
        bioedat2[nrow(bioedat2) + 1,] <- new
      }
    }
  }
}
bioedat2$NRun <- 1:nrow(bioedat2)
view(bioedat2)
hist(bioedat2$fit.to.value - bioedat2$Initial_W)

write_csv(bioedat2, "Bioenergetics/YOYGrowth_BioE_DesignFile_Means.csv")
