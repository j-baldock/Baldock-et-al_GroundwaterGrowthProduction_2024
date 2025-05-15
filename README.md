# Baldock-et-al_GroundwaterGrowthProduction_2024
Repository for code and derived data and model products associated with the manuscript entitled "Groundwater structures fish growth and production across a riverscape". All original data collected for this study can be found in the following U.S. Geological Survey data release: Al-Chokhachy, R., J. R. Baldock, and A. Walters. 2024. Fish sampling data, water temperature data, and groundwater spring location data from the upper Snake River basin, WY, 2021-2023: U.S. Geological Survey data release, https://doi.org/10.5066/P1YJCBPU.

Authors: Jeffrey R. Baldock, Robert Al-Chokhachy, and Annika Walters



"Data manipulation": various scripts and data files used to create "Growth_DataTable_WithCovariates.csv", which is used in the various growth and production R scripts. This file is created within the "Data Table 2.R" script. Note that some of the files called within "Data Table 2.R" are located in other directories (see Density and Temperature folders).


"Density": scripts and files used to calculate fish density across space and time, and to evaluate the relationship between groundwater index and fish density.


"Groundwater modeling": various scripts and files used to generate our index of groundwater influence continuously across the Upper Snake River Watershed. All MaxEnt output is provided. 
  - "GroundwaterIndex_InteractiveMap.html": this is an interactive map of groundwater influence across the study area


"Growth": R scripts and model outputs for modeling of growth as a function of either temperature and density or time and groundwater. "GrowthModel_Letcher.txt" is the JAGS model used in both analyses.


"Network projections": script used to project body size and cumulative production across the study area


"Production": R scripts and model outputs for modeling of production as a function of time and groundwater. "ProductionModel.txt" is the JAGS model used in the analysis.


"Stream Network": this contains the stream network/flowline shapefile generated as in the manuscript


"Temperature": this contains scripts used for the temperature modelling objective. "tsmod2.txt" is the JAGS model used in the analysis.
