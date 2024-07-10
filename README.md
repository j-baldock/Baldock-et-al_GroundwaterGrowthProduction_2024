# Baldock-et-al_GroundwaterGrowthProduction_2024
Data and code repository for manuscript entitled "Groundwater structures fish growth and production across a riverscape" authored by Jeffrey R. Baldock, Robert Al-Chokhachy, and Annika Walters. 

"Growth_DataTable_WithCovariates.csv": this is the data file used in the various growth and production R scripts. This file is created within the "Data Table.R" script, located within the "Data manipulation" folder. This file (and the R script used to created it) contain a number of different calculations and fields that ultimately were not used in the analysis, but were generated in the process of analytical methods development and exploration. The raw fish data can be found in the "GrowthData_working_withAges_YOYonly.csv" file, located within the "Data manipulation" folder.  

"Data manipulation": various scripts and data files used to create "Growth_DataTable_WithCovariates.csv".

"Density": scripts and files used to calculate fish density across space and time, and to evaluate the relationship between groundwater index and fish density.

"Groundwater modeling": various scripts and files used to generate our index of groundwater influence continuously across the Upper Snake River Watershed. All MaxEnt output is provided. 
  - "GroundwaterIndex_InteractiveMap.html": this is an interactive map of groundwater influence across the study area

"Growth": R scripts and model outputs for modeling of growth as a function of either temperature and density or time and groundwater. "GrowthModel_Letcher.txt" is the JAGS model used in both analyses.

"Network projections": script used to project body size and cumulative production across the study area

"Production": R scripts and model outputs for modeling of production as a function of time and groundwater. "ProductionModel.txt" is the JAGS model used in the analysis.

"Stream Network": this contains the stream network/flowline shapefile generated as in the manuscript

"Temperature": this contains scripts and data used in the temperature model. "tsmod2.txt" is the JAGS model used in the analysis
