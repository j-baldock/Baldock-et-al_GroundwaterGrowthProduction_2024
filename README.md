# Baldock-et-al_GroundwaterGrowthProduction_2024
Data and code repository for manuscript entitled "Groundwater structures fish growth and production across a riverscape" authored by Jeffrey R. Baldock, Robert Al-Chokhachy, and Annika Walters. 

"Growth_data_working_withAges_YOYonly.csv": this is raw fish catch data. 

"Growth_DataTable_WithCovariates.csv": this is the data file used in the various growth and production R scripts. This file is created within the "Data Table.R" script, located within the "Data manipulation" folder. The raw fish data can be found in the "GrowthData_working_withAges_YOYonly.csv" file, located within the "Data manipulation" folder.  


"Data manipulation": various scripts and data files used to create "Growth_DataTable_WithCovariates.csv", which is used in the various growth and production R scripts. This file is created within the "Data Table 2.R" script. Note that some of the files called within "Data Table 2.R" are located in other directories (see Density and Temperature folders).


"Density": scripts and files used to calculate fish density across space and time, and to evaluate the relationship between groundwater index and fish density.


"Groundwater modeling": various scripts and files used to generate our index of groundwater influence continuously across the Upper Snake River Watershed. All MaxEnt output is provided. 
  - "GroundwaterIndex_InteractiveMap.html": this is an interactive map of groundwater influence across the study area


"Growth": R scripts and model outputs for modeling of growth as a function of either temperature and density or time and groundwater. "GrowthModel_Letcher.txt" is the JAGS model used in both analyses.


"Network projections": script used to project body size and cumulative production across the study area


"Production": R scripts and model outputs for modeling of production as a function of time and groundwater. "ProductionModel.txt" is the JAGS model used in the analysis.


"Stream Network": this contains the stream network/flowline shapefile generated as in the manuscript


"Temperature": this contains scripts and data used in the temperature model. "tsmod2.txt" is the JAGS model used in the analysis. "TemperatureData_Daily.csv" is the raw temperature data (daily mean, min, and max) used in all temperature analyses. 
