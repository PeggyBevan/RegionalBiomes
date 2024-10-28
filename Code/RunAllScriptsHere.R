###RUN EVERYTHING
#times are calulated on a Macbook pro with 2.3 GHz Quad-Core Intel Core i7 and 16GB RAM


source('Code/01_CreateSiteMetrics.R') #approx 20 mins
source('Code/02_PreProcessingModelData.R') #< 1 min
source('Code/03_ExploreModelData.R') #<1 min
source('Code/04_GlobalModels&Figs.R') #2 mins
source('Code/05_Biome1_RunModels.R') #2 mins
source('Code/06_Biome4_RunModels.R') #4mins
source('Code/07_Biome7_RunModels.R')
