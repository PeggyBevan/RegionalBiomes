# RegionalBiomes
Analysing biodiversity responses to land use change across regional biomes

How to use this repository
- If you want to go straight to running models, it is possible to start in script 3 and load 03_PREDICTSModelData.rds

Data Structure

Data//
PredictsData
- original database from predicts website 
terr-ecoregions-TNX
The TNC egoregion map. It has the same ecoregions as those in the PREDICTS database. I got it from here: http://maps.tnc.org/gis_data.html


02_PREDICTSDivMetrics.csv
  - created in script 01. 
  - baseline dataset to use - do not edit


03_PREDICTSModelData.csv
 - created in script 02
 - urban land use and 'cannot decide' land uses have been removed.
 - extra variables added - land use (1-5), land use:use intensity (1-5), log richness, log abundance, common taxa


Scripts//

01_CreateSiteMetrics.R
 - uses database.rds from rawdata folder in another directory
 - groups observations from a study to one row to give species richness and total abundance
 - choose which columns to retain
 Output: 02_PREDICTSDivMetrics.csv
 
02_PreProcessingModelData
 - remove 'urban' and 'Cannot decide' land use types
 - edit land use & use intensity factor names & add regional biome names, create landuse:use intensity factor for all 5 land use types
 - create taxa variable
 Output: 03_PREDICTSModelData.csv ; 03_PREDICTSModelData.rds
 
03_ExploreModelData
 - use this script to create summary tables of data
 - number of observations in each biome and regional biome
 - number of observations in each land use type by biome and regional biome
 - number of observations in each taxon group by biome and regional biome
 - it would be useful to save these tables as CSVs that can be brought back in later to subset
 
04_Global model and figs
 - run global model of species richness/abundance with land use change 
 - create figures for: biodiversity change with land use
 - biodiversity change with land use and regional biome
 - biodiversity change with land use and regional biome and taxa


05_Biome1
 - Subset biome 1 and do any data tidying (remove data-deficient areas)
 - run model and model selection against species richness and abundance
 - plot models
 

06_Biome4
- repeat 05 but with biome 4


07_Biome7
- repeat 05 but with biome 7


Figs//
- figures produced from scripts
- may still contain some archived figures
