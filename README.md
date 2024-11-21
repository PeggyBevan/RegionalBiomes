# RegionalBiomes

Analysing biodiversity responses to land use change across regional biomes

**How to use this repository**

The scripts are ordered according to the analysis process, and should be run in this order to make sure new data sets are created.

The entire analysis can be run by calling 'RunAllScriptsHere.R'. This should take approximately 2.5 hours depending on compute power.

**Package requirements:**

*The versions listed here are the versions used at publication, but future versions may be usable.*

-   R Version - 4.0.4

-   dplyr Version 1.1.3

-   sf Version 1.0-16

-   devtools Version 2.4.0

-   predictsFunctions Version 1.0

-   knitr Version 1.31

-   tidyr Version 1.3.0

-   reshape2 Version 1.4.4

-   flextable Version 0.6.4

-   officer Version 0.3.17

-   lme4 Version 1.1-28

-   sjPlot Version 2.8.7

-   data.table Version 1.14.2

-   ggplot2 Version 3.4.1

-   performance Version 0.7.0

-   StatisticalModels Version 0.1

-   cowplot Version 1.1.1

-   patchwork Version 1.1.1

### **Data Structure**

### **Data//**

Data must be downloaded from original sources and placed in this folder

#### PREDICTS data

This project uses data from the 2016 release of the PREDICTS database (Hudson et al., 2016). The dataset can be found at <https://doi.org/10.5519/0066354>

This analysis starts with the RDS format of the data, called 'database.rds'. Place within a folder called 'PredictsData'

#### **Terrestrial Ecoregions of the World Map**

Ecoregion data was originally downloaded from the The Nature Conservancy's (TNC) conservation atlas (<http://maps.tnc.org/gis_data.html>). At the time of publishing, the link to this dataset is broken , but the same dataset can be downloaded from ResourceWatch [here](https://resourcewatch.org/data/explore/bio021a-Terrestrial-Ecoregions?section=Discover&selectedCollection=&zoom=3&lat=0&lng=0&pitch=0&bearing=0&basemap=dark&labels=light&layers=%255B%257B%2522dataset%2522%253A%2522d0968f74-f5c1-40a1-b2b5-5bac5de5cb15%2522%252C%2522opacity%2522%253A1%252C%2522layer%2522%253A%252201152647-80b6-41fb-9ebc-48a5f2411327%2522%257D%255D&aoi=&page=1&sort=most-viewed&sortDirection=-1). Note, the labelling of ecoregions in this dataset follows the same format as the PREDICTS database. Other ecoregions maps are available, for example from [WWF](https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world), but will require some editing to fit with the version of the PREDICTS database used here. All ecoregion maps are based on Olson et al., 2001. 

### **Scripts//**

01_CreateSiteMetrics.R - uses database.rds from Data/PredictsData folder - this script uses the mergesites function to group observations from a study to one row to give species richness and total abundance - here we also choose which columns to retain. Output: 02_PREDICTSDivMetrics.csv

02_PreProcessingModelData.R - remove 'urban' and 'Cannot decide' land use types - edit land use & use intensity factor names & add regional biome names, create landuse:use intensity factor for all 5 land use types - create taxa variable Output: 03_PREDICTSModelData.csv ; 03_PREDICTSModelData.rds

03_ExploreModelData.R - use this script to create summary tables of data - number of observations in each biome and regional biome - number of observations in each land use type by biome and regional biome - number of observations in each taxon group by biome and regional biome - it would be useful to save these tables as CSVs that can be brought back in later to subset

04_GlobalModels&Figs.R - run global model of species richness/abundance with land use change - create figures for: biodiversity change with land use - biodiversity change with land use and regional biome - biodiversity change with land use and regional biome and taxa

05_BiomeCaseStudies.R - Subset predicts database to 3 biomes, run regional biome model on species richness and abundance and plot predictions in change in metric over land-use change. Explore the impact of evening out sample sizes. Explore the impact of adding taxon to the species richness model.

PredictGLMERfunction.R - a copy of the PredictGLMER function from the StatisticalModels Package, called in case there are issues loading this package.

### **Figs//**

-   Directory for saving figures produced from scripts (created during runtime)

### **Output//**

-   Directory for any dataframes saved from analysis - primarily model selection tables used in supplementary information. (created during runtime)
