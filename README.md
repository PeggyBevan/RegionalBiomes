---
title: "Regional Biomes Framework"
author: "Anonymous"
output:
  html_document:
    df_print: paged
---

# RegionalBiomes

Analysing biodiversity responses to land use change across regional biomes

**How to use this repository**

The scripts are ordered according to the analysis process, and should be run in this order to make sure new data sets are created.

The entire analysis can be run by calling 'RunAllScriptsHere.R'. This should take approximately 2.5 hours depending on compute power.**Package requirements:**

*The versions listed here are the versions used at publication, but future versions may be useable.*

-   R Version - 4.0.4

-   dplyr Version 1.1.3

-   rgdal Version 1.5-23

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

PredictsData - original database from predicts website terr-ecoregions-TNX The TNC egoregion map. It has the same ecoregions as those in the PREDICTS database. I got it from here: <http://maps.tnc.org/gis_data.html>

02_PREDICTSDivMetrics.csv - created in script 01. - baseline dataset to use - do not edit

03_PREDICTSModelData.csv - created in script 02 - urban land use and 'cannot decide' land uses have been removed. - extra variables added - land use (1-5), land use:use intensity (1-5), log richness, log abundance, common taxa

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
