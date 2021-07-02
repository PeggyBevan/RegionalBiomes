# Author: Peggy Bevan
# Date: 12/01/2021
# Title: 1. Creating site level biodiversity metrics in PREDICTS


# Packages ----------------------------------------------------------------
install.packages("Data/RawData/PredictsData/yarg_0.1-14.tar.gz", repos = NULL, type = "source")

install.packages("glmmADMB", repos = "http://R-Forge.R-project.org")

install.packages("Data/RawData/PredictsData/roquefort_0.1-2.tar.gz", repos = NULL, type = "source")

library(dplyr)  # %>% select() filter() bind_rows()
library(rgdal)
library(yarg)
library(roquefort)


# Data --------------------------------------------------------------------

#read in data
diversity<- readRDS("Data/RawData/PredictsData/database.rds")
dim(diversity)
# 3250404 67
# Code --------------------------------------------------------------------


# We need to edit some columns so they work with the PREDICTS functions (see tutorial)
  
diversity <- mutate(diversity,
                      Measurement = Effort_corrected_measurement,
                      Sampling_effort = Rescaled_sampling_effort)

# merge any sites that are within the same land-use type 
# and that have identical coordinates, start and end dates.
# (this line takes ~5 mins to run)
diversity <- yarg::MergeSites(diversity, silent = TRUE, public = TRUE)


#calculate site level diversity metrics
# this code takes > 10 mins to run
sites <- diversity %>%
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  # calculate SiteMetrics  , including extra columns you want to keep
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Diversity_metric_unit", "Sampling_effort_unit", "Study_common_taxon", "Rank_of_study_common_taxon", "Country", "Ecoregion", "Biome", "Realm", "Wilderness_area", "Hotspot", "Taxon", "Phylum", "Class", "Km_to_nearest_edge_of_habitat", "Years_since_fragmentation_or_conversion")) %>%
  # calculate the total abundance within each study
  group_by(SS) %>%
  mutate(MaxAbundance = ifelse(Diversity_metric_type == "Abundance",
                               max(Total_abundance),
                               NA)) %>%
  ungroup() %>%
  # now calculate the rescaled abundance (abundance divided by the maximum within each study)
  mutate(RescaledAbundance = ifelse(
    Diversity_metric_type == "Abundance",
    Total_abundance/MaxAbundance,
    NA))
dim(sites)
# 22678 40
#The data contain 480 sources, 666 studies and 22678 sites

# save this data frame so you can come back to it #SAVED ON 13/01/2021 
write.csv(sites, 'FinalScriptsAndData/Data/02_PREDICTSDivMetrics.csv', row.names = F)


##see how it changes in a subset to see how what this function is doing to the numbers 

Poveda <- subset(diversity, SS == "SE1_2012__Poveda 2" %>%
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  # calculate SiteMetrics  , including extra columns you want to keep
  yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Diversity_metric_unit", "Sampling_effort_unit", "Study_common_taxon", "Rank_of_study_common_taxon", "Country", "Ecoregion", "Biome", "Realm", "Wilderness_area", "Hotspot", "Taxon", "Phylum", "Class", "Km_to_nearest_edge_of_habitat", "Years_since_fragmentation_or_conversion")) %>%
  # calculate the total abundance within each study
  group_by(SS) %>%
  mutate(MaxAbundance = ifelse(Diversity_metric_type == "Abundance",
                               max(Total_abundance),
                               NA)) %>%
  ungroup() %>%
  # now calculate the rescaled abundance (abundance divided by the maximum within each study)
  mutate(RescaledAbundance = ifelse(
    Diversity_metric_type == "Abundance",
    Total_abundance/MaxAbundance,
    NA))
dim(sites)


