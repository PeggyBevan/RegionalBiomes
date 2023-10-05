# Author: Peggy Bevan
# Date: 12/01/2021
# Title: 1. Creating site level biodiversity metrics in PREDICTS


# Packages ----------------------------------------------------------------
# install.packages(c("geosphere", "ape", "geiger", "wordcloud"))
# install.packages("Data/PredictsData/yarg_0.1-14.tar.gz", repos = NULL, type = "source")
# 
# install.packages("glmmADMB", repos = "http://R-Forge.R-project.org")
# 
# install.packages("Data/PredictsData/roquefort_0.1-2.tar.gz", repos = NULL, type = "source")


library(dplyr)  # %>% select() filter() bind_rows()
library(rgdal)
library(devtools)

install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions")

library(predictsFunctions)


# Data --------------------------------------------------------------------

#read in data
diversity<- readRDS("Data/PredictsData/database.rds")
dim(diversity)
# 3250404 67
# Code --------------------------------------------------------------------


# We need to edit some columns so they work with the PREDICTS functions (see tutorial)
  
diversity <- mutate(diversity,
                      Measurement = Effort_corrected_measurement,
                      Sampling_effort = Rescaled_sampling_effort)

diversity <- diversity[!diversity$Rank_of_study_common_taxon %in% 
                         c('Infraspecies','Species'),] 

# Correct effort-sensitive abundance measures (assumes linear relationship between effort and recorded abundance)
diversity <- predictsFunctions::CorrectSamplingEffort(diversity = diversity)

# merge any sites that are within the same land-use type 
# and that have identical coordinates, start and end dates.
# (this line takes ~5 mins to run)
diversity <- predictsFunctions::MergeSites(diversity, silent = TRUE)


# Site level diversity metrics - no taxa ----------------------------------
#calculate site level diversity metrics - no taxa
# this code takes > 10 mins to run
# sites <- diversity %>%
#   # add Diversity_metric_is_valid column
#   mutate(Diversity_metric_is_valid = TRUE) %>%
#   # calculate SiteMetrics  , including extra columns you want to keep
#   yarg::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Diversity_metric_unit", "Sampling_effort_unit", "Study_common_taxon", "Rank_of_study_common_taxon", "Country", "Ecoregion", "Biome", "Realm", "Wilderness_area", "Hotspot", "Taxon", "Phylum", "Class", "Km_to_nearest_edge_of_habitat", "Years_since_fragmentation_or_conversion")) %>%
#   # calculate the total abundance within each study
#   group_by(SS) %>%
#   mutate(MaxAbundance = ifelse(Diversity_metric_type == "Abundance",
#                                max(Total_abundance),
#                                NA)) %>%
#   ungroup() %>%
#   # now calculate the rescaled abundance (abundance divided by the maximum within each study)
#   mutate(RescaledAbundance = ifelse(
#     Diversity_metric_type == "Abundance",
#     Total_abundance/MaxAbundance,
#     NA))
# dim(sites)
# # 22678 40
# #The data contain 480 sources, 666 studies and 22678 sites
# 
# # save this data frame so you can come back to it #SAVED ON 13/01/2021 
# write.csv(sites, 'Data/02_PREDICTSDivMetrics.csv', row.names = F)
# 
# #again, but splitting by taxa before calculating site level diversity metrics. 
# 
# Including taxa ----------------------------------------------------------

#create my_taxa column with 3 factors: Invertebrate, Vertebrate and Plants
#remove fungi

diversity <- diversity %>%
  mutate(my_taxa = recode(Phylum, 
                          'Arthropoda' = 'Invertebrate',
                          'Tracheophyta'  = 'Plant',
                          'Chordata'  = 'Vertebrate',
                          'Ascomycota' = 'Fungi',
                          'Mollusca' = 'Invertebrate',
                          'Platyhelminthes' = 'Invertebrate',
                          'Annelida' = 'Invertebrate',
                          'Bryophyta' = 'Plant',
                          'Nematoda' = 'Invertebrate',
                          'Basidiomycota' = 'Fungi',
                          'Mycetozoa' = 'Fungi',
                          'Glomeromycota' = 'Fungi',
                          'Onychophora' = 'Invertebrate')
  )


table(diversity$my_taxa)

# In some cases, phylum is a blank cell so filling these in using Kingdom. 
# If Kingdom = Plantae, its a plant
# If = Fungi or Protozoa is Fungi
# If = Animalia, not sure. 

diversity <- diversity %>% mutate(my_taxa = ifelse(my_taxa == "", paste(Kingdom), paste(my_taxa)))

table(diversity$my_taxa)

# diversity <- diversity %>%
diversity <- diversity %>%
  mutate(my_taxa = recode(my_taxa,
                          "Plantae" = 'Plant'
  )
  ) 

#Unknown animalias: Hayward 2009 is all on vertebrates - 'unknown predator'
diversity$my_taxa[diversity$my_taxa == 'Animalia' & diversity$Reference == 'Hayward 2009'] <- 'Vertebrate'

# Gaigher & Samways 2010 - these are all 'Bdellodes sp.' which are a type of mite.
diversity$my_taxa[diversity$my_taxa == 'Animalia' & diversity$Reference == 'Gaigher and Samways 2010'] <- 'Invertebrate'

table(diversity$my_taxa)

# now i want to check if the list of sites for each taxon group is distinct or if they overlap. 
Plant_sites <- unique(diversity$SSS[diversity$my_taxa == 'Plant'])
Vert_sites <- unique(diversity$SSS[diversity$my_taxa == 'Vertebrate'])
Invert_sites <- unique(diversity$SSS[diversity$my_taxa == 'Invertebrate'])
Fungi_sites <- unique(diversity$SSS[diversity$my_taxa == 'Fungi'])

#number of sites for each group
length(Plant_sites) + length(Vert_sites) + length(Invert_sites) + length(Fungi_sites)
# 23187
length(unique(diversity$SSS))
# 22678

# this means that there are some sites overlapping multiple kingdoms, so we have to 
# separate them into different dataframes before merging sites.

plants <- diversity[diversity$my_taxa =='Plant',] %>% droplevels()
verts <- diversity[diversity$my_taxa =='Vertebrate',] %>% droplevels()
inverts <- diversity[diversity$my_taxa =='Invertebrate',] %>% droplevels()
fungi <- diversity[diversity$my_taxa =='Fungi',] %>% droplevels()

Plant_sites <- predictsFunctions::SiteMetrics(plants, extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Diversity_metric_unit", "Sampling_effort_unit", "Study_common_taxon", "Rank_of_study_common_taxon", "Country", "Ecoregion", "Biome", "Realm", "Wilderness_area", "Hotspot", 'my_taxa', "Km_to_nearest_edge_of_habitat", "Years_since_fragmentation_or_conversion"))


Vert_sites <- predictsFunctions::SiteMetrics(verts, extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Diversity_metric_unit", "Sampling_effort_unit", "Study_common_taxon", "Rank_of_study_common_taxon", "Country", "Ecoregion", "Biome", "Realm", "Wilderness_area", "Hotspot", 'my_taxa', "Km_to_nearest_edge_of_habitat", "Years_since_fragmentation_or_conversion"))


Invert_sites <- predictsFunctions::SiteMetrics(inverts, extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Diversity_metric_unit", "Sampling_effort_unit", "Study_common_taxon", "Rank_of_study_common_taxon", "Country", "Ecoregion", "Biome", "Realm", "Wilderness_area", "Hotspot", 'my_taxa', "Km_to_nearest_edge_of_habitat", "Years_since_fragmentation_or_conversion"))


Fungi_sites <- predictsFunctions::SiteMetrics(fungi, extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Diversity_metric_unit", "Sampling_effort_unit", "Study_common_taxon", "Rank_of_study_common_taxon", "Country", "Ecoregion", "Biome", "Realm", "Wilderness_area", "Hotspot", 'my_taxa', "Km_to_nearest_edge_of_habitat", "Years_since_fragmentation_or_conversion"))

sites_taxa <- rbind(Plant_sites, Vert_sites, Invert_sites, Fungi_sites)


# calculate the rescaled abundance within each study
sites_taxa <- sites_taxa %>% 
  group_by(SS, my_taxa) %>%
  mutate(MaxAbundance = ifelse(Diversity_metric_type == "Abundance",
                               max(Total_abundance),
                               NA)) %>%
  ungroup() %>%
  # now calculate the rescaled abundance (abundance divided by the maximum within each study)
  mutate(RescaledAbundance = ifelse(
    Diversity_metric_type == "Abundance",
    Total_abundance/MaxAbundance,
    NA))
dim(sites_taxa)
# 23187 34

write.csv(sites_taxa, 'Data/02_PREDICTSDivMetrics_taxa.csv', row.names = F)

