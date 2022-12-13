#02 - Prep model data
# Author: Peggy Bevan
# Date: 20/01/2021

# **Part1:**
#   
#   * What regional biomes are particularly underrepresented in the PREDICTS database?
#   
#   * Which regional biomes have enough data that they can be included in my models?

# Order of code:
# 1. Packages
# 2. Load data
# 3. Data prep
#  - some data manipulation of predicts (BD)
#  - refactoring variables
# 
# - SAVE data for model - 'Data/03_PREDICTSModelData.csv'

# 1. Packages ----------------------------------------------------------------

#library(rgdal)  # readOGR() spTransform()
library(dplyr) #mutate()
#library(kableExtra) #kable() (R markdown only)
library(knitr)
library(tidyr) #drop_na()
#library(lme4) #glmer()
#install.packages("../MetaAnalysis/Data/RawData/PredictsData/roquefort_0.1-2.tar.gz", repos = NULL, type = "source")
#library(roquefort)
#library(sjPlot)
#library(flextable)
#library(data.table)
library(reshape2)

# 2. Load Data --------------------------------------------------------------------

# PREDICTS: Load dataframe from script 01. Each row is a site or block of a study, with biodiversity metrics e.g. total abundance. 
BD<- read.csv("/Data/02_PREDICTSDivMetrics.csv", stringsAsFactors = T)
dim(BD)
# 22678 40

# 3. Prep Data --------------------------------------------------------------------

# create regional biome information
BD <- BD %>%
  mutate(
    Biome_num = recode_factor(Biome,
                              'Tropical & Subtropical Moist Broadleaf Forests' = 1,
                              'Temperate Broadleaf & Mixed Forests' = 4,
                              'Mediterranean Forests, Woodlands & Scrub' = 12,
                              'Tropical & Subtropical Grasslands, Savannas & Shrublands' = 7,
                              'Temperate Conifer Forests' = 5,
                              'Temperate Grasslands, Savannas & Shrublands' = 8,
                              'Deserts & Xeric Shrublands' = 13,
                              'Mangroves' = 14,
                              'Tropical & Subtropical Dry Broadleaf Forests' = 2,
                              'Montane Grasslands & Shrublands' = 10,
                              'Tropical & Subtropical Coniferous Forests' = 3,
                              'Boreal Forests/Taiga' = 6,
                              'Flooded Grasslands & Savannas' = 9,
                              'Tundra' = 11
    ),
    Realm_code = recode_factor(Realm,
                               'Indo-Malay' = 'IM',
                               'Palearctic' = 'PA',
                               'Australasia' = 'AA',
                               'Nearctic' = 'NA',
                               'Neotropic' = 'NT',
                               'Afrotropic' = 'AT',
                               'Oceania' = 'OC',
                               'Antarctica' = 'AN'
    )
  ) %>%
  #remove rows where land use is unknown.
  filter(Predominant_land_use!='Cannot decide') %>%
  # remove sites on urban land - not enough data
  filter(Predominant_land_use!='Urban') %>%
  droplevels()
# tnc stands for The Nature Conservancy
BD$RB_tnc = paste0(BD$Realm_code, BD$Biome_num)
unique(BD$RB_tnc)

# create land use variable
BD <- BD %>%
  dplyr::mutate(
    LandUse = recode(Predominant_land_use, #grouping secondary veg into one category
                     'Primary vegetation' = 
                       'Primary Vegetation',
                     'Intermediate secondary vegetation' = 
                       'Secondary Vegetation',
                     'Mature secondary vegetation' = 
                       'Secondary Vegetation',
                     'Secondary vegetation (indeterminate age)' = 
                       'Secondary Vegetation',
                     'Young secondary vegetation' = 
                       'Secondary Vegetation'
    ),
    LandUse = na_if(LandUse, "Cannot decide"),
    Use_intensity = na_if(Use_intensity, "Cannot decide")
  )
# put into order so Primary Veg is reference level
BD$LandUse <- factor(BD$LandUse, c('Primary Vegetation', 'Secondary Vegetation', 'Plantation forest', 'Cropland', 'Pasture')) %>%
  droplevels()


##create other land use combinations, and land use//intensity variables

BD <- BD %>%
  mutate(LandUse2 = recode_factor(LandUse, # binary landuses
                                  'Primary Vegetation' = 
                                    'Natural Vegetation',
                                  'Secondary Vegetation' = 
                                    'Natural Vegetation',
                                  'Cropland' = 'Agriculture',
                                  'Pasture' = 'Agriculture',
                                  'Plantation forest' = 'Agriculture'
  ),
  LandUse3 = recode_factor(LandUse, # PV, SV and Agriculture
                           'Primary Vegetation' = 'Primary Vegetation',
                           'Secondary Vegetation' = 'Secondary Vegetation',
                           'Cropland' = 'Agriculture',
                           'Pasture' = 'Agriculture',
                           'Plantation forest' = 'Agriculture'
  ),
  LandUse4 = recode_factor(LandUse, # PV, SV, Harvested and Pasture
                           'Primary Vegetation' = 'Primary Vegetation',
                           'Secondary Vegetation' = 'Secondary Vegetation',
                           'Plantation forest' = 'Harvested',
                           'Cropland' = 'Harvested',
                           'Pasture' = 'Pasture'
  ),
  LandUse5 = recode_factor(LandUse, # PV, SV, commercial forest & agriculture
                           'Primary Vegetation' = 'Primary Vegetation',
                           'Secondary Vegetation' = 'Secondary Vegetation',
                           'Plantation forest' = 'Plantation forest',
                           'Cropland' = 'Agriculture',
                           'Pasture' = 'Agriculture'
  ),
  #log response variables
  LogAbund = log(Total_abundance+1),
  LogRichness = log(Species_richness+1)
  ) %>%
  #class 'Light' use intensity as 'minimal'
  mutate(Use_intensity = recode_factor(Use_intensity,
                                         'Light use' = 'Minimal use')) %>%
  #combine land use and use intensity variables
  mutate(LU_UI = ifelse(is.na(Use_intensity), NA, paste0(LandUse, '_', Use_intensity)),
         LU_UI_2 = ifelse(is.na(Use_intensity), NA, paste0(LandUse2, '_', Use_intensity)),
         LU_UI_3 = ifelse(is.na(Use_intensity), NA, paste0(LandUse3, '_', Use_intensity)),
         LU_UI_4 = ifelse(is.na(Use_intensity), NA, paste0(LandUse4, '_', Use_intensity)),
         LU_UI_5 = ifelse(is.na(Use_intensity), NA, paste0(LandUse5, '_', Use_intensity))
         ) 


unique(BD$Use_intensity)

BD$LU_UI <- as.factor(BD$LU_UI)
BD$LU_UI_2 <- as.factor(BD$LU_UI_2)
BD$LU_UI_3 <- as.factor(BD$LU_UI_3)
BD$LU_UI_4 <- as.factor(BD$LU_UI_4)
BD$LU_UI_5 <- as.factor(BD$LU_UI_5)

#re-order levels to Primary minimal is reference level
BD$LU_UI <- factor(BD$LU_UI, c("Primary Vegetation_Minimal use", "Primary Vegetation_Intense use", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Plantation forest_Minimal use", "Plantation forest_Intense use","Cropland_Minimal use","Cropland_Intense use", "Pasture_Minimal use", "Pasture_Intense use")) %>%
  droplevels()
BD$LU_UI_2 <- factor(BD$LU_UI_2, c("Natural Vegetation_Minimal use", "Natural Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use")) %>%
  droplevels()
BD$LU_UI_3 <- factor(BD$LU_UI_3, c("Primary Vegetation_Minimal use", "Primary Vegetation_Intense use", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use")) %>%
  droplevels()
BD$LU_UI_4 <- factor(BD$LU_UI_4, c("Primary Vegetation_Minimal use", "Primary Vegetation_Intense use", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Harvested_Minimal use", "Harvested_Intense use", "Pasture_Minimal use", "Pasture_Intense use")) %>%
  droplevels()
BD$LU_UI_5 <- factor(BD$LU_UI_5, c("Primary Vegetation_Minimal use", "Primary Vegetation_Intense use", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use", "Plantation forest_Minimal use","Plantation forest_Intense use")) %>%
  droplevels()

# create taxa variables

#change blank cells to NA
BD$Class <- as.character(BD$Class)
BD$Class[BD$Class==""] <- NA
BD$Class <- as.factor(BD$Class)

BD <- BD %>%
  mutate(CommonTaxon_Phylum = recode(Phylum, 
                                     "Annelida" = 'Invertebrate',
                                     "Arthropoda" = 'Invertebrate',
                                     "Ascomycota" = 'Fungi',
                                     "Basidiomycota" = 'Fungi',
                                     "Bryophyta" = 'Plant',      
                                     "Chordata" = 'Vertebrate',
                                     "Glomeromycota" = 'Fungi',
                                     "Mollusca" = 'Invertebrate',
                                     "Nematoda" = 'Invertebrate', 
                                     "Platyhelminthes" = 'Invertebrate',
                                     "Tracheophyta" = 'Plant')
         ,
         CommonTaxon = recode(Class,
                              "Adenophorea" = 'Plant',
                              "Agaricomycetes" = 'Fungi',
                              "Amphibia" = 'Herptile',
                              "Arachnida" = 'Invertebrate',
                              "Aves" = "Bird",
                              "Bryopsida" = 'Plant',"Chilopoda" = 'Invertebrate',
                              "Clitellata" = 'Invertebrate',
                              "Entognatha" = 'Invertebrate',
                              "Eurotiomycetes" = 'Fungi',
                              "Gastropoda" = 'Invertebrate',
                              "Glomeromycetes" = "Fungi",
                              "Insecta" = "Invertebrate",
                              "Jungermanniopsida" = "Plant",
                              "Lecanoromycetes" = "Fungi",
                              "Liliopsida" = "Plant", 
                              "Magnoliopsida" = "Plant",
                              "Malacostraca"  = "Invertebrate",
                              "Mammalia" = "Mammal",         
                              "Pinopsida" = "Plant",
                              "Polypodiopsida" = "Plant",
                              "Reptilia" = "Herptile",
                              "Secernentea" = "Invertebrate")
  )
#In some cases, the highest common class could not be named. Using Coalesce to fill these
# empties based on my edited version of Phylum
BD$CommonTaxon <- coalesce(BD$CommonTaxon, BD$CommonTaxon_Phylum) %>%
  droplevels()

table(BD$CommonTaxon)


#For now, lets save the dataset. 
write.csv(BD, 'Data/03_PREDICTSModelData.csv', row.names = F)
#RDS is better as it saves level orders of factors
saveRDS(BD, 'Data/03_PREDICTSModelData.rds')
