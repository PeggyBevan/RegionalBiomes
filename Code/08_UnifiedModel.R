# 04B - GLobal Models and Figs - 3 biome case study

#Aim
# Run same model as in 04. but reduce dataset to three biomes
# even out sample sizes so that there is not a massive bias to one regional biome
#Order of code
# 1. Packages
# 2. Functions
#runmodels()
# 3. Data
#Load predicts database & regional biome summary
# 4. Script
#run model with no sample size adjustment
#run with sample size adjustment
#SPecies Richness
#Abundance
#Land-Use intensity (not presented in final results)
#Taxon-Regional biome (species richness)
#Plot figures for supplementary info.



# Packages ----------------------------------------------------------------
library(lme4) #glmer()
library(sjPlot)
library(data.table) #rbindlist()
library(ggplot2)
library(performance)
library(devtools)
library(dplyr)
library(flextable) #saving output
library(tidyr)
library(rgdal)  # readOGR()
# load PREDICTGLMER function 
source('Code/PredictGLMERfunction.R')

install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions")
library(predictsFunctions)

install_github("timnewbold/StatisticalModels")
library(StatisticalModels)

# Data --------------------------------------------------------------------

data <- readRDS('Data/03_PREDICTSModelData_taxa.rds')
dim(data)
ecoregs <- read.csv("Data/04_RBsummary_taxa.csv", stringsAsFactors = T)
names(data)
tnc_ecoregions <- readOGR(dsn = 'Data/terr-ecoregions-TNC', layer = 'tnc_terr_ecoregions')
# Code --------------------------------------------------------------------
#prep data
#subset to the three most data rich biomes
unique(data$Biome)
Biomes = c('Temperate Broadleaf & Mixed Forests', 'Tropical & Subtropical Moist Broadleaf Forests','Tropical & Subtropical Grasslands, Savannas & Shrublands')
data <- subset(data, Biome %in% Biomes)
data <- droplevels(data)

# Selecting regional biomes ----------------------------------------------------

#using 25 site threshold. 
#regional biomes that have >25 observations in PV and Ag
RB_LUth.25 <- ecoregs[ecoregs$PV.th25 == 1 & ecoregs$Ag.th25 == 1,]

#to create a subset in data:
data_LUth.25 <- filter(data, RB_tnc %in% RB_LUth.25$RB_tnc)

ecoregs <- ecoregs %>%
  subset(RB_tnc %in% unique(data_LUth.25$RB_tnc))


# Model - no sample size adjustment ---------------------------------------

#the best fitting model from previous experiment
m1 = StatisticalModels::GLMER(modelData = data_LUth.25, responseVar = 'LogRichness',
                         fitFamily = 'gaussian', fixedStruct = "LandUse*RB_tnc", 
                         randomStruct = "(1|SS)+(1|SSB)", REML = T)
summary(m1$model)

data_LUth.25$RB_tnc = as.factor(data_LUth.25$RB_tnc)
RBs <- data.frame(RB = unique(data_LUth.25$RB_tnc))

test <- sapply(RBs, FUN = function(rb) {
  cat(paste0(rb,'\n'))
})
#predict species richness change
preds <- apply(RBs, 1, FUN = function(rb) {
  nd <- data.frame(LandUse=factor(
    c('Primary Vegetation','Secondary Vegetation',"Plantation forest", "Cropland", "Pasture"),
    levels=levels(data_LUth.25$LandUse)))
  nd$RB_tnc <- factor(rb, levels = levels(data_LUth.25$RB_tnc))
  nd$LogRichness = 0
  preds <- StatisticalModels::PredictGLMER(model = m1$model, data = nd,
                                           se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(RB = nd$RB_tnc, LU=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
} )

preds2 <- data.table::rbindlist(preds)  

#include sample size
for (i in 1:nrow(preds2)) {
  preds2$n[i] <- nrow(subset(data_LUth.25, RB_tnc == preds2$RB[i] & LandUse == preds2$LU[i])) 
}

preds2$gt10 = ifelse(test = preds2$n > 10, yes = 1, no = 0)
preds2$gt10 <- as.factor(preds2$gt10)

preds2 <- preds2 %>%
  mutate(Realm_code = substr(RB, 1,2),
         Biome_num = substr(RB, 3,3),
         Realm = recode_factor(Realm_code,
                               'IM' = 'Indo-Malay',
                               'PA' = 'Palearctic',
                               'AA' = 'Australasia',
                               'NA' = 'Nearctic',
                               'NT' = 'Neotropic',
                               'AT' = 'Afrotropic'),
         Biome = recode_factor(Biome_num,
                               "1" = 'Tropical & Subtropical Moist Broadleaf Forests',
                               "4" = 'Temperate Broadleaf & Mixed Forests',
                               "7" = 'Tropical & Subtropical Grasslands, Savannas & Shrublands'
                               ))

# preds2$Biome_num <- as.numeric(preds2$Biome_num) 
preds2<-arrange(preds2, Biome_num, Realm)
preds2= subset(preds2, n > 20)
preds2$Realm = factor(preds2$Realm, levels = c('Afrotropic', 'Australasia','Indo-Malay', 'Nearctic', 'Neotropic', 'Palearctic'))

#plot richness predictions
figb <- ggplot(preds2, aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, group = Realm, colour = Realm)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 2) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#FC8D62","#A6D854","#8DA0CB", "#FFD92F")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6)) +
  #geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), linewidth = 1) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', "SV", "PF", "Pa", "Cr")) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100), hjust = 1.1, angle = 45,  data = preds2[preds2$upper > 100,], colour = 'black', position = position_dodge(width = 0.6), size = 2) +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    text = element_text(size = 9, colour = 'black'),
    legend.position = 'top',
    legend.title = element_blank(),
    strip.text = element_text(hjust = 0)) +
  facet_wrap(~Biome, nrow = 3)


figb

# Model - sample size adjustment ------------------------------------------

#need to combine study and taxa
data_LUth.25$Study_name = gsub("-", "_", data_LUth.25$Study_name)
data_LUth.25$Study_name_taxa = paste0(data_LUth.25$Study_name, "-", data_LUth.25$my_taxa)

#the average number of sites is
# mean(ecoregs$n_sites)
# #1293
# median(ecoregs$n_sites)
# #856
# range(ecoregs$n_sites)
#178 - 3941
#the average number of studies is also important because these represent different locations within a regional biome and represents variety more than sites - one study could have 100 sites.
# mean(ecoregs$n_studies)
# #42
# median(ecoregs$n_studies)
# #24
# range(ecoregs$n_studies)
# #13 - 127
# hist(data_LUth.25$Species_richness)

# Sample size adjustment ---------------------------------

#having looked at the sample sizes, there is actually not that much discrepancy between biomes, as the reviewer suggested. 
#in fact, there are two regional biomes which are much higher than the others.
mean(ecoregs$n_studies)
median(ecoregs$n_studies)
#there is a big difference between mean and median, showing a skew to ridiculously high numbers.
#if i halve the two extreme regional biomes, how does mean/median change?

#remove fungi 
data_LUth.25 = subset(data_LUth.25, my_taxa!= 'Fungi')
#edit ecoregs n_taxa so it doesn't include fungi. 

set.seed(9)
#create empty dataframe
data_samp_max = NULL
for (i in 1:nrow(ecoregs)) {
  rb = as.character(ecoregs$RB_tnc[i])
  if(ecoregs$n_studies[i] > 100) {
    #get list of studies in that rb
    rbstudies = unique(data_LUth.25$Study_name_taxa[data_LUth.25$RB_tnc==rb])
    n_taxa = n_distinct(data_LUth.25$my_taxa[data_LUth.25$RB_tnc==rb])
    studies = ecoregs$n_studies[i]
    newstudies = studies/2
    studiessample = sample(rbstudies, newstudies)
    samp = data_LUth.25 %>%
      subset(RB_tnc == rb) %>%
      subset(Study_name_taxa %in% studiessample)
    n_taxa_s = n_distinct(samp$my_taxa)
    if (n_taxa_s < n_taxa) {
      print(paste0(rb, ' lost a taxon'))
    }
    data_samp_max = rbind(data_samp_max, samp)
  } else {
    samp =  data_LUth.25 %>%
      subset(RB_tnc == rb)
    data_samp_max = rbind(data_samp_max, samp)
  } #end else
} #end for loop

#make new summary table
ecoregs_max <- data.frame('RB_tnc' = unique(data_samp_max$RB_tnc)) %>%
  mutate(Realm_code = substr(RB_tnc, 1,2),
         Biome_num = substr(RB_tnc, 3,3),
         Realm = recode_factor(Realm_code,
                               'IM' = 'Indo-Malay',
                               'PA' = 'Palearctic',
                               'AA' = 'Australasia',
                               'NA' = 'Nearctic',
                               'NT' = 'Neotropic',
                               'AT' = 'Afrotropic'
         ),
         Biome = recode_factor(Biome_num,
                               "1" = 'Tropical & Subtropical Moist Broadleaf Forests',
                               "4" = 'Temperate Broadleaf & Mixed Forests',
                               "7" = 'Tropical & Subtropical Grasslands, Savannas & Shrublands')
  )

ecoregs_max$Biome_num<-as.numeric(ecoregs_max$Biome_num)
ecoregs_max<-arrange(ecoregs_max, Biome_num, Realm)
# calculate number of studies etc. for each RB
for (i in (1:nrow(ecoregs_max))) {
  rb = ecoregs_max$RB_tnc[i]
  tnc_code = paste(ecoregs_max$RB_tnc[i])
  ecoregs_max$n_studies[i] = n_distinct(data_samp_max$Study_name_taxa[data_samp_max$RB_tnc==rb])
  ecoregs_max$n_sources[i] = n_distinct(data_samp_max$Source_ID[data_samp_max$RB_tnc==rb])
  ecoregs_max$n_sites[i] = n_distinct(data_samp_max$SSBS[data_samp_max$RB_tnc==rb]) 
  ecoregs_max$n_ecoregions[i] = n_distinct(data_samp_max$Ecoregion[data_samp_max$RB_tnc==rb])
  ecoregs_max$n_taxa[i] = n_distinct(data_samp_max$my_taxa[data_samp_max$RB_tnc==rb])
  ecoregs_max$n_country[i] =n_distinct(data_samp_max$Country[data_samp_max$RB_tnc==rb])
  ecoregs_max$t_ecoregions[i] = n_distinct(tnc_ecoregions@data$ECO_NAME[tnc_ecoregions@data$RealmMHT==tnc_code])
  ecoregs_max$p_ecoregions[i] = round(ecoregs_max$n_ecoregions[i]/ecoregs_max$t_ecoregions[i], 2)
}

rm(tnc_ecoregions)
LU <- data_samp_max %>% count(RB_tnc, LandUse) %>% 
  spread(LandUse, n)
ecoregs_max<- merge(ecoregs_max, LU, by.x = 'RB_tnc', all.x = TRUE)
#create agriculture count
ecoregs_max$Agriculture <- rowSums(ecoregs_max[,14:16], na.rm = T)

Taxa = data_samp_max %>% count(RB_tnc, my_taxa) %>%
  spread(my_taxa, n)

ecoregs_max = merge(ecoregs_max, Taxa, by.x = 'RB_tnc', all.x = TRUE)


ecoregs_max$`Primary Vegetation`
ecoregs_max$Agriculture


#now the mean is much closer to the median
mean(ecoregs_max$n_studies)
median(ecoregs_max$n_studies)
range(ecoregs_max$n_sites)
mean(ecoregs_max$n_sites)
median(ecoregs_max$n_sites)
# #there is now a less than factor of 10 change
range(ecoregs_max$n_studies)

#no taxonomic groups were lost from a regional biome during downsampling

## Species richness model --------------------------------------------------



m2 = StatisticalModels::GLMER(modelData = data_samp_max, responseVar = 'LogRichness',
                              fitFamily = 'gaussian', fixedStruct = "LandUse*RB_tnc", 
                              randomStruct = "(1|SS)+(1|SSB)", REML = F)
summary(m2$model)
data_samp_max$RB_tnc = as.factor(data_samp_max$RB_tnc)
RBs <- data.frame(RB = unique(data_samp_max$RB_tnc))

preds_max <- apply(RBs, 1, FUN = function(rb) {
  nd <- data.frame(LandUse=factor(
    c('Primary Vegetation','Secondary Vegetation',"Plantation forest", "Cropland", "Pasture"),
    levels=levels(data_samp_min$LandUse)))
  nd$RB_tnc <- factor(rb, levels = levels(data_samp_min$RB_tnc))
  nd$LogRichness = 0
  preds <- StatisticalModels::PredictGLMER(model = m2$model, data = nd,
                                           se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  #75% CIs
  preds75 <- StatisticalModels::PredictGLMER(model = m2$model, data = nd,
                                             se.fit = TRUE, seMultiplier = 1.15, randEffs = F)
  preds75$yplus75 <- ((exp(preds75$yplus)/exp(preds75$y[1]))*100)-100
  preds75$yminus75 <- ((exp(preds75$yminus)/exp(preds75$y[1]))*100)-100
  preds75$y75 <- ((exp(preds75$y)/exp(preds75$y[1]))*100)-100
  
  temp.df <- data.frame(RB = nd$RB_tnc, LU=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus, y75 = preds75$y75, lower75=preds75$yminus75, upper75=preds75$yplus75)
} )

preds2_max <- data.table::rbindlist(preds_max)  

#include sample size
for (i in 1:nrow(preds2_max)) {
  preds2_max$n[i] <- nrow(subset(data_samp_max, RB_tnc == preds2_max$RB[i] & LandUse == preds2_max$LU[i])) 
}

preds2_max <- preds2_max %>%
  mutate(Realm_code = substr(RB, 1,2),
         Biome_num = substr(RB, 3,3),
         Realm = recode_factor(Realm_code,
                               'IM' = 'Indo-Malay',
                               'PA' = 'Palearctic',
                               'AA' = 'Australasia',
                               'NA' = 'Nearctic',
                               'NT' = 'Neotropic',
                               'AT' = 'Afrotropic'),
         Biome = recode_factor(Biome_num,
                               "1" = 'a) Tropical Forest',
                               "4" = 'b) Temperate Forest',
                               "7" = 'c) Tropical Grassland'
         ))

# preds2$Biome_num <- as.numeric(preds2$Biome_num) 
preds2_max<-arrange(preds2_max, Biome_num, Realm)
#remove ranks with less than 25 detections
preds2_max= subset(preds2_max, n > 25)
preds2_max$Realm = factor(preds2_max$Realm, levels = c('Afrotropic', 'Australasia','Indo-Malay', 'Nearctic', 'Neotropic', 'Palearctic'))

level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")

## Species richness plot ---------------------------------------------------


figx <- ggplot(preds2_max, aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, group = Realm, colour = Realm)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 2) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#FC8D62","#A6D854","#8DA0CB", "#FFD92F")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6)) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), linewidth = 1) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', "SV", "PF", "Pa", "Cr")) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100), hjust = 1.1, angle = 45,  data = preds2_max[preds2_max$upper > 100,], colour = 'black', position = position_dodge(width = 0.6), size = 2) +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    text = element_text(size = 9, colour = 'black'),
    legend.position = 'top',
    legend.title = element_blank(),
    strip.text = element_text(hjust = 0)) +
  facet_wrap(~Biome, nrow = 3)

figx

ggsave('output/Eco_response/sampledmax_speciesrichness_casestudy.pdf', width = 7.5, height = 14, units = 'cm')


# Abundance model ---------------------------------------------------------

#remove all values over 10,000 
data_samp_max$LogAbund[data_samp_max$Total_abundance >= 10000] = NA
hist(data_samp_max$LogAbund)

ma1 = StatisticalModels::GLMER(modelData = data_samp_max, responseVar = 'LogAbund',
                               fitFamily = 'gaussian', fixedStruct = "LandUse*RB_tnc", 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
summary(ma1$model)

preds_mxa <- apply(RBs, 1, FUN = function(rb) {
  nd <- data.frame(LandUse=factor(
    c('Primary Vegetation','Secondary Vegetation',"Plantation forest", "Cropland", "Pasture"),
    levels=levels(data_samp_max$LandUse)))
  nd$RB_tnc <- factor(rb, levels = levels(data_samp_max$RB_tnc))
  nd$LogAbund = 0
  preds <- StatisticalModels::PredictGLMER(model = ma1$model, data = nd,
                                           se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  #75% CIs
  preds75 <- StatisticalModels::PredictGLMER(model = ma1$model, data = nd,
                                             se.fit = TRUE, seMultiplier = 1.15, randEffs = F)
  
  preds75$yplus75 <- ((exp(preds75$yplus)/exp(preds75$y[1]))*100)-100
  preds75$yminus75 <- ((exp(preds75$yminus)/exp(preds75$y[1]))*100)-100
  preds75$y75 <- ((exp(preds75$y)/exp(preds75$y[1]))*100)-100
  
  temp.df <- data.frame(RB = nd$RB_tnc, LU=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus, y75 = preds75$y75, lower75=preds75$yminus75, upper75=preds75$yplus75)
  
} )

preds2_mxa <- data.table::rbindlist(preds_mxa)  

#include sample size
for (i in 1:nrow(preds2_mxa)) {
  preds2_mxa$n[i] <- nrow(subset(data_samp_max[!is.na(data_samp_max$Total_abundance),], RB_tnc == preds2_mxa$RB[i] & LandUse == preds2_mxa$LU[i])) 
}

#add regional biome information for plotting
preds2_mxa = left_join(preds2_mxa, preds2_max[,c(1,2,10:13)], by = c('RB', 'LU'))
preds2_mxa$Biome_num <- as.numeric(preds2_mxa$Biome_num) 
preds2_mxa<-arrange(preds2_mxa, Biome_num, Realm)
#remove ranks with less than 25 detections
preds2_mxa= subset(preds2_mxa, n > 25)
preds2_mxa$Realm = factor(preds2_mxa$Realm, levels = c('Afrotropic', 'Australasia','Indo-Malay', 'Nearctic', 'Neotropic', 'Palearctic'))


figmxa <- ggplot(preds2_mxa, aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, group = Realm, colour = Realm)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.65), size = 2) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#FC8D62", "#A6D854", "#8DA0CB", "#FFD92F")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.65)) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.65), linewidth = 1) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100), hjust = 1.1, angle = 45,  data = preds2_mxa[preds2_mxa$upper > 100,], colour = 'black', position = position_dodge(width = 0.7), size = 1.5) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', "SV", "PF", "Pa", "Cr")) +
  scale_y_continuous(name = 'Change in Total abundance (%)') +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    text = element_text(size = 9, colour = 'black'),
    legend.position = 'top',
    legend.title = element_blank(),
    strip.text = element_text(hjust = 0)) +
  facet_wrap(~Biome, nrow = 3)

figmxa

ggsave('output/Eco_response/sampledmax_abundance_casestudy.pdf', width = 7.5, height = 14, units = 'cm')



# Land-Use intensity ------------------------------------------------------

#Comment from reviewer:
#I’m fairly concerned with grouping together plantations, pastures, and croplands and then dividing them solely based on their use intensity for a few reasons. First, I think there are inherently major differences in which species will be able to use each class. In some senses there is logic to grouping together the natural vegetation classes because the same species will likely prefer them both, but just some species will be more sensitive and need either older growth secondary or primary vegetation, but the directionality of each species’ response to natural cover is likely to be the same. I don’t think the same can be said of plantations, pasture, and cropland. 
# levels(data_samp_min$LU_UI)
# levels(data_samp_min$LU_UI_4) #pasture kept sep, PF and Cr together - this is used by Newbold
levels(data_samp_min$LU_UI_3) #this is the one that perhaps doesn't make sense
# levels(data_samp_min$LU_UI_5) #cropland and pasture together, forest sep. 
# levels(data_samp_min$Use_intensity)
#to try: 
#LU_UI
#LU_UI_5 - this makes sense as pasture/cropland are open habitats, forests closed.

data_samp_min <- data_samp_min %>%
  mutate(LU_UI_3 = recode_factor(LU_UI_3, 
                                 "Primary Vegetation_Minimal use" = "Primary Vegetation",
                                 "Primary Vegetation_Intense use" = "Primary Vegetation"),
         LU_UI = recode_factor(LU_UI,
                               "Primary Vegetation_Minimal use" = "Primary Vegetation",
                               "Primary Vegetation_Intense use" = "Primary Vegetation"),
         LU_UI_5 = recode_factor(LU_UI_5,
                                 "Primary Vegetation_Minimal use" = "Primary Vegetation",
                                 "Primary Vegetation_Intense use" = "Primary Vegetation"),
         LU_UI_4 = recode_factor(LU_UI_5,
                                 "Primary Vegetation_Minimal use" = "Primary Vegetation",
                                 "Primary Vegetation_Intense use" = "Primary Vegetation")
  )

#note - i tried with LU-UI. The results are not useful as the data is spread so thin, the impact of use intensity within regional biomes cannot be compared.

m4 = StatisticalModels::GLMER(modelData = data_samp_min, responseVar = 'LogRichness',
                              fitFamily = 'gaussian', fixedStruct = "LU_UI_3*RB_tnc", 
                              randomStruct = "(1|SS)+(1|SSB)", REML = F)
summary(m4$model)

preds_m4 <- apply(RBs, 1, FUN = function(rb) {
  nd <- data.frame(LU_UI_3=factor(
    c("Primary Vegetation", "Secondary Vegetation_Minimal use","Secondary Vegetation_Intense use",
      "Agriculture_Minimal use","Agriculture_Intense use"),
    levels=levels(data_samp_min$LU_UI_3)))
  nd$RB_tnc <- factor(rb, levels = levels(data_samp_min$RB_tnc))
  nd$LogRichness = 0
  preds <- StatisticalModels::PredictGLMER(model = m4$model, data = nd,
                                           se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(RB = nd$RB_tnc, LU_UI_3=nd$LU_UI_3, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
} )

preds2_m4 <- data.table::rbindlist(preds_m4)  

#include sample size
for (i in 1:nrow(preds2_m4)) {
  preds2_m4$n[i] <- nrow(subset(data_samp_min, RB_tnc == preds2_m4$RB[i] & LU_UI_3 == preds2_m4$LU_UI_3[i])) 
}

preds2_m4 <- preds2_m4 %>%
  mutate(Realm_code = substr(RB, 1,2),
         Biome_num = substr(RB, 3,3),
         Realm = recode_factor(Realm_code,
                               'IM' = 'Indo-Malay',
                               'PA' = 'Palearctic',
                               'AA' = 'Australasia',
                               'NA' = 'Nearctic',
                               'NT' = 'Neotropic',
                               'AT' = 'Afrotropic'),
         Biome = recode_factor(Biome_num,
                               "1" = 'Tropical & Subtropical Moist Broadleaf Forests',
                               "4" = 'Temperate Broadleaf & Mixed Forests',
                               "7" = 'Tropical & Subtropical Grasslands, Savannas & Shrublands'
         ))

preds2_m4$Biome_num <- as.numeric(preds2_m4$Biome_num) 
preds2_m4<-arrange(preds2_m4, Biome_num, Realm)
#remove ranks with less than 25 detections
preds2_m4= subset(preds2_m4, n > 25)
preds2_m4$Realm = factor(preds2_m4$Realm, levels = c('Afrotropic', 'Australasia','Indo-Malay', 'Nearctic', 'Neotropic', 'Palearctic'))


# figm4 <- ggplot(preds2_m4, aes(x = LU_UI_3, y = y, ymax = upper, ymin = lower, group = Realm, colour = Realm)) +
#   theme_bw() +
#   geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
#   geom_point(position = position_dodge(width = 0.6), size = 4) +
#   scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#FC8D62",, "#A6D854","#8DA0CB", "#FFD92F")) +
#   geom_errorbar(width = 0.2, position = position_dodge(width = 0.6)) +
#   coord_cartesian(ylim = c(-100,100)) +
#   scale_x_discrete(name = 'Land Use - Use Intensity', labels = c('PV', "SV-M", "SV-I", "Ag-M", "Ag-I")) +
#   scale_y_continuous(name = 'Change in Species Richness (%)') +
#   theme_bw() +
#   theme(#panel.grid.major = element_blank(),
#     #panel.grid.minor = element_blank(),
#     strip.background = element_blank()) +
#   facet_wrap(~Biome, nrow = 3)

#figm4



# Taxon model -------------------------------------------------------------

data_samp_max = droplevels(data_samp_max)

m1tx = StatisticalModels::GLMER(modelData = data_samp_max, responseVar = 'LogRichness',
                               fitFamily = 'gaussian', fixedStruct = "LandUse:RB_tnc:my_taxa", 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
summary(m1tx$model)

#vertebrate
preds_tx <- apply(RBs, 1, FUN = function(rb) {
  nd <- data.frame(LandUse=factor(rep(c('Primary Vegetation','Secondary Vegetation',"Plantation forest", "Cropland", "Pasture"), 3),
                                  levels=levels(data_samp_max$LandUse)))
  nd$RB_tnc <- factor(rb, levels = levels(data_samp_max$RB_tnc))
  nd$my_taxa = factor(c("Invertebrate", "Plant", "Vertebrate"), levels = levels(data_samp_max$my_taxa))
  nd$LogRichness = 0
  preds <- StatisticalModels::PredictGLMER(model = m1tx$model, data = nd,
                                           se.fit = TRUE, seMultiplier = 1.96)
  preds = cbind(nd, preds)
  preds = arrange(preds, LandUse, my_taxa)
  preds = preds %>%
    group_by(my_taxa) %>%
    mutate(yplus95 = ((exp(yplus)/exp(y[1]))*100)-100,
           yminus95 = ((exp(yminus)/exp(y[1]))*100)-100,
           y95 = ((exp(y)/exp(y[1]))*100)-100)
  #75% CIs
  preds75 <- StatisticalModels::PredictGLMER(model = m1t$model, data = nd,
                                             se.fit = TRUE, seMultiplier = 1.15, randEffs = F)
  preds75 = cbind(nd, preds75)
  preds75 = arrange(preds75, LandUse, my_taxa)
  preds75 = preds75 %>%
    group_by(my_taxa) %>%
    mutate(yplus75 = ((exp(yplus)/exp(y[1]))*100)-100,
           yminus75 = ((exp(yminus)/exp(y[1]))*100)-100,
           y75 = ((exp(y)/exp(y[1]))*100)-100)
  
  temp.df <- data.frame(my_taxa = preds$my_taxa, RB = preds$RB_tnc, LU=preds$LandUse,
                        y=preds$y95,lower=preds$yminus95,upper=preds$yplus95, 
                        y75 = preds75$y75, lower75=preds75$yminus75, upper75=preds75$yplus75)
  
} )

preds2_tx <- data.table::rbindlist(preds_tx)  

#include sample size
for (i in 1:nrow(preds2_tx)) {
  preds2_tx$n[i] <- nrow(subset(data_samp_max[!is.na(data_samp_max$Total_abundance),], RB_tnc == preds2_tx$RB[i] & LandUse == preds2_tx$LU[i])) 
}

#add regional biome information for plotting
preds2_tx <- preds2_tx %>%
  mutate(Realm_code = substr(RB, 1,2),
         Biome_num = substr(RB, 3,3),
         Realm = recode_factor(Realm_code,
                               'IM' = 'Indo-Malay',
                               'PA' = 'Palearctic',
                               'AA' = 'Australasia',
                               'NA' = 'Nearctic',
                               'NT' = 'Neotropic',
                               'AT' = 'Afrotropic'),
         Biome = recode_factor(Biome_num,
                               "1" = 'Tropical Forest',
                               "4" = 'Temperate Forest',
                               "7" = 'Tropical Grassland'
         ))


preds2_tx$Biome_num <- as.numeric(preds2_tx$Biome_num) 
preds2_tx<-arrange(preds2_tx, Biome_num, Realm)
#remove ranks with less than 25 detections
preds2_tx= subset(preds2_tx, n > 25)
preds2_tx$Realm = factor(preds2_tx$Realm, levels = c('Afrotropic', 'Australasia','Indo-Malay', 'Nearctic', 'Neotropic', 'Palearctic'))


figtx <- ggplot(preds2_tx, aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, group = Realm, colour = Realm)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.7), size = 1) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#FC8D62", "#A6D854", "#8DA0CB", "#FFD92F")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.7)) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), linewidth = 0.75) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100), hjust = 1.1, angle = 45,  data = preds2_tx[preds2_tx$upper > 100,], colour = 'black', position = position_dodge(width = 0.75), size = 1) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', "SV", "PF", "Pa", "Cr")) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    text = element_text(size = 9, colour = 'black'),
    legend.position = 'top', 
    legend.title = element_blank()) +
  facet_grid(Biome~my_taxa)


figtx

ggsave('output/Eco_response/downsampledmax_speciesrichness_taxon.pdf', width = 11.5, height = 11.5, units = 'cm')

# Supplementary tables & figures ------------------------------------------

#Table of down-sampled regional biomes

ecoregs_t = arrange(ecoregs_max, Biome, Realm)
TS4 = flextable(ecoregs_t[,c(5,4,6,8, 20,21,22, 9,13)])
TS4 = merge_v(TS4, j = ~ Biome)
TS4

typology <- data.frame(
  col_keys = TS4$col_keys,
  type = c("Regional Biome", "Regional Biome",
           'Sampling Effort', 'Sampling Effort', 
           'Taxonomic Sampling Effort', 'Taxonomic Sampling Effort', 'Taxonomic Sampling Effort',
           'Proportion of Ecoregions sampled', 'Proportion of Ecoregions sampled'),
  what = c("Biome", "Realm", "Studies", "Sites", "Invertebrates", "Plants", "Vertebrates", "Total ecoregions",
           "Proportion of Ecoregions sampled"),
  stringsAsFactors = FALSE )
TS4
formatflext <- function(ftable) {
  ftable <- set_header_df(ftable, mapping = typology, key = 'col_keys')
  ftable <- merge_h(ftable, part = "header")
  ftable <- merge_v(ftable, part = "header")
  ftable <- align(ftable, part = 'header', align = 'left')
  ftable <- theme_vanilla(ftable)
}

TS4 <- formatflext(TS4)
TS4 = align(TS4, part = 'header', align = 'left')
TS4 <- fix_border_issues(TS4)

#save_as_docx(f1, path = 'testtable.docx')
save_as_image(TS4, path = 'output/Eco_response/TableS5.png')