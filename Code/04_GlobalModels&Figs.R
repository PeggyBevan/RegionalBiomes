# Author: Peggy Bevan
# Date: 26/01/2021

#in this script, I will run a model (selected in 03_ModelSelection.r)
# i create visualisations that show how species richness varies with 
# regional biome & taxa.

# create figures for: biodiversity change with land use
# biodiversity change with land use and regional biome
# biodiversity change with land use and regional biome and taxa

# Packages ----------------------------------------------------------------

library(lme4) #glmer()
library(sjPlot)
library(data.table) #rbindlist()
library(ggplot2)
library(performance)
library(devtools)
library(dplyr)
# load PREDICTGLMER function 
source('FinalScriptsAndData/Code/05_PredictGLMERfunction.R')


install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions")

library(predictsFunctions)

install_github("timnewbold/StatisticalModels")

library(StatisticalModels)


# Data --------------------------------------------------------------------

#change this to RDS so it keeps stuff

data <- readRDS('FinalScriptsAndData/Data/03_PREDICTSModelData.rds')
dim(data)

# Code --------------------------------------------------------------------
#prep data
#remove data deficient biomes
data <- subset(data, Biome != "Tundra")
data <- subset(data, Biome != "Flooded Grasslands & Savannas")
data <- subset(data, Biome != "Mangroves")
data <- droplevels(data)

# Summary stats -----------------------------------------------------------

#number of studies
n_studies <- n_distinct(data$Study_name)
n_sources <- n_distinct(data$Source_ID)
n_sites <- n_distinct(data$SSBS)
#number of sources
#number of sites


# Selecting ecoregions ----------------------------------------------------

ecoregs <- read.csv("FinalScriptsAndData/Data/04_RBsummary.csv", stringsAsFactors = T)

ecoregs <- subset(ecoregs, Biome != "Tundra")
ecoregs <- subset(ecoregs, Biome != "Flooded Grasslands & Savannas")
ecoregs <- subset(ecoregs, Biome != "Mangroves")
ecoregs <- subset(ecoregs, n_sites > 1)

#thresholds 

#entire dataset, with mangroves, tundra, flooded grasslands removed. = data

#regional biomes that have <1 observations in PV and Ag
RB_LUth.1 <- ecoregs[ecoregs$PV.th1 == 1 & ecoregs$Ag.th1 == 1,]
RB_LUth.5 <- ecoregs[ecoregs$PV.th5 == 1 & ecoregs$Ag.th5 == 1,]
RB_LUth.25 <- ecoregs[ecoregs$PV.th25 == 1 & ecoregs$Ag.th25 == 1,]


#to create a subset in data:

data_LUth.1 <- filter(data, RB_tnc %in% RB_LUth.1$RB_tnc)
data_LUth.5 <- filter(data, RB_tnc %in% RB_LUth.5$RB_tnc)
data_LUth.25 <- filter(data, RB_tnc %in% RB_LUth.25$RB_tnc)


names(ecoregs)
# Run Models --------------------------------------------------------------

#check random effects structure
r0 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
r1 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)", REML = F)
r2 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SSB)", REML = F)
r3 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc'
                             , randomStruct = 0, REML = F)

AIC(r0$model, r1$model, r2$model)
#SS & SSB is best

#check best land use structure

L0 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l1 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse2*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l2 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse3*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l3 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse4*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
l4 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse5*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
  

aic <- AIC(L0$model, l1$model, l2$model, l3$model, l4$model)
R2 <- NULL
R2[[1]] <- R2GLMER(L0$model)
R2[[2]] <- R2GLMER(l1$model)
R2[[3]] <- R2GLMER(l2$model)
R2[[4]] <- R2GLMER(l3$model)
R2[[5]] <- R2GLMER(l4$model)
R2 <- rbindlist(R2)

GlobalLUsel <- data.frame(Dataset = "Global_LUth1", 
                      Response = 'SpeciesRichness', 
                      Fixef = 'LandUse*RB', 
                      LandUseVar = c("LandUse", "LandUse2", "LandUse3", "LandUse4", "LandUse5"), 
                      AIC = aic$AIC, 
                      R2Marginal = R2$marginal, 
                      R2Conditional = R2$conditional, 
                      n = nrow(data_LUth.1),
                      n_RB = n_distinct(data_LUth.1$RB_tnc)
)
GlobalLUsel$deltaAIC <- GlobalLUsel$AIC - GlobalLUsel$AIC[1]



m0 <- StatisticalModels::GLMER(modelData = data, responseVar = "LogRichness",
                                   fitFamily = 'gaussian', fixedStruct = 'LandUse', 
                                   randomStruct = "(1|SS)+(1|SSB)", REML = F)

m1 <- StatisticalModels::GLMER(modelData = data, responseVar = "LogRichness",
                                       fitFamily = 'gaussian', fixedStruct = 'LandUse*Biome', 
                                       randomStruct = "(1|SS)+(1|SSB)", REML = F)

m2 <- StatisticalModels::GLMER(modelData = data, responseVar = "LogRichness",
                                   fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                                   randomStruct = "(1|SS)+(1|SSB)", REML = F)

m3 <- StatisticalModels::GLMER(modelData = data, responseVar = "LogRichness",
                                   fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                                   randomStruct = "(1|SS)+(1|SSB)", REML = F)


#test if the final two models are overdispersed
GLMEROverdispersion(m1$model)
GLMEROverdispersion(m3$model)
#looks good

aic <- AIC(m0$model, m1$model, m2$model, m3$model)

R2 <- NULL
R2[[1]] <- R2GLMER(m0$model)
R2[[2]] <- R2GLMER(m1$model)
R2[[3]] <- R2GLMER(m2$model)
R2[[4]] <- R2GLMER(m3$model)
R2 <- rbindlist(R2)

names(m0$data[1])

modelresuts <- data.frame("Response" = c(names(m0$data[1]), names(m1$data[1]), names(m2$data[1]), names(m3$data[1])),
                          'Dataset' = 'data',
                          "AIC" = aic$AIC,
                          "R2Marginal" = R2$marginal, 
                          "R2Conditional" = R2$conditional)
responseVar = 'LogRichness'

runmodels1 <- function(data, responseVar, LandUseVar) {
  m <- NULL
  m[[1]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste(LandUseVar), 
                                     randomStruct = "(1|SS)+(1|SSB)", REML = F)
  
  m[[2]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste0(LandUseVar, "*Biome"), 
                                     randomStruct = "(1|SS)+(1|SSB)", REML = F)
  
  m[[3]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste0(LandUseVar, "*Realm"), 
                                     randomStruct = "(1|SS)+(1|SSB)", REML = F)
  
  m[[4]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste0(LandUseVar, "*RB_tnc"), 
                                     randomStruct = "(1|SS)+(1|SSB)", REML = F)
  
  aic <- AIC(m[[1]]$model, m[[2]]$model, m[[3]]$model, m[[4]]$model)
  nrow <- NULL
  nrow[1] <- nrow(m[[1]]$data)
  nrow[2] <- nrow(m[[2]]$data)
  nrow[3] <- nrow(m[[3]]$data)
  nrow[4] <- nrow(m[[4]]$data)
  
  R2 <- NULL
  R2[[1]] <- R2GLMER(m[[1]]$model)
  R2[[2]] <- R2GLMER(m[[2]]$model)
  R2[[3]] <- R2GLMER(m[[3]]$model)
  R2[[4]] <- R2GLMER(m[[4]]$model)
  R2 <- rbindlist(R2)
  
  modelresults <- data.frame('Dataset' = deparse(substitute(data)),
                             "Response" = responseVar,
                             "Fixef" = c("LandUse", "LandUse:Biome", "LandUse:Realm", "LandUse:RegionalBiome"),
                             "LandUseVar" = LandUseVar,
                             "AIC" = aic$AIC,
                             "R2Marginal" = R2$marginal, 
                             "R2Conditional" = R2$conditional,
                             "n" = nrow,
                             "n_RB" = n_distinct(data$RB_tnc))
  modelresults$deltaAIC = modelresults$AIC - modelresults$AIC[1]
  
  return(modelresults)
}

alldataM<- runmodels1(data = data, responseVar = 'LogRichness', LandUseVar = 'LandUse')
alldataMA <- runmodels1(data = data, responseVar = 'LogAbund', LandUseVar = 'LandUse')

#this is the dataset that should be used for everything
RB_LUth.1_Sresults <- runmodels1(data = data_LUth.1, responseVar = 'LogRichness', LandUseVar = 'LandUse') 
RB_LUth.1_Aresults <- runmodels1(data = data_LUth.1, responseVar = 'LogAbund', LandUseVar = 'LandUse')

#for upgrade appendix:
Globalmodsel <- rbind(GlobalLUsel, RB_LUth.1_Sresults, RB_LUth.1_Aresults)
write.csv(Globalmodsel, "FinalScriptsAndData/Figs/GlobalModSel.csv", row.names = F)

RB_LUth.5_Sresults <- runmodels1(data = data_LUth.5, responseVar = 'LogRichness') 
RB_LUth.5_Aresults <- runmodels1(data = data_LUth.5, responseVar = 'LogAbund')

RB_LUth.25_Sresults <- runmodels1(data = data_LUth.25, responseVar = 'LogRichness') 
RB_LUth.25_Aresults <- runmodels1(data = data_LUth.25, responseVar = 'LogAbund')


globalmods <- rbind(alldataM, alldataMA, RB_LUth.1_Sresults, RB_LUth.1_Aresults, RB_LUth.5_Sresults, RB_LUth.25_Aresults, RB_LUth.25_Sresults, RB_LUth.25_Aresults)

#are biome and regional biome a significant variables?

m1s <- GLMERSelect(modelData = data, responseVar = "LogRichness",
                            fitFamily = 'gaussian', fixedFactors = c('LandUse', 'Biome'), fixedInteractions = ("LandUse:Biome"), 
                            randomStruct = "(1|SS)+(1|SSB)")

m1s$stats

#all fixed effects are significant

m3s <- GLMERSelect(modelData = data, responseVar = "LogRichness",
                   fitFamily = 'gaussian', fixedFactors = c('LandUse', 'RB_tnc'), fixedInteractions = ("LandUse:RB_tnc"), 
                   randomStruct = "(1|SS)+(1|SSB)")

m3s$stats

#model selection suggests that Regional Biome is a stronger predictor than Biome, as it causes a stronger change in AIC.


# how i did the model in lme4
# model3<- lmer(LogRichness ~ LandUse*RB_tnc
#               + (1|SS) + (1|SSB), 
#               data = data, REML = F, 
#               control = lmerControl(optimizer = 'bobyqa', 
#                                     optCtrl = list(maxfun = 2e5)))



#are biome and regional biome a significant variables?

m1as <- GLMERSelect(modelData = data, responseVar = "LogAbund",
                   fitFamily = 'gaussian', fixedFactors = c('LandUse', 'Biome'), fixedInteractions = ("LandUse:Biome"), 
                   randomStruct = "(1|SS)+(1|SSB)")

m1as$stats

#all fixed effects are significant

m3as <- GLMERSelect(modelData = data, responseVar = "LogAbund",
                   fitFamily = 'gaussian', fixedFactors = c('LandUse', 'RB_tnc'), fixedInteractions = ("LandUse:RB_tnc"), 
                   randomStruct = "(1|SS)+(1|SSB)")

m3as$stats

#all fixed effects significant, the interaction between LU and RB is the strongest



# model3<- lmer(LogRichness ~ LandUse*CommonTaxon
#      + (1|SS) + (1|SSB), 
#      data = data, REML = F, 
#      control = lmerControl(optimizer = 'bobyqa', 
#                            optCtrl = list(maxfun = 2e5)))

# model4 <- lmer(LogRichness ~ LandUse*RB_tnc*CommonTaxon
#                      + (1|SS) + (1|SSB), 
#                      data = data, REML = F, 
#                      control = lmerControl(optimizer = 'bobyqa', 
#                                            optCtrl = list(maxfun = 2e5)))

#AIC selection


# Plot AIC ----------------------------------------------------------------
set.seed(20210430); # Set the randon number generator seed to ensure replicatable results

N = 100
sample_size = floor(nrow(data)*0.90) # Make our sample size 90% of the original data
studies = unique(data$SSB)
SSB_samplesize = floor(length(studies)*0.90)

sample_results = list() # Somewhere to put results for each sample
for (i in 1:N) {
  # For each attempt, randomly sample the data
  data_sample = sample_n(data, sample_size)
  #alternate way of sampling where you remove 90% of studies: 
  #study_sample = sample(studies, SSB_samplesize)
  #data_sample = subset(data, SSB %in% study_sample)
  #might have to run droplevels()
  # Run models using this random sample
  sample_dataM <- runmodels1(data = data_sample, responseVar = 'LogRichness', LandUseVar = "LandUse")
  sample_dataM$sample = i # Store sample number with these results
  
  sample_dataMA <- runmodels1(data = data_sample, responseVar = 'LogAbund',LandUseVar = "LandUse")
  sample_dataMA$sample = i # Store sample number with these results
  
  sample_results[[i]] = rbind(sample_dataM, sample_dataMA) # Store in the list
}

# Convert results list to data frame
sample_results_df = rbindlist(sample_results)

# Save results
write.csv(sample_results_df, file = "sample_results_df.csv")


# Fixing order of factors
sample_results_df$Fixef = factor(sample_results_df$Fixef, levels = c("LandUse", "LandUse:Realm", "LandUse:Biome", "LandUse:RegionalBiome"))

# Boxplot of deltaAIC across random samples...
facet_labels = c("A", "B")
ggplot(sample_results_df, aes(x = Fixef, y = deltaAIC, group = Fixef)) +
  geom_point(aes(shape = Fixef), position = "jitter") + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~Response, scales = "free_y") +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "none",
        text = element_text(size = 20, colour = 'black')) +
  scale_x_discrete(name = 'Fixed Effects', labels = c('LU', 'LU:Realm', 'LU:Biome', 'LU:RB')) +
  ylab(expression(paste(Delta, "AIC"))) +
  scale_y_continuous(limits = c(-1200,0))

ggsave("FinalScriptsAndData/Figs/dAIC2_sites.png")



# Plotting models ---------------------------------------------------------

# Species Richness ~ LandUse * Biome --------------------------------------------------------

# nd should be LandUse X Biome
Biome <- data.frame(Biome = unique(data$Biome))

test <- sapply(Biome, FUN = function(biome) {
  cat(paste0(biome,'\n'))
})

preds1 <- apply(Biome, 1, FUN = function(biome) {
  nd <- data.frame(LandUse=factor(
    c('Primary Vegetation','Secondary Vegetation',"Plantation forest", "Cropland", "Pasture"),
    levels=levels(data$LandUse)))
  nd$Biome <- factor(biome, levels = levels(data$Biome))
  nd$LogRichness = 0
  preds <- StatisticalModels::PredictGLMER(model = m1$model, data = nd,
                        se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(Biome = nd$Biome, LU=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
  } )

preds1 <- data.table::rbindlist(preds1)

for ( i in 1:nrow(preds1)) {
  preds1$n[i] <- nrow(subset(data, Biome == preds1$Biome[i] & LandUse == preds1$LU[i])) 
}


level_order <- c("Primary Vegetation", "Secondary Vegetation", "Cropland", "Pasture", "Plantation forest")

swr = function(string, nwrap=35) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

# Create line breaks in Year
preds1$Biome = swr(preds1$Biome)
preds1$gt10 = ifelse(test = preds1$n > 10, yes = 1, no = 0)
preds1$gt10 <- as.factor(preds1$gt10)


fig1 <- ggplot(preds1, aes(x = factor(LU, levels = level_order), y = y, ymax = upper, ymin = lower, group = Biome)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.5)) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use type', labels = c('PV', 'SV', 'Cr', 'Pa', 'PF')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  facet_wrap(~Biome, scales = 'free') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
fig1
ggsave("FinalScriptsAndData/Figs/attempt2/Fig07_LandUseBiome.png")

# Figure 2//Model 2 -------------------------------------------------------
##figure 2 / model 2
# Land Use * Regional Biome

data$RB_tnc <- as.factor(data$RB_tnc)
RBs <- data.frame(RB_tnc = unique(data$RB_tnc))

preds2 <- apply(RBs, 1, FUN = function(rb) {
  nd <- data.frame(LandUse=factor(
      c('Primary Vegetation','Secondary Vegetation',"Plantation forest", "Cropland", "Pasture"),
      levels=levels(data$LandUse)))
    nd$RB_tnc <- factor(rb, levels = levels(data$RB_tnc))
    nd$LogRichness = 0
    preds <- StatisticalModels::PredictGLMER(model = m3$model, data = nd,
                          se.fit = TRUE, seMultiplier = 1.96)
    
    preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
    preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
    preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
    
    temp.df <- data.frame(RB = nd$RB_tnc, LU=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus)
    
  } )

preds2 <- data.table::rbindlist(preds2)  

#add realm and biome to df  

for (i in 1:nrow(preds2)) {
  preds2$n[i] <- nrow(subset(data, RB_tnc == preds2$RB[i] & LandUse == preds2$LU[i])) 
}

preds2$gt10 = ifelse(test = preds2$n > 10, yes = 1, no = 0)
preds2$gt10 <- as.factor(preds2$gt10)

preds2 <- preds2 %>%
  mutate(Realm_code = substr(RB, 1,2),
                 Biome_num = ifelse(nchar(as.character(RB))<=3, 
                            substr(RB, 3,3),
                            substr(RB,3,4)
         ),
         Realm = recode_factor(Realm_code,
                               'IM' = 'Indo-Malay',
                               'PA' = 'Palearctic',
                               'AA' = 'Australasia',
                               'NA' = 'Nearctic',
                               'NT' = 'Neotropic',
                               'AT' = 'Afrotropic',
                               "OC" = 'Oceania'),
         Biome = recode_factor(Biome_num,
                               "1" = 'Tropical & Subtropical Moist Broadleaf Forests',
                               "2" = 'Tropical & Subtropical Dry Broadleaf Forests',
                               "3" = 'Tropical & Subtropical Coniferous Forests',
                               "4" = 'Temperate Broadleaf & Mixed Forests',
                               "5" = 'Temperate Conifer Forests',
                               "6" = "Boreal Forests/Taiga",
                               "7" = 'Tropical & Subtropical Grasslands, Savannas & Shrublands',
                               "8" = 'Temperate Grasslands, Savannas & Shrublands',
                               "10" = 'Montane Grasslands & Shrublands',
                               "11" = 'Tundra',
                               "12" = 'Mediterranean Forests, Woodlands & Scrub',
                               "13" = 'Deserts & Xeric Shrublands',
                               "14" = 'Mangroves'
                               ))

# preds2$Biome_num <- as.numeric(preds2$Biome_num) 
preds2<-arrange(preds2, Biome_num, Realm)

preds2$Biome = swr(preds2$Biome)
  
fig2 <- ggplot(preds2, aes(x = factor(LU, levels = level_order), y = y, ymax = upper, ymin = lower, group = Realm, color = Realm)) +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  coord_cartesian(ylim = c(-100,100)) +
  theme_classic() +
  theme() +
  scale_x_discrete(name = 'Land Use type', labels = c('PV', 'SV', 'Cr', 'Pa', 'PF')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  facet_wrap(~Biome, scales = 'free') + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

fig2


#STOP HERE --------------------
###fig 3 // model 3
# LandUse ** CommonTaxon
Taxa <- data.frame(Taxa = unique(data$CommonTaxon))

preds3 <- apply(Taxa, 1, FUN = function(taxa) {
  nd <- data.frame(LandUse=factor(
    c('Primary Vegetation','Secondary Vegetation',"Plantation forest", "Cropland", "Pasture"),
    levels=levels(data$LandUse)))
  nd$CommonTaxon <- factor(taxa, levels = levels(data$CommonTaxon))
  nd$LogRichness = 0
  preds <- PredictGLMER(model = model3, data = nd,
                        se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(Taxa = nd$CommonTaxon, LU=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
} )

preds3 <- data.table::rbindlist(preds3)

for ( i in 1:nrow(preds3)) {
  preds3$n[i] <- nrow(subset(data, CommonTaxon == preds3$Taxa[i] & LandUse == preds3$LU[i])) 
}

# #remove fungi as its over-fitting
# preds3 <- subset(preds3, Taxa != 'Fungi')
preds3$gt10 = ifelse(test = preds3$n > 10, yes = 1, no = 0)
preds3$gt10 <- as.factor(preds3$gt10)

preds3$Biome = swr(preds3$Biome)


fig3 <- ggplot(preds3, aes(x = factor(LU, levels = level_order), y = y, ymax = upper, ymin = lower, group = Taxa, color = Taxa)) +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  coord_cartesian(ylim = c(-100,100)) +
  theme_classic() +
  scale_x_discrete(name = 'Land Use type', labels = c('PV', 'SV', 'Cr', 'Pa', 'PF')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  theme() +
  facet_wrap(~Taxa, scales = 'free') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())


fig3


### fig 4 // model4
## LandUse * RB * Taxa

RBT <- crossing(RB_tnc = data$RB_tnc, CommonTaxon = data$CommonTaxon) 

test <- apply(RBT, 1, FUN = function(RB) {
  cat(RB)
})

preds4 <- apply(RBT, 1, FUN = function(cross) {
  
  RB = cross[1]
  taxa = cross[2]
  #cat(paste0(realm,'\n'))
  
  nd <- data.frame(LandUse=factor(
    c('Primary Vegetation','Secondary Vegetation',"Plantation forest", "Cropland", "Pasture"),
    levels=levels(data$LandUse)))
  
  nd$RB_tnc <- factor(RB,levels=levels(data$RB_tnc))
  
  nd$CommonTaxon <- factor(taxa, levels= levels(data$CommonTaxon))
  
  nd$LogRichness <- 0
  
  preds <- PredictGLMER(model = model4, data = nd,
                        se.fit = TRUE,seMultiplier = 1.96)
  
  #create confidence intervals, with back-transformation
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(RB = RB, Taxa = taxa, LU=nd$LandUse,y=preds$y,lower=preds$yminus,upper=preds$yplus)
}
)

preds4 <- data.table::rbindlist(preds4)

for ( i in 1:nrow(preds4)) {
  preds4$n[i] <- nrow(subset(data, CommonTaxon == preds4$Taxa[i] & RB_tnc == preds4$RB[i] & LandUse == preds4$LU[i])) 
}

preds4 <- preds4 %>%
  mutate(Realm_code = substr(RB, 1,2),
         Biome_num = ifelse(nchar(as.character(RB))<=3, 
                          substr(RB, 3,3),
                          substr(RB,3,4)
       ),
       Realm = recode_factor(Realm_code,
                             'IM' = 'Indo-Malay',
                             'PA' = 'Palearctic',
                             'AA' = 'Australasia',
                             'NA' = 'Nearctic',
                             'NT' = 'Neotropic',
                             'AT' = 'Afrotropic',
                             'OC' = 'Oceania'),
       Biome = recode_factor(Biome_num,
                             "1" = 'Tropical & Subtropical Moist Broadleaf Forests',
                             "4" = 'Temperate Broadleaf & Mixed Forests',
                             "12" = 'Mediterranean Forests, Woodlands & Scrub',
                             "7" = 'Tropical & Subtropical Grasslands, Savannas & Shrublands',
                             "5" = 'Temperate Conifer Forests',
                             "8" = 'Temperate Grasslands, Savannas & Shrublands',
                             "13" = 'Deserts & Xeric Shrublands',
                             "14" = 'Mangroves',
                             "2" = 'Tropical & Subtropical Dry Broadleaf Forests',
                             "10" = 'Montane Grasslands & Shrublands',
                             "3" = 'Tropical & Subtropical Coniferous Forests',
                             "6" = 'Boreal Forests/Taiga',
                             "9" = 'Flooded Grasslands & Savannas',
                             "11" = 'Tundra',
                             "98" = "Inland Water", 
                             "99" = "Rock and Ice"
       ))

preds4$Biome = swr(preds4$Biome)

#just plot biomes that are having more sever reactions - 
# meditteranean (12), tropical forest(1), tropical grasslands (7),
#badbiomes$biomes <- c("Tropical & Subtropical Moist\nBroadleaf Forests", "Mediterranean Forests, Woodlands &\nScrub", "Tropical & Subtropical Grasslands,\nSavannas & Shrublands")

preds4plot <- filter(preds4, Biome == "Tropical & Subtropical Moist\nBroadleaf Forests" & n > 10)


fig4 <- ggplot(preds4plot[preds4plot$n > 10], aes(x = factor(LU, levels = level_order), y = y, ymax = upper, ymin = lower, group = Taxa, color = Taxa)) +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  facet_grid(~Realm) +
  coord_cartesian(ylim = c(-100,100)) +
  theme_classic() +
  scale_x_discrete(name = 'Land Use type', labels = c('PV', 'SV', 'Cr', 'Pa', 'PF')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  theme() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

fig4

ggsave('FinalScriptsAndData/Figs/06_GraphBiome.png', plot = fig1, height = 10, width = 12)
# Fig X. Predicted species richnes change compared to the reference level at primary vegetation (PV) under different land use scenarios for each biome. Each point is the predicted value with error bars. Red points are those where n < 1 so should be discounted. Land Use scenarios are Primary vegetation (PV), Secondary Vegetation (SV), Cropland (Cr), Pasture (Pa) and Plantation Forest (PF)

ggsave('FinalScriptsAndData/Figs/07_GraphRegBiome.png', plot = fig2, width = 12)
# Fig X. Predicted species richness change compared to the reference level at primary vegetation for each region in each biome. Land Use scenarios are Primary vegetation (PV), Secondary Vegetation (SV), Cropland (Cr), Pasture (Pa) and Plantation Forest (PF)

ggsave('FinalScriptsAndData/Figs/08_GraphTaxaBiome.png', plot = fig3, width = 10)
# Fig X. Predicted species richness change compared to the reference level at primary vegetation for each taxa. Land Use scenarios are Primary vegetation (PV), Secondary Vegetation (SV), Cropland (Cr), Pasture (Pa) and Plantation Forest (PF)

ggsave('FinalScriptsAndData/Figs/09_GraphBiome01Taxa.png', plot = fig4, width = 12, height = 6)
# Fig X. Predicted species richness change compared to the reference level at primary vegetation for each taxon group in every regional biome of Tropical Moist Forest (Biome 01. Land Use scenarios are Primary vegetation (PV), Secondary Vegetation (SV), Cropland (Cr), Pasture (Pa) and Plantation Forest (PF)
