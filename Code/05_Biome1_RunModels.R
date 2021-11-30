## 03_Biome1RunModels

##Here I will subset the main dataset to biome 1 and run a model to look at variation in responses to land use change by realm. 


# 1. Packages
# 2. Data
# 3. Script
# 3.1 Prep and subset data
# 3.2 Species richness
#   Run models
#     Check AIC
#   Plot models
#      Overall biome LU and LU*UI
#      by realm LU and LU*UI
# 3.3 Abundance
#   Run Models
#     Check AIC
#   Plot models
#     Overall biome LU and LU*UI
#     by realm LU and LU*UI


# 1. Packages -------------------------------------------------------------
library(dplyr)
library(flextable)
library(ggplot2)
library(devtools)
library(StatisticalModels)
library(data.table) #rbindlist()
#function for getting predicted values from all models





# 2. DATA -----------------------------------------------------------------
data <- readRDS('Data/03_PREDICTSModelData.rds')
 dim(data)
# 3. SCRIPT ---------------------------------------------------------------
# 3.1 Prep and subset data ------------------------------------------------

#subset data to biome 1

Biome1 <- data[data$Biome=="Tropical & Subtropical Moist Broadleaf Forests",]
Biome1<- droplevels(Biome1)

#check distribution - make a table of realm/LU_UI

table(Biome1$Realm, Biome1$LandUse)
table(Biome1$Realm, Biome1$LandUse2)
table(Biome1$Realm, Biome1$LandUse3)
table(Biome1$Realm, Biome1$LandUse4)

table(Biome1$Realm, Biome1$LU_UI)
table(Biome1$Realm, Biome1$LU_UI_3, useNA = 'always')
table(Biome1$Realm, Biome1$LU_UI_4)
table(Biome1$LU_UI_3[Biome1$Realm=='Oceania'], useNA = 'always')
# I am choosing land use 3 because it gives 6 outputs and it allows easy comparison between different biome types 


#make PV useintensities the same
Biome1 <- Biome1 %>%
  mutate(LU_UI_3 = recode_factor(LU_UI_3, 
                                 "Primary Vegetation_Minimal use" = "Primary Vegetation",
                                 "Primary Vegetation_Intense use" = "Primary Vegetation"))

#remove oceania & australasia due to low sample size (going on 50 observation threshold)
Biome1 <- Biome1[Biome1$Realm!='Oceania',] %>%
  subset(Realm!='Australasia') %>%
  droplevels()

#subset rows with abundance values
Biome1_abund <- Biome1[!is.na(Biome1$Total_abundance),]
table(Biome1_abund$Realm, Biome1_abund$LU_UI_3, useNA = 'always')


#test random effect structure 
#check random effects structure
r0 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
r1 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)", REML = F)
r2 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SSB)", REML = F)
r3 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm'
                               , randomStruct = 0, REML = F)

AIC(r0$model, r1$model, r2$model)
#SS & SSB is best

#check best land use structure

L10 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l11 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse2*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l12 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse3*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l13 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse4*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
l14 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse5*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)



aic <- AIC(L10$model, l11$model, l12$model, l13$model, l14$model)
R2 <- NULL
R2[[1]] <- R2GLMER(L10$model)
R2[[2]] <- R2GLMER(l11$model)
R2[[3]] <- R2GLMER(l12$model)
R2[[4]] <- R2GLMER(l13$model)
R2[[5]] <- R2GLMER(l14$model)
R2 <- rbindlist(R2)

B1LUsel <- data.frame(Dataset = "Biome1", 
                      Response = 'SpeciesRichness', 
                      Fixef = 'LandUse*Realm', 
                      LandUseVar = c("LandUse", "LandUse2", "LandUse3", "LandUse4", "LandUse5"), 
                      AIC = aic$AIC, 
                      R2Marginal = R2$marginal, 
                      R2Conditional = R2$conditional, 
                      n = nrow(Biome1[!is.na(Biome1$Use_intensity),]),
                      n_RB = n_distinct(Biome1$RB_tnc)
                      )
B1LUsel$deltaAIC <-  B1LUsel$AIC - B1LUsel$AIC[1]

#Landuse 5 is marginally better than landuse 1, but lowest AIC is landuse 1.

# 3.2 Species Richness --------------------------------------------------
# Run Models: Species Richness ----------------------------------------------------------------
Biomemodels1 <- function(data, responseVar, LandUseVar, UseIntensityVar) {
  m <- NULL
  m[[1]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste(LandUseVar), 
                                     randomStruct = "(1|SS)+(1|SSB)", REML = F)
  
  m[[2]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste0(LandUseVar, "*Realm"), 
                                     randomStruct = "(1|SS)+(1|SSB)", REML = F)
  
  m[[3]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste0(UseIntensityVar), 
                                     randomStruct = "(1|SS)+(1|SSB)", REML = F)
  
  m[[4]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste0(UseIntensityVar, "*Realm"), 
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
  
  LandUseVardf <- c(rep(LandUseVar, 2), rep(UseIntensityVar, 2))
  
  modelresults <- data.frame('Dataset' = deparse(substitute(data)),
                             "Response" = responseVar,
                             "Fixef" = c("LandUse", "LandUse:Realm", "UseIntensity", "UseIntensity:Realm"),
                             "LandUseVar" = LandUseVardf,
                             "AIC" = aic$AIC,
                             "R2Marginal" = R2$marginal, 
                             "R2Conditional" = R2$conditional,
                             "n" = nrow,
                             "n_RB" = n_distinct(data$RB_tnc))
  modelresults$deltaAIC = modelresults$AIC - modelresults$AIC[1]
  
  return(modelresults)
}

B1sr<- Biomemodels1(data = Biome1[!is.na(Biome1$Use_intensity),], responseVar = 'LogRichness', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')
B1A <- Biomemodels1(data = Biome1[!is.na(Biome1$Use_intensity),], responseVar = 'LogAbund', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')

B1 <- rbind(B1sr, B1A)
B1 <- rbind(B1LUsel, B1)
write.csv(B1, "Figs/B1ModSel.csv", row.names = F)


# b0 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),], responseVar = "LogRichness",
#                                fitFamily = 'gaussian', fixedStruct = 'LandUse', 
#                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
# 
# 
# b1 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),], responseVar = "LogRichness",
#                                fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
#                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
# 
# b2 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "LogRichness",
#                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3', 
#                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
# 
# b3 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "LogRichness",
#                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3*Realm', 
#                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
# 
# 
# AIC(b0$model, b1$model, b2$model, b3$model)
# #adding realm improves model, land use/ use intensity is best combo
# 
# R2GLMER(b0$model)
# R2GLMER(b1$model)
# R2GLMER(b2$model)
# R2GLMER(b3$model)



# Plot models: Species Richness ----------------------------------------------------------



# Functions ---------------------------------------------------------------
#creating function to make predicted values from model - specific to land use variable

realmPredsLU_UI <- function(realm, model, data) {
  nd <- data.frame(LU_UI_3=factor(
    c("Primary Vegetation", 
      "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", 
      "Agriculture_Minimal use", "Agriculture_Intense use"),
    levels=levels(data$LU_UI_3)))
  nd$Realm <- factor(realm, levels = levels(data$Realm))
  nd$LogRichness = 0
  preds <- StatisticalModels::PredictGLMER(model = model, data = nd,
                                           se.fit = TRUE, seMultiplier = 1.96, randEffs = F)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(Realm = nd$Realm, LU=nd$LU_UI_3, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
}

realmPredsLU <- function(realm, model, data) {
  nd <- data.frame(LandUse=factor(
    c("Primary Vegetation", "Secondary Vegetation", "Cropland", "Pasture","Plantation forest"),
    levels=levels(data$LandUse)))
  nd$Realm <- factor(realm, levels = levels(data$Realm))
  nd$LogRichness = 0
  preds <- PredictGLMER(model = model, data = nd,
                        se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(Realm = nd$Realm, LU=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
}



#factors and order of factors in this model
LandUseFact <- c("Primary Vegetation", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use")

# overall biome - Land Use

#model:
b0 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),], 
                               responseVar = "LogRichness",
                               fitFamily = 'gaussian', 
                               fixedStruct = 'LandUse', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

nd <- data.frame(LandUse = factor(levels(b0$data$LandUse),
                                  levels = levels(b0$data$LandUse)),
                 LogRichness = 0)

preds <- StatisticalModels::PredictGLMER(model = b0$model, data = nd,
                                         se.fit = TRUE, seMultiplier = 1.96, randEffs = F)
#back transform
preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
temp.df <- data.frame(LandUse=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus)


figa <- ggplot(temp.df, aes(x = factor(LandUse), y = y, ymax = upper, ymin = lower)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.5), size = 3.5) + 
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.5)) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV', 'Cr', 'Pa' ,'PF')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  #geom_text(aes(x = LU, y = -100, label=paste('n = ', n)), 
  #          position=position_dodge(width=1.0), color = 'black'
  #) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  ggtitle("Biome 01: Tropical and Subtropical Moist Broadleaf Forest")

figa
ggsave(filename = 'Figs/attempt2/Biome1LU.png')


#Overall biome Land Use * Use intensity
#whole biome, land use * use intensity

#model
b2 <- StatisticalModels::GLMER(modelData = Biome1, 
                               responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LU_UI_3', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

nd <- data.frame(LU_UI_3 = factor(levels(b2$data$LU_UI_3),
                                  levels = levels(b2$data$LU_UI_3)),
                 LogRichness = 0)
preds <- StatisticalModels::PredictGLMER(model = b2$model, data = nd,
                                         se.fit = TRUE, seMultiplier = 1.96, randEffs = F)
preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
temp.df <- data.frame(LU_UI_3=nd$LU_UI_3, y=preds$y,lower=preds$yminus,upper=preds$yplus)
temp.df <- cbind(temp.df, LandUse, UseIntensity)

figb <- ggplot(temp.df, aes(x = factor(LU_UI_3), y = y, ymax = upper, ymin = lower)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(aes(shape = UseIntensity), position = position_dodge(width = 0.5), size = 3.5) + 
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.5)) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use/ Use Intensity', labels = c('PV_m', "PV_I", 'SV_M', "SV_I", 'Agr_M', "Agr_I")) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  #geom_text(aes(x = LU, y = -100, label=paste('n = ', n)), 
  #          position=position_dodge(width=1.0), color = 'black'
  #) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  ggtitle("Biome 01: Tropical and Subtropical Moist Broadleaf Forest")

figb
ggsave(filename = 'Figs/attempt2/Biome1LU_UI.png')


# by realm, Land Use

#model 

b1 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),], 
                               responseVar = "LogRichness",
                               fitFamily = 'gaussian', 
                               fixedStruct = 'LandUse*Realm',
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

#list realms in this biome
RealmB1 <- data.frame(Realm = unique(Biome1$Realm))

#test I've got the right realms
test <- sapply(RealmB1, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})


B1_LU<- apply(RealmB1, 1, FUN = realmPredsLU, model = b1$model, data = Biome1[!is.na(Biome1$Use_intensity),])
B1_LU <- data.table::rbindlist(B1_LU)
for ( i in 1:nrow(B1_LU)) {
  B1_LU$n[i] <- nrow(subset(Biome1[!is.na(Biome1$Use_intensity),], Realm == B1_LU$Realm[i] & LandUse == B1_LU$LU[i]))
}

level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Cropland", "Pasture","Plantation forest")

install.packages('wesanderson')
library(wesanderson)
library(RColorBrewer)
brewer.pal(n = 12, "Set2")


fig2 <- ggplot(B1_LU[B1_LU$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), size = 1) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV', 'Cr', 'Pa' ,'PF')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  #geom_text(aes(y = lower ), position="identity", color = 'black', size = 3, angle = 45) +        
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 25, colour = 'black'))
  
fig2
ggsave("Figs/attempt2/Biome1_RealmLU.png")

#colours = Afrotropic = '#DD8D29', Indo-Malay = '#E2D200', Neotropic = '#46ACC8'


  # + ggtitle("Biome 01: Tropical and Subtropical Moist Broadleaf Forest")



install.packages("ggrepel")
library(ggrepel)

# by realm, Land Use * Use intensity

#MODEL
b3 <- StatisticalModels::GLMER(modelData = Biome1,
                               responseVar = "LogRichness",
                               fitFamily = 'gaussian', 
                               fixedStruct = 'LU_UI_3*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

#find predicted values for b3 - use intensity/realm
B1_LUUI3 <- apply(RealmB1, 1, FUN = realmPredsLU_UI, model = b3$model, data = Biome1[!is.na(Biome1$Use_intensity),])
B1_LUUI3 <- data.table::rbindlist(B1_LUUI3)
#find n observations in each realm/land use combo
for ( i in 1:nrow(B1_LUUI3)) {
  B1_LUUI3$n[i] <- nrow(subset(Biome1[!is.na(Biome1$Use_intensity),], Realm == B1_LUUI3$Realm[i] & LU_UI_3 == B1_LUUI3$LU[i])) 
}

#list land use types & use intensity types, add to pred values df
LandUse = c("Primary Vegetation", 'Secondary Vegetation', 'Secondary Vegetation', 'Agriculture', 'Agriculture')
UseIntensity <- c("PV", 'Minimal Use', 'Instense Use', 'Minimal Use', 'Instense Use')
B1_LUUI3 <- cbind(B1_LUUI3, LandUse, UseIntensity)

level_order <- c("Primary Vegetation", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use")

fig1 <- ggplot(B1_LUUI3[B1_LUUI3$n > 25,], aes(x = factor(LU, levels = level_order), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(aes(shape = UseIntensity), position = position_dodge(width = 0.5), size = 5) + 
  geom_errorbar(width = 0, position = position_dodge(width = 0.5), size = 1) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 25, colour = 'black'))

fig1
ggsave("Figs/attempt2/Biome1_LU_UI3Realm.png")





realmPredsLU3 <- function(realm, model, data) {
  nd <- data.frame(LandUse3=factor(
    c("Primary Vegetation", "Secondary Vegetation", "Agriculture"),
    levels=levels(data$LandUse3)))
  nd$Realm <- factor(realm, levels = levels(data$Realm))
  nd$LogRichness = 0
  preds <- PredictGLMER(model = model, data = nd,
                        se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(Realm = nd$Realm, LU=nd$LandUse3, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
}

B1_LU<- apply(RealmB1, 1, FUN = realmPredsLU3, model = b1$model, data = Biome1[!is.na(Biome1$Use_intensity),])
B1_LU <- data.table::rbindlist(B1_LU)
for ( i in 1:nrow(B1_LU)) {
  B1_LU$n[i] <- nrow(subset(Biome1[!is.na(Biome1$Use_intensity),], Realm == B1_LU$Realm[i] & LandUse3 == B1_LU$LU[i]))
}

level_order_LU3 <- c("Primary Vegetation", "Secondary Vegetation", "Agriculture")



# Abundance ---------------------------------------------------------------
# Run Models: Abundance ---------------------------------------------------------


b1a0 <- StatisticalModels::GLMER(modelData = Biome1_abund[!is.na(Biome1_abund$Use_intensity),], responseVar = "LogAbund",
                         fitFamily = 'gaussian', fixedStruct = 'LandUse', 
                         randomStruct = "(1|SS)+(1|SSB)", REML = F)

b1a1 <- StatisticalModels::GLMER(modelData = Biome1_abund[!is.na(Biome1_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

b1a2 <- StatisticalModels::GLMER(modelData = Biome1_abund[!is.na(Biome1_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

b1a3 <- StatisticalModels::GLMER(modelData = Biome1_abund[!is.na(Biome1_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

AIC(ba0$model, ba1$model, ba2$model, ba3$model)

R2GLMER(ba0$model)
R2GLMER(ba1$model)
R2GLMER(ba2$model)
R2GLMER(ba3$model)

# Plot models: Abundance --------------------------------------------

#total biome: Land Use

# total biome: LandUse * Use intensity




# LandUse*UI --------------------------------------------------------------

realmPredsLU_UI_a <- function(realm, model, data) {
  nd <- data.frame(LU_UI_3=factor(
    c("Primary Vegetation", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use"),
    levels=levels(data$LU_UI_3)))
  nd$Realm <- factor(realm, levels = levels(data$Realm))
  nd$LogAbund = 0
  preds <- PredictGLMER(model = model, data = nd,
                        se.fit = TRUE, seMultiplier = 1.96, randEffs = F)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(Realm = nd$Realm, LU=nd$LU_UI_3, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
}

B1_LUUI3_a <- apply(RealmB1, 1, FUN = realmPredsLU_UI_a, model = b1a3$model, data = Biome1_abund[!is.na(Biome1_abund$Use_intensity),])

B1_LUUI3_a <- data.table::rbindlist(B1_LUUI3_a)
for ( i in 1:nrow(B1_LUUI3_a)) {
  B1_LUUI3_a$n[i] <- nrow(subset(Biome1, Realm == B1_LUUI3_a$Realm[i] & LU_UI_3 == B1_LUUI3_a$LU[i])) 
}


B1_LUUI3_a <- cbind(B1_LUUI3_a, LandUse, UseIntensity)

fig7 <- ggplot(B1_LUUI3_a[B1_LUUI3_a$n > 25,], aes(x = factor(LU, levels = level_order), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(aes(shape = UseIntensity), position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), size = 1) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        text = element_text(size = 25, colour = 'black'))

fig7
ggsave("Figs/attempt2/Biome1_LUUI3_abRealm.png")

#Biome 1 LUUI

nd <- data.frame(LU_UI_3 = factor(levels(ba2$data$LU_UI_3),
                                  levels = levels(ba2$data$LU_UI_3)),
                 LogAbund = 0)

preds <- StatisticalModels::PredictGLMER(model = ba2$model, data = nd,
                                         se.fit = TRUE, seMultiplier = 1.96, randEffs = F)


preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100

temp.df <- data.frame(LU_UI_3=nd$LU_UI_3, y=preds$y,lower=preds$yminus,upper=preds$yplus)

temp.df <- cbind(temp.df, LandUse, UseIntensity)



figc <- ggplot(temp.df, aes(x = factor(LU_UI_3), y = y, ymax = upper, ymin = lower)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(aes(shape = UseIntensity), position = position_dodge(width = 0.5), size = 3.5) + 
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.5)) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use/ Use Intensity', labels = c('PV_m', "PV_I", 'SV_M', "SV_I", 'Agr_M', "Agr_I")) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  #geom_text(aes(x = LU, y = -100, label=paste('n = ', n)), 
  #          position=position_dodge(width=1.0), color = 'black'
  #) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  ggtitle("Biome 01: Tropical and Subtropical Moist Broadleaf Forest")

figc
ggsave(filename = 'Figs/attempt2/Biome1LU_UI_abund.png')


##Land Use - Abundance 


# LandUse*Realm -----------------------------------------------------------

realmPredsLU_a <- function(realm, model, data) {
  nd <- data.frame(LandUse=factor(
    c("Primary Vegetation", "Secondary Vegetation", "Cropland", "Pasture","Plantation forest"),
    levels=levels(data$LandUse)))
  nd$Realm <- factor(realm, levels = levels(data$Realm))
  nd$LogAbund = 0
  preds <- PredictGLMER(model = model, data = nd,
                        se.fit = TRUE, seMultiplier = 1.96)  

  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(Realm = nd$Realm, LU=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
}

# LandUse*Realm

B1_LU_a <- apply(RealmB1, 1, FUN = realmPredsLU_a, model = ba1$model, data = Biome1)
B1_LU_a <- data.table::rbindlist(B1_LU_a)
for ( i in 1:nrow(B1_LU_a)) {
  B1_LU_a$n[i] <- nrow(subset(Biome1, Realm == B1_LU_a$Realm[i] & LandUse == B1_LU_a$LU[i])) 
}


fig8 <- ggplot(B1_LU_a[B1_LU_a$n > 20,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), size = 1) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV', 'Cr', 'Pa' ,'PF')) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        text = element_text(size = 25, colour = 'black')
        )

fig8
ggsave("Figs/attempt2/Biome1_LU_aRealm.png")





# Biome1:richness//realm//taxa -------------------------------------------

#remove fungi and change land use in to PV, SV, AGR_M, AGR_I

Biome1_a <- subset(Biome1, CommonTaxon_Phylum != 'Fungi')

Biome1_a <- droplevels(Biome1_a)

Biome1_a <- mutate(Biome1_a, LU_UI_6 = recode_factor(LU_UI_3,
           "Primary Vegetation_Minimal use" = "Primary Vegetation",
           "Primary Vegetation_Intense use" = "Primary Vegetation",
           "Secondary Vegetation_Minimal use" = "Secondary Vegetation",
           "Secondary Vegetation_Intense use" = 'Secondary Vegetation',
           "Agriculture_Minimal use" = "Agriculture_Minimal use", 
           "Agriculture_Intense use" = "Agriculture_Intense use"))
  


b5 <- GLMER(modelData = Biome1_a, responseVar = 'LogRichness', fitFamily = 'gaussian', fixedStruct = "LandUse3*CommonTaxon_Phylum*Realm", randomStruct = "(1|SS) + (1|SSB)", REML = F)
  

with(Biome1_a, table(LU_UI_6, Realm, CommonTaxon_Phylum))
library(sjPlot)

plot_model(b5$model, type = 'pred', terms = c('LandUse3',"CommonTaxon_Phylum", "Realm"))


realmPredsLU_UI_a <- function(realm, model, data) {
  nd <- data.frame(LU_UI_6=factor(
    c("Primary Vegetation", "Secondary Vegetation", "Agriculture_Minimal use", "Agriculture_Intense use"),
    levels=levels(data$LU_UI_6)))
  nd$Realm <- factor(realm, levels = levels(data$Realm))
  nd$CommonTaxon_Phylum <- factor(c("Invertebrate", "Plant", "Vertebrate"), levels = levels(data$CommonTaxon_Phylum))
  nd$LogAbund = 0
  preds <- PredictGLMER(model = model, data = nd,
                        se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(Realm = nd$Realm, Phylum = nd$CommonTaxon_Phylum, LU=nd$LU_UI_3, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
}

library(tidyr)
realmtaxa <- crossing(Realm = Biome1_a$Realm, CommonTaxon_Phylum = Biome1_a$CommonTaxon_Phylum) 


B1_LU3_tr <- apply(realmtaxa, 1, FUN = function(cross) {
  
  realm = cross[1]
  taxa = cross[2]
  
  nd <- data.frame(LandUse3=factor(
    c("Primary Vegetation", "Secondary Vegetation", "Agriculture"),
    levels=levels(b5$data$LandUse3)))
  nd$Realm <- factor(realm, levels = levels(b5$data$Realm))
  nd$CommonTaxon_Phylum <- factor(taxa, levels = levels(b5$data$CommonTaxon_Phylum))
  nd$LogRichness = 0
  
  preds <- PredictGLMER(model = b5$model, data = nd,
                        se.fit = TRUE,seMultiplier = 1.96)
  
  #create confidence intervals, with back-transformation
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(Realm = realm, Taxa = taxa, LU=nd$LandUse3,y=preds$y,lower=preds$yminus,upper=preds$yplus)
})


B1_LU3_tr <- data.table::rbindlist(B1_LU3_tr)

for ( i in 1:nrow(B1_LU3_tr)) {
  B1_LU3_tr$n[i] <- nrow(subset(Biome1_a, CommonTaxon_Phylum == B1_LU3_tr$Taxa[i] & Realm == B1_LU3_tr$Realm[i] & LandUse3 == B1_LU3_tr$LU[i])) 
}

fig8 <- ggplot(B1_LU3_tr[B1_LU3_tr$n > 5,], aes(x = factor(LU, levels = c("Primary Vegetation", "Secondary Vegetation", "Agriculture")), y = y, ymax = upper, ymin = lower, colour = Taxa)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.5), size = 3.5) + 
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.5)) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV' ,'Agr')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  geom_text(aes(x = LU, y = -100, label=paste('n = ', n)), 
            position=position_dodge(width=1.0), color = 'black'
  ) +
  facet_wrap(~Realm) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  ggtitle("Biome 01: Tropical and Subtropical Moist Broadleaf Forest")

fig8
ggsave("Figs/attempt2/Biome1_LU3_RealmTaxa.png")

write.csv(Biome1, "Figs/Maps/Biome1.csv", row.names = F)
