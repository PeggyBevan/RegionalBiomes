##04_Biome4_RunModels

##Here I will subset the main dataset to biome 1 and run a model to look at variation in responses to land use change by realm. 


# 1. Packages
# 2. Data
# 3. Script
# 3.1 Prep and subset data
# Species richness
#   Run models
#     Check AIC
#   Plot models
#      Overall biome LU and LU*UI
#      by realm LU and LU*UI
# Abundance
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
library(wesanderson)
# 2. DATA -----------------------------------------------------------------
data <- readRDS('FinalScriptsAndData/Data/03_PREDICTSModelData.rds')
dim(data)
# 3. SCRIPT ---------------------------------------------------------------


# 3.1 Prep and subset data ------------------------------------------------

#subset data to biome 4

Biome4 <- data[data$Biome=="Temperate Broadleaf & Mixed Forests",]

#check distribution - make a table of realm/LU_UI



Biome4 <- droplevels(Biome4)


#subset data to biome 4

#check distribution - make a table of realm/LU_UI

table(Biome4$Realm, Biome4$LandUse)
table(Biome4$Realm, Biome4$LandUse2)
table(Biome4$Realm, Biome4$LandUse3)
table(Biome4$Realm, Biome4$LandUse4)

table(Biome4$Realm, Biome4$LU_UI)
table(Biome4$Realm, Biome4$LU_UI_3, useNA = 'always')
table(Biome4$Realm, Biome4$LU_UI_4)




table(Biome4_abund$Realm, Biome4_abund$LU_UI_3, useNA = 'always')

#make PV I & M the same
Biome4 <- Biome4 %>%
  mutate(LU_UI_3 = recode_factor(LU_UI_3, 
                                 "Primary Vegetation_Minimal use" = "Primary Vegetation",
                                 "Primary Vegetation_Intense use" = "Primary Vegetation"))

#remove Indo-Malay for species richness models as there are no obs in Agriculture
Biome4 <- Biome4[Biome4$Realm!='Indo-Malay',] %>%
  droplevels()

Biome4_abund <- Biome4[!is.na(Biome4$Total_abundance),]

#test random effects
r0 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
r1 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)", REML = F)
r2 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SSB)", REML = F)
r3 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm'
                               , randomStruct = 0, REML = F)

AIC(r0$model, r1$model, r2$model)
#SS & SSB is best

#check best land use structure

L40 <- StatisticalModels::GLMER(modelData = Biome4[!is.na(Biome4$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l41 <- StatisticalModels::GLMER(modelData = Biome4[!is.na(Biome4$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse2*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l42 <- StatisticalModels::GLMER(modelData = Biome4[!is.na(Biome4$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse3*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l43 <- StatisticalModels::GLMER(modelData = Biome4[!is.na(Biome4$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse4*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
l44 <- StatisticalModels::GLMER(modelData = Biome4[!is.na(Biome4$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse5*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)



aic <- AIC(L40$model, l41$model, l42$model, l43$model, l44$model)
R2 <- NULL
R2[[1]] <- R2GLMER(L40$model)
R2[[2]] <- R2GLMER(l41$model)
R2[[3]] <- R2GLMER(l42$model)
R2[[4]] <- R2GLMER(l43$model)
R2[[5]] <- R2GLMER(l44$model)
R2 <- rbindlist(R2)

B4LUsel <- data.frame(Dataset = "Biome4", 
                      Response = 'SpeciesRichness', 
                      Fixef = 'LandUse*Realm', 
                      LandUseVar = c("LandUse", "LandUse2", "LandUse3", "LandUse4", "LandUse5"), 
                      AIC = aic$AIC, 
                      R2Marginal = R2$marginal, 
                      R2Conditional = R2$conditional, 
                      n = nrow(Biome4[!is.na(Biome4$Use_intensity),]),
                      n_RB = n_distinct(Biome4$RB_tnc)
)
B4LUsel$deltaAIC <- B4LUsel$AIC - B4LUsel$AIC[1]



# 4.1 Species Richness: Run Models ----------------------------------------------------------

#model selection
#using biome models function from biome1 script
B4sr<- Biomemodels1(data = Biome4[!is.na(Biome4$Use_intensity),], responseVar = 'LogRichness', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')
B4A <- Biomemodels1(data = Biome4[!is.na(Biome4$Use_intensity),], responseVar = 'LogAbund', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')

B4 <- rbind(B4sr, B4A)
B4 <- rbind(B4LUsel, B4)
write.csv(B4, "FinalScriptsAndData/Figs/B4ModSel.csv", row.names = F)


b0 <- StatisticalModels::GLMER(modelData = Biome4[!is.na(Biome4$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)


b1 <- StatisticalModels::GLMER(modelData = Biome4[!is.na(Biome4$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

b2 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LU_UI_3', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

b3 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LU_UI_3*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)


AIC(b0$model, b1$model, b2$model, b3$model)
#adding realm improves model, land use/ use intensity is best combo

R2GLMER(b0$model)
R2GLMER(b1$model)
R2GLMER(b2$model)
R2GLMER(b3$model)



# 4.3 Plot model ----------------------------------------------------------

RealmB4 <- data.frame(Realm = unique(Biome4$Realm))

test <- sapply(RealmB4, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})

B4_LUUI3 <- apply(RealmB4, 1, FUN = realmPredsLU_UI, model = b3$model, data = Biome4)
B4_LUUI3 <- data.table::rbindlist(B4_LUUI3)
for ( i in 1:nrow(B4_LUUI3)) {
  B4_LUUI3$n[i] <- nrow(subset(Biome4, Realm == B4_LUUI3$Realm[i] & LU_UI_3 == B4_LUUI3$LU[i])) 
}
LandUse = c("Primary Vegetation", 'Secondary Vegetation', 'Secondary Vegetation', 'Agriculture', 'Agriculture')
UseIntensity <- c('PV', 'Minimal Use', 'Instense Use', 'Minimal Use', 'Instense Use')
B4_LUUI3 <- cbind(B4_LUUI3, LandUse, UseIntensity)

level_order <- c("Primary Vegetation", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use")

fig3 <- ggplot(B4_LUUI3[B4_LUUI3$n > 25,], aes(x = factor(LU, levels = level_order), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  scale_colour_manual(values = c("#E78AC3", "#A6D854", '#8DA0CB', "#FFD92F")) +
  geom_point(aes(shape = UseIntensity), position = position_dodge(width = 0.6), size = 5) + 
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), size = 1) +
  coord_cartesian(ylim = c(-100,300)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        text = element_text(size = 25, colour = 'black'))

fig3
ggsave("FinalScriptsAndData/Figs/attempt2/Biome4_LU_UI3Realm.png")


## Biome 4: Land Use * Realm


B4_LU<- apply(RealmB4, 1, FUN = realmPredsLU, model = b1$model, data = Biome4[!is.na(Biome4$Use_intensity),])
B4_LU <- data.table::rbindlist(B4_LU)
for ( i in 1:nrow(B4_LU)) {
  B4_LU$n[i] <- nrow(subset(Biome4[!is.na(Biome4$Use_intensity),], Realm == B4_LU$Realm[i] & LandUse == B4_LU$LU[i]))
}

level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Cropland", "Pasture","Plantation forest")

brewer.pal(n = 12, "Set2")

fig4 <- ggplot(B4_LU[B4_LU$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#E78AC3", "#A6D854", '#8DA0CB', "#FFD92F")) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), size = 1) +
  coord_cartesian(ylim = c(-100,200)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV', 'Cr', 'Pa' ,'PF')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  #geom_text(aes(x = LU, y = -100, label=paste('n = ', n)), 
  #          position=position_dodge(width=1.0), color = 'black'
  #) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 25, colour = 'black'))
  #+ggtitle("Biome 04: Temperate Broadleaf & Mixed Forests")

fig4
ggsave("FinalScriptsAndData/Figs/attempt2/Biome4_RealmLU.png")


#whole biome, just land use 

nd <- data.frame(LandUse = factor(levels(b0$data$LandUse),
                                  levels = levels(b0$data$LandUse)),
                 LogRichness = 0)

preds <- StatisticalModels::PredictGLMER(model = b0$model, data = nd,
                                         se.fit = TRUE, seMultiplier = 1.96, randEffs = F)


preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100

temp.df <- data.frame(LandUse=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus)


figa <- ggplot(temp.df, aes(x = factor(LandUse, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = LandUse)) +
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
  ggtitle("Biome 04: Temperate Broadleaf and Mixed Forests")

figa
ggsave(filename = 'FinalScriptsAndData/Figs/attempt2/Biome4LU.png')


#whole biome, land use * use intensity

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



figb <- ggplot(temp.df, aes(x = factor(LU_UI_3), y = y, ymax = upper, ymin = lower, colour = LandUse)) +
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
  ggtitle("Biome 04: Temperate Broadleaf and Mixed Forests")

figb
ggsave(filename = 'FinalScriptsAndData/Figs/attempt2/Biome4LU_UI.png')


# Abundance ---------------------------------------------------------------
# Run Models: Abundance ---------------------------------------------------------


ba0 <- StatisticalModels::GLMER(modelData = Biome4_abund[!is.na(Biome4_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LandUse', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

ba1 <- StatisticalModels::GLMER(modelData = Biome4_abund[!is.na(Biome4_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

ba2 <- StatisticalModels::GLMER(modelData = Biome4_abund[!is.na(Biome4_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

b4a3 <- StatisticalModels::GLMER(modelData = Biome4_abund[!is.na(Biome4_abund$Use_intensity),], responseVar = "LogAbund",
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

#by realm Land Use
#by realm Land Use * Use intensity


realmPredsLU_UI_a <- function(realm, model, data) {
  nd <- data.frame(LU_UI_3=factor(
    c("Primary Vegetation", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use"),
    levels=levels(data$LU_UI_3)))
  nd$Realm <- factor(realm, levels = levels(data$Realm))
  nd$LogAbund = 0
  preds <- PredictGLMER(model = model, data = nd,
                        se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  temp.df <- data.frame(Realm = nd$Realm, LU=nd$LU_UI_3, y=preds$y,lower=preds$yminus,upper=preds$yplus)
  
}

B4_LUUI3_a <- apply(RealmB4, 1, FUN = realmPredsLU_UI_a, model = b4a3$model, data = Biome4_abund[!is.na(Biome4_abund$Use_intensity),])
B4_LUUI3_a <- data.table::rbindlist(B4_LUUI3_a)
for ( i in 1:nrow(B4_LUUI3_a)) {
  B4_LUUI3_a$n[i] <- nrow(subset(Biome4, Realm == B4_LUUI3_a$Realm[i] & LU_UI_3 == B4_LUUI3_a$LU[i])) 
}


B4_LUUI3_a <- cbind(B4_LUUI3_a, LandUse, UseIntensity)



# Abundance: UseIntensity*Realm - NOT SIGNIFICANT -------------------------
fig7a <- ggplot(B4_LUUI3_a[B4_LUUI3_a$n > 25,], aes(x = factor(LU, levels = level_order), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(aes(shape = UseIntensity), position = position_dodge(width = 0.6), size = 5) +
  scale_colour_manual(values = c("#E78AC3", "#A6D854", '#8DA0CB', "#FFD92F")) +
  geom_point(aes(shape = UseIntensity), position = position_dodge(width = 0.6), size = 5) + 
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), size = 1) +
  coord_cartesian(ylim = c(-100,300)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        text = element_text(size = 25, colour = 'black'))

fig7a
ggsave("FinalScriptsAndData/Figs/attempt2/Fig07_Biome4_LU_UI_3_abundanceRealm.png")



# LandUse&Realm -----------------------------------------------------------

B4_LU_a <- apply(RealmB4, 1, FUN = realmPredsLU_a, model = ba1$model, data = Biome4_abund[!is.na(Biome4_abund$Use_intensity),])
B4_LU_a <- data.table::rbindlist(B4_LU_a)
for ( i in 1:nrow(B4_LU_a)) {
  B4_LU_a$n[i] <- nrow(subset(Biome4, Realm == B4_LU_a$Realm[i] & LandUse == B4_LU_a$LU[i])) 
}


figB4LU <- ggplot(B4_LU_a[B4_LU_a$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#E78AC3", "#A6D854", '#8DA0CB', "#FFD92F")) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), size = 1) +
  coord_cartesian(ylim = c(-100,650)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV', 'Cr', 'Pa' ,'PF')) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        text = element_text(size = 25, colour = 'black')
  )

figB4LU
ggsave("FinalScriptsAndData/Figs/attempt2/Biome4_LU_aRealm.png")



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



figc <- ggplot(temp.df, aes(x = factor(LU_UI_3), y = y, ymax = upper, ymin = lower, colour = LandUse)) +
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
ggsave(filename = 'FinalScriptsAndData/Figs/attempt2/Biome4LU_UI_abund.png')



# Biome4:richness//realm//taxa -------------------------------------------

#remove fungi and change land use in to PV, SV, AGR_M, AGR_I

Biome4_a <- subset(Biome4, CommonTaxon_Phylum != 'Fungi')

Biome4_a <- droplevels(Biome4_a)

Biome4_a <- mutate(Biome4_a, LU_UI_6 = recode_factor(LU_UI_3,
                                                     "Primary Vegetation_Minimal use" = "Primary Vegetation",
                                                     "Primary Vegetation_Intense use" = "Primary Vegetation",
                                                     "Secondary Vegetation_Minimal use" = "Secondary Vegetation",
                                                     "Secondary Vegetation_Intense use" = 'Secondary Vegetation',
                                                     "Agriculture_Minimal use" = "Agriculture_Minimal use", 
                                                     "Agriculture_Intense use" = "Agriculture_Intense use"))



b5 <- GLMER(modelData = Biome4_a, responseVar = 'LogRichness', fitFamily = 'gaussian', fixedStruct = "LandUse3*CommonTaxon_Phylum*Realm", randomStruct = "(1|SS) + (1|SSB)", REML = F)


with(Biome4_a, table(LU_UI_6, Realm, CommonTaxon_Phylum))
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
realmtaxa <- crossing(Realm = Biome4_a$Realm, CommonTaxon_Phylum = Biome4_a$CommonTaxon_Phylum) 


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
  B1_LU3_tr$n[i] <- nrow(subset(Biome4_a, CommonTaxon_Phylum == B1_LU3_tr$Taxa[i] & Realm == B1_LU3_tr$Realm[i] & LandUse3 == B1_LU3_tr$LU[i])) 
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
ggsave("FinalScriptsAndData/Figs/attempt2/Biome4_LU3_RealmTaxa.png")



write.csv(Biome4, "FinalScriptsAndData/Figs/Maps/Biome4.csv", row.names = F)



