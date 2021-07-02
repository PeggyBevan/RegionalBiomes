### 05_Biome7_RunModels


##Here I will subset the main dataset to biome 7 and run a model to look at variation in responses to land use change by realm. 


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
# 2. DATA -----------------------------------------------------------------
data <- readRDS('FinalScriptsAndData/Data/03_PREDICTSModelData.rds')
dim(data)
# 3. SCRIPT ---------------------------------------------------------------


# 3.1 Prep and subset data ------------------------------------------------

#subset data to biome 7

Biome7 <- data[data$Biome=="Tropical & Subtropical Grasslands, Savannas & Shrublands",] %>%
  droplevels()

#check distribution - make a table of realm/LU_UI

table(Biome7$Realm, Biome7$LandUse)
table(Biome7$Realm, Biome7$LandUse2)
table(Biome7$Realm, Biome7$LandUse3)
table(Biome7$Realm, Biome7$LandUse4)

table(Biome7$Realm, Biome7$LU_UI)
table(Biome7$Realm, Biome7$LU_UI_3, useNA = 'always')
table(Biome7$Realm, Biome7$LU_UI_4)

Biome7 <- Biome7 %>%
  mutate(LU_UI_3 = recode_factor(LU_UI_3, 
                                 "Primary Vegetation_Minimal use" = "Primary Vegetation",
                                 "Primary Vegetation_Intense use" = "Primary Vegetation"))

Biome7_abund <- Biome7[!is.na(Biome7$Total_abundance),]
table(Biome7_abund$Realm, Biome7_abund$LU_UI_3, useNA = 'always')

#check random effects & land use structure

#test random effects
r0 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
r1 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)", REML = F)
r2 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SSB)", REML = F)
r3 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm'
                               , randomStruct = 0, REML = F)

AIC(r0$model, r1$model, r2$model)
#SS & SSB is best

#check best land use structure

L70 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l71 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse2*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l72 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse3*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

l73 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse4*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
l74 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse5*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)



aic <- AIC(L70$model, l71$model, l72$model, l73$model, l74$model)
R2 <- NULL
R2[[1]] <- R2GLMER(L70$model)
R2[[2]] <- R2GLMER(l71$model)
R2[[3]] <- R2GLMER(l72$model)
R2[[4]] <- R2GLMER(l73$model)
R2[[5]] <- R2GLMER(l74$model)
R2 <- rbindlist(R2)

B7LUsel <- data.frame(Dataset = "Biome7", 
                      Response = 'SpeciesRichness', 
                      Fixef = 'LandUse*Realm', 
                      LandUseVar = c("LandUse", "LandUse2", "LandUse3", "LandUse4", "LandUse5"), 
                      AIC = aic$AIC, 
                      R2Marginal = R2$marginal, 
                      R2Conditional = R2$conditional, 
                      n = nrow(Biome7[!is.na(Biome7$Use_intensity),]),
                      n_RB = n_distinct(Biome7$RB_tnc)
)
B7LUsel$deltaAIC <- B7LUsel$AIC - B7LUsel$AIC[1]



# 4.1 Species Richness: Run Models ----------------------------------------------------------

#model selection
#using biome models function from biome1 script
B7sr<- Biomemodels1(data = Biome7[!is.na(Biome7$Use_intensity),], responseVar = 'LogRichness', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')
B7A <- Biomemodels1(data = Biome7[!is.na(Biome7$Use_intensity),], responseVar = 'LogAbund', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')

B7 <- rbind(B7sr, B7A)
B7 <- rbind(B7LUsel, B7)
write.csv(B7, "FinalScriptsAndData/Figs/B7ModSel.csv", row.names = F)





# Species richness --------------------------------------------------------
# 3.2 Run Models ----------------------------------------------------------





# 5.2 Run Models ----------------------------------------------------------

#model selection


b70 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)


b71 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

b72 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LU_UI_3', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

b73 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LU_UI_3*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)


AIC(b70$model, b71$model, b72$model, b73$model)
#adding realm improves model, land use/ use intensity is best combo

R2GLMER(b0$model)
R2GLMER(b1$model)
R2GLMER(b2$model)
R2GLMER(b3$model)

# 3.3 Plot model ----------------------------------------------------------

RealmB7 <- data.frame(Realm = unique(Biome7$Realm))

test <- sapply(RealmB7, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})


B7_LUUI3 <- apply(RealmB7, 1, FUN = realmPredsLU_UI, model = b73$model, data = Biome7[!is.na(Biome7$Use_intensity),])
B7_LUUI3 <- data.table::rbindlist(B7_LUUI3)
for ( i in 1:nrow(B7_LUUI3)) {
  B7_LUUI3$n[i] <- nrow(subset(Biome7, Realm == B7_LUUI3$Realm[i] & LU_UI_3 == B7_LUUI3$LU[i])) 
}
B7_LUUI3 <- cbind(B7_LUUI3, LandUse, UseIntensity)


fig5 <- ggplot(B7_LUUI3[B7_LUUI3$n > 25,], aes(x = factor(LU, levels = level_order), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(aes(shape = UseIntensity), position = position_dodge(width = 0.5), size = 5) + 
  geom_errorbar(width = 0, position = position_dodge(width = 0.5), size = 1) +
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 25, colour = 'black'))

fig5
ggsave("FinalScriptsAndData/Figs/attempt2/Biome7_LUUI3Realm.png")

##Biome7: Land Use * Realm (model 2)

B7_LU<- apply(RealmB7, 1, FUN = realmPredsLU, model = b71$model, data = Biome7[!is.na(Biome7$Use_intensity),])
B7_LU <- data.table::rbindlist(B7_LU)
for ( i in 1:nrow(B7_LU)) {
  B7_LU$n[i] <- nrow(subset(Biome7[!is.na(Biome7$Use_intensity),], Realm == B7_LU$Realm[i] & LandUse == B7_LU$LU[i]))
}


fig6 <- ggplot(B7_LU[B7_LU$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), size = 1) +
  coord_cartesian(ylim = c(-100,150)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV', 'Cr', 'Pa' ,'PF')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         legend.position = "none",
         text = element_text(size = 25, colour = 'black'))
  #+ ggtitle("Biome 07: Tropical & Subtropical Grasslands, Savannas & Shrublands")

fig6
ggsave("FinalScriptsAndData/Figs/attempt2/Biome7_RealmLU.png")


# whole biome, just land use 

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
  ggtitle("Biome 07: Tropical & Subtropical Grasslands, Savannahs & Shrublands")

figa
ggsave(filename = 'FinalScriptsAndData/Figs/attempt2/Biome7LU.png')


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
  ggtitle("Biome 07: Tropical & Subtropical Grasslands, Savannahs & Shrublands")

figb
ggsave(filename = 'FinalScriptsAndData/Figs/attempt2/Biome7LU_UI.png')



# Abundance ---------------------------------------------------------------
Biome7_abund <- Biome7[!is.na(Biome7$Total_abundance),]

b7a0 <- StatisticalModels::GLMER(modelData = Biome7_abund[!is.na(Biome7_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LandUse', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

b7a1 <- StatisticalModels::GLMER(modelData = Biome7_abund[!is.na(Biome7_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

b7a2 <- StatisticalModels::GLMER(modelData = Biome7_abund[!is.na(Biome7_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

b7a3 <- StatisticalModels::GLMER(modelData = Biome7_abund[!is.na(Biome7_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

AIC(b7a0$model, b7a1$model, b7a2$model, b7a3$model)

R2GLMER(ba0$model)
R2GLMER(ba1$model)
R2GLMER(ba2$model)
R2GLMER(ba3$model)


# Plot Models: Abundance --------------------------------------------------
B7_LUUI3_a <- apply(RealmB7, 1, FUN = realmPredsLU_UI_a, model = b7a3$model, data = Biome7_abund[!is.na(Biome7_abund$Use_intensity),])
B7_LUUI3_a <- data.table::rbindlist(B7_LUUI3_a)
for ( i in 1:nrow(B7_LUUI3_a)) {
  B7_LUUI3_a$n[i] <- nrow(subset(Biome7, Realm == B7_LUUI3_a$Realm[i] & LU_UI_3 == B7_LUUI3_a$LU[i])) 
}


B7_LUUI3_a <- cbind(B7_LUUI3_a, LandUse, UseIntensity)



# Abundance: UseIntensity*Realm -------------------------
fig9a <- ggplot(B7_LUUI3_a[B7_LUUI3_a$n > 25,], aes(x = factor(LU, levels = level_order), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(aes(shape = UseIntensity), position = position_dodge(width = 0.6), size = 5) +
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_point(aes(shape = UseIntensity), position = position_dodge(width = 0.6), size = 5) + 
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), size = 1) +
  coord_cartesian(ylim = c(-200,400)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        text = element_text(size = 25, colour = 'black'))

fig9a
ggsave("FinalScriptsAndData/Figs/attempt2/Biome7_LUUI3_abRealm.png")



# LandUse*Realm -----------------------------------------------------------

B7_LU_a <- apply(RealmB7, 1, FUN = realmPredsLU_a, model = b7a1$model, data = Biome7_abund[!is.na(Biome7_abund$Use_intensity),])
B7_LU_a <- data.table::rbindlist(B7_LU_a)
for ( i in 1:nrow(B7_LU_a)) {
  B7_LU_a$n[i] <- nrow(subset(Biome7, Realm == B7_LU_a$Realm[i] & LandUse == B7_LU_a$LU[i])) 
}


figB7LU <- ggplot(B7_LU_a[B7_LU_a$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), size = 1) +
  coord_cartesian(ylim = c(-100,300)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV', 'Cr', 'Pa' ,'PF')) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        text = element_text(size = 25, colour = 'black')
  )

figB7LU
ggsave("FinalScriptsAndData/Figs/attempt2/Biome7_LU_aRealm.png")


write.csv(Biome7, "FinalScriptsAndData/Figs/Maps/Biome7.csv")
