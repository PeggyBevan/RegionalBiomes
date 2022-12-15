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
library(data.table)
library(cowplot)
library(patchwork)


# 2. FUNCTIONS ------------------------------------------------------------
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
  modelresults$deltaAIC = modelresults$AIC - max(modelresults$AIC)
  modelresults = arrange(modelresults, AIC)
  return(modelresults)
}

#creating function to make predicted values from model - specific to land use variable
realmPredsLU_UI <- function(realm, model, data) {
  nd <- data.frame(LU_UI_3=factor(
    c("Primary Vegetation", 
      "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", 
      "Agriculture_Minimal use", "Agriculture_Intense use"),
    levels=levels(data$LU_UI_3)))
  nd$Realm <- factor(realm, levels = levels(data$Realm))
  nd$LogRichness = 0
  #95% confidence intervals
  preds <- StatisticalModels::PredictGLMER(model = model, data = nd,
                                           se.fit = TRUE, seMultiplier = 1.96, randEffs = F)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  #75% confidence intervals
  preds2 <- StatisticalModels::PredictGLMER(model = model, data = nd,
                                            se.fit = TRUE, seMultiplier = 1.15, randEffs = F)
  
  preds2$yplus75 <- ((exp(preds2$yplus)/exp(preds2$y[1]))*100)-100
  preds2$yminus75 <- ((exp(preds2$yminus)/exp(preds2$y[1]))*100)-100
  preds2$y75 <- ((exp(preds2$y)/exp(preds2$y[1]))*100)-100
  
  
  temp.df <- data.frame(Realm = nd$Realm, LU=nd$LU_UI_3, y=preds$y,lower=preds$yminus,upper=preds$yplus, y75 = preds2$y75, lower75=preds2$yminus75, upper75=preds2$yplus75)
  
}

realmPredsLU <- function(realm, model, data) {
  nd <- data.frame(LandUse=factor(
    c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland"),
    levels=levels(data$LandUse)))
  nd$Realm <- factor(realm, levels = levels(data$Realm))
  nd$LogRichness = 0
  #95% confidence intervals
  preds <- PredictGLMER(model = model, data = nd,
                        se.fit = TRUE, seMultiplier = 1.96)
  
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  #75% confidence intervals
  preds2 <- StatisticalModels::PredictGLMER(model = model, data = nd,
                                            se.fit = TRUE, seMultiplier = 1.15, randEffs = F)
  
  preds2$yplus75 <- ((exp(preds2$yplus)/exp(preds2$y[1]))*100)-100
  preds2$yminus75 <- ((exp(preds2$yminus)/exp(preds2$y[1]))*100)-100
  preds2$y75 <- ((exp(preds2$y)/exp(preds2$y[1]))*100)-100
  
  
  temp.df <- data.frame(Realm = nd$Realm, LU=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus, y75 = preds2$y75, lower75=preds2$yminus75, upper75=preds2$yplus75)
  
}

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
  
  #75% confidence
  #75% confidence intervals
  preds2 <- StatisticalModels::PredictGLMER(model = model, data = nd,
                                            se.fit = TRUE, seMultiplier = 1.15, randEffs = F)
  
  preds2$yplus75 <- ((exp(preds2$yplus)/exp(preds2$y[1]))*100)-100
  preds2$yminus75 <- ((exp(preds2$yminus)/exp(preds2$y[1]))*100)-100
  preds2$y75 <- ((exp(preds2$y)/exp(preds2$y[1]))*100)-100
  
  
  temp.df <- data.frame(Realm = nd$Realm, LU=nd$LU_UI_3, y=preds$y,lower=preds$yminus,upper=preds$yplus, y75 = preds2$y75, lower75=preds2$yminus75, upper75=preds2$yplus75)
  
}

realmPredsLU_a <- function(realm, model, data) {
  nd <- data.frame(LandUse=factor(
    c("Primary Vegetation", "Secondary Vegetation", "Plantation forest",  "Pasture","Cropland"),
    levels=levels(data$LandUse)))
  nd$Realm <- factor(realm, levels = levels(data$Realm))
  nd$LogAbund = 0
  preds <- PredictGLMER(model = model, data = nd,
                        se.fit = TRUE, seMultiplier = 1.96)  
  #95% CI
  preds$yplus <- ((exp(preds$yplus)/exp(preds$y[1]))*100)-100
  preds$yminus <- ((exp(preds$yminus)/exp(preds$y[1]))*100)-100
  preds$y <- ((exp(preds$y)/exp(preds$y[1]))*100)-100
  
  #75% confidence intervals
  preds2 <- StatisticalModels::PredictGLMER(model = model, data = nd,
                                            se.fit = TRUE, seMultiplier = 1.15, randEffs = F)
  
  preds2$yplus75 <- ((exp(preds2$yplus)/exp(preds2$y[1]))*100)-100
  preds2$yminus75 <- ((exp(preds2$yminus)/exp(preds2$y[1]))*100)-100
  preds2$y75 <- ((exp(preds2$y)/exp(preds2$y[1]))*100)-100
  
  
  temp.df <- data.frame(Realm = nd$Realm, LU=nd$LandUse, y=preds$y,lower=preds$yminus,upper=preds$yplus, y75 = preds2$y75, lower75=preds2$yminus75, upper75=preds2$yplus75)
  
}

# 3. DATA -----------------------------------------------------------------
data <- readRDS('Data/03_PREDICTSModelData.rds')
dim(data)
# 4. SCRIPT ---------------------------------------------------------------
## 4.1 Prep and subset data ------------------------------------------------
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


## 4.2 Selecting model structure -------------------------------------------
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
B7LUsel$deltaAIC <- B7LUsel$AIC - max(B7LUsel$AIC)
B7LUsel = arrange(B7LUsel, AIC)


## 4.3 Run Models ----------------------------------------------------------

#model selection
#using biome models function from biome1 script
B7sr<- Biomemodels1(data = Biome7[!is.na(Biome7$Use_intensity),], responseVar = 'LogRichness', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')
B7A <- Biomemodels1(data = Biome7[!is.na(Biome7$Use_intensity),], responseVar = 'LogAbund', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')

B7 <- rbind(B7sr, B7A)
B7$Dataset = 'Biome7'
B7 <- rbind(B7LUsel, B7)
write.csv(B7, "output/B7ModSel.csv", row.names = F)

TS7 = flextable(B7[,c(1:4,8,9,6,7,5,10)])
small_border = fp_border(color="black", width = 2)
tiny_border = fp_border(color="black", width = 1.5)
TS7 <- merge_v(TS7, j = ~ Dataset + Response)
TS7 <- theme_vanilla(TS7)
TS7 <- fix_border_issues(TS7)
TS7 <- hline(TS7, border = small_border, i = c(5,9))
TS7
save_as_image(TS7, 'Output/TableS7_Biome7Models.png')
save_as_docx(TS7, path = 'Output/TableS7_Biome7Models.docx')

b70 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "LogRichness",
                                fitFamily = 'gaussian', fixedStruct = 'LandUse', 
                                randomStruct = "(1|SS) ", REML = F)


b71 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "LogRichness",
                                fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                                randomStruct = "(1|SS) ", REML = F)

b72 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "LogRichness",
                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3', 
                                randomStruct = "(1|SS) ", REML = F)

b73 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "LogRichness",
                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3*Realm', 
                                randomStruct = "(1|SS) ", REML = F)


AIC(b70$model, b71$model, b72$model, b73$model)

#Abundance
b7a0 <- StatisticalModels::GLMER(modelData = Biome7_abund[!is.na(Biome7_abund$Use_intensity),], responseVar = "LogAbund",
                                 fitFamily = 'gaussian', fixedStruct = 'LandUse', 
                                 randomStruct = "(1|SS) ", REML = F)

b7a1 <- StatisticalModels::GLMER(modelData = Biome7_abund[!is.na(Biome7_abund$Use_intensity),], responseVar = "LogAbund",
                                 fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                                 randomStruct = "(1|SS) ", REML = F)

b7a2 <- StatisticalModels::GLMER(modelData = Biome7_abund[!is.na(Biome7_abund$Use_intensity),], responseVar = "LogAbund",
                                 fitFamily = 'gaussian', fixedStruct = 'LU_UI_3', 
                                 randomStruct = "(1|SS) ", REML = F)

b7a3 <- StatisticalModels::GLMER(modelData = Biome7_abund[!is.na(Biome7_abund$Use_intensity),], responseVar = "LogAbund",
                                 fitFamily = 'gaussian', fixedStruct = 'LU_UI_3*Realm', 
                                 randomStruct = "(1|SS) ", REML = F)

AIC(b7a0$model, b7a1$model, b7a2$model, b7a3$model)

## 4.3 Plot Models ----------------------------------------------------------

RealmB7 <- data.frame(Realm = unique(Biome7$Realm))

test <- sapply(RealmB7, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})

# Fig 6.a -----------------------------------------------------------------

B7_LU<- apply(RealmB7, 1, FUN = realmPredsLU, model = b71$model, data = Biome7[!is.na(Biome7$Use_intensity),])
B7_LU <- data.table::rbindlist(B7_LU)
for ( i in 1:nrow(B7_LU)) {
  B7_LU$n[i] <- nrow(subset(Biome7[!is.na(Biome7$Use_intensity),], Realm == B7_LU$Realm[i] & LandUse == B7_LU$LU[i]))
}


figB7.a <- ggplot(B7_LU[B7_LU$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  #geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1.1,  data = B7_LUUI3[B7_LUUI3$n > 25 & B7_LUUI3$upper > 100,], colour = 'black') +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black')) + 
  ggtitle("a")

figB7.a

# Fig 6.b -----------------------------------------------------------------


B7_LUUI3 <- apply(RealmB7, 1, FUN = realmPredsLU_UI, model = b73$model, data = Biome7[!is.na(Biome7$Use_intensity),])
B7_LUUI3 <- data.table::rbindlist(B7_LUUI3)
for ( i in 1:nrow(B7_LUUI3)) {
  B7_LUUI3$n[i] <- nrow(subset(Biome7, Realm == B7_LUUI3$Realm[i] & LU_UI_3 == B7_LUUI3$LU[i])) 
}
B7_LUUI3 <- cbind(B7_LUUI3, LandUse, UseIntensity)


figB7.b <- ggplot(B7_LUUI3[B7_LUUI3$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  #geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1.1,  data = B7_LUUI3[B7_LUUI3$n > 25 & B7_LUUI3$upper > 100,], colour = 'black') +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = '') +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black')) +
  ggtitle('b')

figB7.b



# Fig 6.c -----------------------------------------------------------------

B7_LU_a <- apply(RealmB7, 1, FUN = realmPredsLU_a, model = b7a1$model, data = Biome7_abund[!is.na(Biome7_abund$Use_intensity),])
B7_LU_a <- data.table::rbindlist(B7_LU_a)
for ( i in 1:nrow(B7_LU_a)) {
  B7_LU_a$n[i] <- nrow(subset(Biome7, Realm == B7_LU_a$Realm[i] & LandUse == B7_LU_a$LU[i])) 
}


figB7.c <- ggplot(B7_LU_a[B7_LU_a$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100), hjust = 1.1, angle = 45,  data = B7_LU_a[B7_LU_a$n > 25 & B7_LU_a$upper > 100,], colour = 'black', position = position_dodge(width = 0.6)) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr')) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  theme_bw()+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black')) + 
  ggtitle('c')

figB7.c

# Fig 6.d -----------------------------------------------------------------

B7_LUUI3_a <- apply(RealmB7, 1, FUN = realmPredsLU_UI_a, model = b7a3$model, data = Biome7_abund[!is.na(Biome7_abund$Use_intensity),])
B7_LUUI3_a <- data.table::rbindlist(B7_LUUI3_a)
for ( i in 1:nrow(B7_LUUI3_a)) {
  B7_LUUI3_a$n[i] <- nrow(subset(Biome7, Realm == B7_LUUI3_a$Realm[i] & LU_UI_3 == B7_LUUI3_a$LU[i])) 
}

B7_LUUI3_a <- cbind(B7_LUUI3_a, LandUse, UseIntensity)

figB7.d <- ggplot(B7_LUUI3_a[B7_LUUI3_a$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) +
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100), hjust = 1.1, angle = 45,  data = B7_LUUI3_a[B7_LUUI3_a$n > 25 & B7_LUUI3_a$upper > 100,], position = position_dodge(width = 0.6), colour = 'black') +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use - Use Intensity', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = '') +
  theme_bw()+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black')) +
  ggtitle('d')

figB7.d

l <- cowplot::get_legend(figB7.a + theme(legend.position = "bottom"))

cowplot::plot_grid(l)/(figB7.a + figB7.b + figB7.c + figB7.d) + 
  plot_layout(heights = unit(c(1,20), "cm")) & theme(legend.position = "none")

write.csv(Biome7, "Figs/Maps/Biome7.csv")
