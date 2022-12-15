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

#subset data to biome 4

Biome4 <- data[data$Biome=="Temperate Broadleaf & Mixed Forests",] %>%
  droplevels()

#check distribution - make a table of realm/LU_UI

table(Biome4$Realm, Biome4$LandUse)
table(Biome4$Realm, Biome4$LandUse2)
table(Biome4$Realm, Biome4$LandUse3)
table(Biome4$Realm, Biome4$LandUse4)

table(Biome4$Realm, Biome4$LU_UI)
table(Biome4$Realm, Biome4$LU_UI_3, useNA = 'always')
table(Biome4$Realm, Biome4$LU_UI_4)

#make PV Intense & Minimal the same
Biome4 <- Biome4 %>%
  mutate(LU_UI_3 = recode_factor(LU_UI_3, 
                                 "Primary Vegetation_Minimal use" = "Primary Vegetation",
                                 "Primary Vegetation_Intense use" = "Primary Vegetation"))

#remove Indo-Malay for species richness models as there are no obs in Agriculture
Biome4 <- Biome4[Biome4$Realm!='Indo-Malay',] %>%
  droplevels()

Biome4_abund <- Biome4[!is.na(Biome4$Total_abundance),]
table(Biome4_abund$Realm, Biome4_abund$LU_UI_3, useNA = 'always')


## 4.2 Selecting model structure -------------------------------------------

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
B4LUsel$deltaAIC <- B4LUsel$AIC - max(B4LUsel$AIC)
B4LUsel = arrange(B4LUsel, AIC)


# 4.3: Run Models ----------------------------------------------------------

#model selection
#using biome models function from biome1 script
B4sr<- Biomemodels1(data = Biome4[!is.na(Biome4$Use_intensity),], responseVar = 'LogRichness', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')
B4A <- Biomemodels1(data = Biome4[!is.na(Biome4$Use_intensity),], responseVar = 'LogAbund', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')

B4 <- rbind(B4sr, B4A)
B4$Dataset = 'Biome4'
B4 <- rbind(B4LUsel, B4)
write.csv(B4, "output/B4ModSel.csv", row.names = F)


TS6 = flextable(B4[,c(1:4,8,9,6,7,5,10)])
small_border = fp_border(color="black", width = 2)
tiny_border = fp_border(color="black", width = 1.5)
TS6 <- merge_v(TS6, j = ~ Dataset + Response)
TS6 <- theme_vanilla(TS6)
TS6 <- fix_border_issues(TS6)
TS6 <- hline(TS6, border = small_border, i = c(5,9))
TS6
save_as_image(TS6, 'Output/TableS6_Biome4Models.png')
save_as_docx(TS6, path = 'Output/TableS6_Biome4Models.docx')

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

#Abundance models
ba0 <- StatisticalModels::GLMER(modelData = Biome4_abund[!is.na(Biome4_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LandUse', 
                                randomStruct = "(1|SS)", REML = F)

ba1 <- StatisticalModels::GLMER(modelData = Biome4_abund[!is.na(Biome4_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                                randomStruct = "(1|SS) ", REML = F)

ba2 <- StatisticalModels::GLMER(modelData = Biome4_abund[!is.na(Biome4_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3', 
                                randomStruct = "(1|SS) ", REML = F)

ba3 <- StatisticalModels::GLMER(modelData = Biome4_abund[!is.na(Biome4_abund$Use_intensity),], responseVar = "LogAbund",
                                fitFamily = 'gaussian', fixedStruct = 'LU_UI_3*Realm', 
                                randomStruct = "(1|SS) ", REML = F)


## 4.4 Plot model ----------------------------------------------------------
RealmB4 <- data.frame(Realm = unique(Biome4$Realm))

test <- sapply(RealmB4, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})
### Fig 5a -----------------------------------------------------------------

B4_LU<- apply(RealmB4, 1, FUN = realmPredsLU, model = b1$model, data = Biome4[!is.na(Biome4$Use_intensity),])
B4_LU <- data.table::rbindlist(B4_LU)
for ( i in 1:nrow(B4_LU)) {
  B4_LU$n[i] <- nrow(subset(Biome4[!is.na(Biome4$Use_intensity),], Realm == B4_LU$Realm[i] & LandUse == B4_LU$LU[i]))
}

level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Cropland", "Pasture","Plantation forest")

figB4.a <- ggplot(B4_LU[B4_LU$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#E78AC3", "#A6D854", '#8DA0CB', "#FFD92F")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1.1,  data = B4_LU[B4_LU$n > 25 & B4_LU$upper > 100,], colour = 'black') +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  #geom_text(aes(x = LU, y = -100, label=paste('n = ', n)), 
  #          position=position_dodge(width=1.0), color = 'black'
  #) +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black')) + 
  ggtitle("a")

figB4.a


### Fig 5b -----------------------------------------------------------------

B4_LUUI3 <- apply(RealmB4, 1, FUN = realmPredsLU_UI, model = b3$model, data = Biome4)
B4_LUUI3 <- data.table::rbindlist(B4_LUUI3)
for ( i in 1:nrow(B4_LUUI3)) {
  B4_LUUI3$n[i] <- nrow(subset(Biome4, Realm == B4_LUUI3$Realm[i] & LU_UI_3 == B4_LUUI3$LU[i])) 
}
LandUse = c("Primary Vegetation", 'Secondary Vegetation', 'Secondary Vegetation', 'Agriculture', 'Agriculture')
UseIntensity <- c('PV', 'Minimal Use', 'Instense Use', 'Minimal Use', 'Instense Use')
B4_LUUI3 <- cbind(B4_LUUI3, LandUse, UseIntensity)

level_order <- c("Primary Vegetation", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use")

figB4.b <- ggplot(B4_LUUI3[B4_LUUI3$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  scale_colour_manual(values = c("#E78AC3", "#A6D854", '#8DA0CB', "#FFD92F")) +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1,  data = B4_LUUI3[B4_LUUI3$n > 25 & B4_LUUI3$upper > 130,], colour = 'black') +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use - Use Intensity', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = '') +
  theme_bw()+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black'))+
  ggtitle('b')

figB4.b


### Fig 5.c -----------------------------------------------------------------

B4_LU_a <- apply(RealmB4, 1, FUN = realmPredsLU_a, model = ba1$model, data = Biome4_abund[!is.na(Biome4_abund$Use_intensity),])
B4_LU_a <- data.table::rbindlist(B4_LU_a)
for ( i in 1:nrow(B4_LU_a)) {
  B4_LU_a$n[i] <- nrow(subset(Biome4, Realm == B4_LU_a$Realm[i] & LandUse == B4_LU_a$LU[i])) 
}


figB4.c <- ggplot(B4_LU_a[B4_LU_a$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#E78AC3", "#A6D854", '#8DA0CB', "#FFD92F")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1, angle = 45, data = B4_LU_a[B4_LU_a$n > 25 & B4_LU_a$upper > 100,], colour = 'black', position = position_dodge(width = 0.6)) +
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

figB4.c
#ggsave("FinalScriptsAndData/Figs/attempt2/Biome4_LU_aRealm.png")


# Fig 5.d -----------------------------------------------------------------

B4_LUUI3_a <- apply(RealmB4, 1, FUN = realmPredsLU_UI_a, model = ba3$model, data = Biome4_abund[!is.na(Biome4_abund$Use_intensity),])
B4_LUUI3_a <- data.table::rbindlist(B4_LUUI3_a)
for ( i in 1:nrow(B4_LUUI3_a)) {
  B4_LUUI3_a$n[i] <- nrow(subset(Biome4, Realm == B4_LUUI3_a$Realm[i] & LU_UI_3 == B4_LUUI3_a$LU[i])) 
}

B4_LUUI3_a <- cbind(B4_LUUI3_a, LandUse, UseIntensity)

figB4.d <- ggplot(B4_LUUI3_a[B4_LUUI3_a$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) +
  scale_colour_manual(values = c("#E78AC3", "#A6D854", '#8DA0CB', "#FFD92F")) + 
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1,  angle = 45, data = B4_LUUI3_a[B4_LUUI3_a$n > 25 & B4_LUUI3_a$upper > 100,], colour = 'black', position = position_dodge(width = 0.6)) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  theme_bw()+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black')) +
  ggtitle('d')

figB4.d

blank <- ggplot() + theme_void()
l <- cowplot::get_legend(figB4.a + theme(legend.position = "bottom"))

cowplot::plot_grid(l)/(figB4.a + figB4.b + figB4.c + figB4.d) + 
  plot_layout(heights = unit(c(1,20), "cm")) & theme(legend.position = "none")


write.csv(Biome4, "Figs/Maps/Biome4.csv", row.names = F)


#Not in manuscript, just for visualisation & data checking:
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
ggsave(filename = 'Figs/Biome4LU.png')


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
ggsave(filename = 'Figs/Biome4LU_UI.png')




