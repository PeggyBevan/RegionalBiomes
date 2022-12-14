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
library(cowplot) #get_legend, grid_plot()
library(patchwork) #join plots together


# FUNCTIONS ---------------------------------------------------------------

#this function runs 4 models: Land use, land use*realm, landuseintensity, useintensity*realm
# and gives AIC and R2 values
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


# 3.2 Selecting Model Structure -------------------------------------------

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
                      Response = 'LogRichness', 
                      Fixef = 'LandUse*Realm', 
                      LandUseVar = c("LandUse", "LandUse2", "LandUse3", "LandUse4", "LandUse5"), 
                      AIC = aic$AIC, 
                      R2Marginal = R2$marginal, 
                      R2Conditional = R2$conditional, 
                      n = nrow(Biome1[!is.na(Biome1$Use_intensity),]),
                      n_RB = n_distinct(Biome1$RB_tnc)
                      )
B1LUsel$deltaAIC <-  B1LUsel$AIC - max(B1LUsel$AIC)

B1LUsel = arrange(B1LUsel, AIC)

#lowest AIC & highest R2 values is landuse 1.

# 3.3 Run Models --------------------------------------------------

B1sr<- Biomemodels1(data = Biome1[!is.na(Biome1$Use_intensity),], responseVar = 'LogRichness', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')
B1A <- Biomemodels1(data = Biome1[!is.na(Biome1$Use_intensity),], responseVar = 'LogAbund', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3')

B1 <- rbind(B1sr, B1A)
B1$Dataset = 'Biome1'
B1 <- rbind(B1LUsel, B1)
write.csv(B1, "output/B1ModSel.csv", row.names = F)

TS5 = flextable(B1[,c(1:4,8,9,6,7,5,10)])
small_border = fp_border(color="black", width = 2)
tiny_border = fp_border(color="black", width = 1.5)
TS5 <- merge_v(TS5, j = ~ Dataset + Response)
TS5 <- theme_vanilla(TS5)
TS5 <- fix_border_issues(TS5)
TS5 <- hline(TS5, border = small_border, i = c(5,9))
TS5
save_as_image(TS5, 'Output/TableS5_Biome1Models.png')
save_as_docx(TS4, path = 'Output/TableS5_Biome1Models.docx')

# 3.4 Plot models ----------------------------------------------------------

### Overall Biome  ----------------------------------------------------------


#factors and order of factors in this model
LandUseFact <- c("Primary Vegetation", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use")

# overall biome - Land Use
#FOR CHECKING - FIGURE NOT USED IN MANUSCRIPT
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

#plot
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
ggsave(filename = 'Figs/Biome1LU.png')


#Overall biome Land Use * Use intensity
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
  geom_point(position = position_dodge(width = 0.5), size = 3.5) + 
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
ggsave(filename = 'Figs/Biome1LU_UI.png')


## By regional biome -------------------------------------------------------

### Fig4A -------------------------------------------------------------------



#model 

b1 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),], 
                               responseVar = "LogRichness",
                               fitFamily = 'gaussian', 
                               fixedStruct = 'LandUse*Realm',
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

#list realms in this biome
RealmB1 <- data.frame(Realm = unique(Biome1$Realm))

#check its the right realms
test <- sapply(RealmB1, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})


B1_LU<- apply(RealmB1, 1, FUN = realmPredsLU, model = b1$model, data = Biome1[!is.na(Biome1$Use_intensity),])
B1_LU <- data.table::rbindlist(B1_LU)
for ( i in 1:nrow(B1_LU)) {
  B1_LU$n[i] <- nrow(subset(Biome1[!is.na(Biome1$Use_intensity),], Realm == B1_LU$Realm[i] & LandUse == B1_LU$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Cropland", "Pasture","Plantation forest")


figB1.a <- ggplot(B1_LU[B1_LU$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr')) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  #geom_text(aes(y = lower ), position="identity", color = 'black', size = 3, angle = 45) +       
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    #legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black'))+
  ggtitle("a")

figB1.a
ggsave("Figs/attempt2/Biome1_RealmLU.png")

#colours = Afrotropic = '#DD8D29', Indo-Malay = '#E2D200', Neotropic = '#46ACC8'

### Fig4.B ------------------------------------------------------------------

b3 <- StatisticalModels::GLMER(modelData = Biome1[!is.na(Biome1$Use_intensity),],
                               responseVar = "LogRichness",
                               fitFamily = 'gaussian', 
                               fixedStruct = 'LU_UI_3*Realm', 
                               randomStruct = "(1|SS)", REML = F)

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


figB1.b <- ggplot(B1_LUUI3[B1_LUUI3$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.5), size = 5) + 
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.5), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.5), size = 2) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = '') +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black'))+
  ggtitle('b')

figB1.b


### Fig 4C ------------------------------------------------------------------

b1a1 <- StatisticalModels::GLMER(modelData = Biome1_abund[!is.na(Biome1_abund$Use_intensity),], responseVar = "LogAbund",
                                 fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                                 randomStruct = "(1|SS)", REML = F)

B1_LU_a <- apply(RealmB1, 1, FUN = realmPredsLU_a, model = b1a1$model, data = Biome1)
B1_LU_a <- data.table::rbindlist(B1_LU_a)
for ( i in 1:nrow(B1_LU_a)) {
  B1_LU_a$n[i] <- nrow(subset(Biome1, Realm == B1_LU_a$Realm[i] & LandUse == B1_LU_a$LU[i])) 
}


figB1.c <- ggplot(B1_LU_a[B1_LU_a$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =115, ),  data = B1_LU_a[B1_LU_a$n > 25 & B1_LU_a$upper > 110,], colour = 'black', position = position_dodge(width = 0.6) ) +
  coord_cartesian(ylim = c(-110,110)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr')) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  theme_bw()+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black')
  ) +
  ggtitle('c')

figB1.c


### Fig 4.d -----------------------------------------------------------------

b1a3 <- StatisticalModels::GLMER(modelData = Biome1_abund[!is.na(Biome1_abund$Use_intensity),], responseVar = "LogAbund",
                                 fitFamily = 'gaussian', fixedStruct = 'LU_UI_3*Realm', 
                                 randomStruct = "(1|SS)", REML = F)

B1_LUUI3_a <- apply(RealmB1, 1, FUN = realmPredsLU_UI_a, model = b1a3$model, data = Biome1_abund[!is.na(Biome1_abund$Use_intensity),])

B1_LUUI3_a <- data.table::rbindlist(B1_LUUI3_a)
for ( i in 1:nrow(B1_LUUI3_a)) {
  B1_LUUI3_a$n[i] <- nrow(subset(Biome1, Realm == B1_LUUI3_a$Realm[i] & LU_UI_3 == B1_LUUI3_a$LU[i])) 
}


B1_LUUI3_a <- cbind(B1_LUUI3_a, LandUse, UseIntensity)

figB1.d <- ggplot(B1_LUUI3_a[B1_LUUI3_a$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  #geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =115, ),  data = B1_LUUI3_a[B1_LUUI3_a$n > 25 & B1_LUUI3_a$upper > 100,], colour = 'black', position = position_dodge(width = 0.6) ) +
  coord_cartesian(ylim = c(-110,110)) +
  scale_x_discrete(name = 'Land Use - Use Intensity', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = '') +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    #legend.position = 'top',
    text = element_text(size = 20, colour = 'black',),
    legend.title = element_blank())+
  ggtitle('d')

figB1.d



l <- cowplot::get_legend(figB1.a + theme(legend.position = "bottom"))

cowplot::plot_grid(l)/(figB1.a + figB1.b + figB1.c + figB1.d) + 
  plot_layout(heights = unit(c(1,20), "cm")) & theme(legend.position = "none")
#plot_annotation(title = "main title") 

#adding figures together
figB1.a + figB1.b + figB1.c +figB1.d + (plot_layout(guides = 'collect')) & theme(legend.position = 'top')



write.csv(Biome1, "Figs/Maps/Biome1.csv", row.names = F)
