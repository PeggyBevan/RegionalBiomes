#############################
#Biome case study functions
#############################

#Functions used in calculating model outputs for the biome case studies in scripts 5 - 7. 

# 1. FUNCTIONS ---------------------------------------------------------------

#this function runs 4 models: Land use, land use*realm, landuseintensity, useintensity*realm
# and gives AIC and R2 values
Biomemodels1 <- function(data, responseVar, LandUseVar, UseIntensityVar, fitfamily) {
  m <- NULL
  m[[1]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'poisson', fixedStruct = paste(LandUseVar), 
                                     randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
  
  m[[2]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = fitfamily, fixedStruct = paste0(LandUseVar, "*Realm"), 
                                     randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
  
  m[[3]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'poisson', fixedStruct = paste0(UseIntensityVar), 
                                     randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
  
  m[[4]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'poisson', fixedStruct = paste0(UseIntensityVar, "*Realm"), 
                                     randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
  
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
  modelresults$deltaAIC = modelresults$AIC - min(modelresults$AIC)
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
  nd$Species_richness = 0
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
  #nd$Species_richness = 0
  nd$Species_richness = 0
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