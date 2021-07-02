# Selecting models
# Author: Peggy Bevan
# Date: 20/01/2021

#Running models on the predicts database to ask
#Is the response of biodiversity in each region of a biome significantly different?

# Here I will run some models on the PREDICTS database to ask whether the response of biodiversity in each region of a biome is significantly different? 
#   
#   I will test this by firstly asking whether **land use, realm, biome** and **taxa** contribute to explaining variation in species richness and abundance.
# 

# For next script?:
# I will then look at the predicted values for change in species richness and abundance under each land use type and visualise how the responses change by biome. 

# Packages ----------------------------------------------------------------
library(dplyr) #mutate()
library(kableExtra) #kable() (R markdown only)
library(knitr)
library(tidyr) #drop_na()
library(lme4) #glmer()
library(sjPlot)
library(performance)

# Data --------------------------------------------------------------------

data <- read.csv("FinalScriptsAndData/Data/03_PREDICTSModelData.csv")
dim(data)
# 20609    46

# Code --------------------------------------------------------------------

###Fixed-effect & random effect selection: Does land use, realm, biome and taxa explain variation in biodiversity?

# Build the model with all possible fixed effects
#(I have tested these models with different groupings on Land Use categories and found that the orginal (5 categories) gives the lowest AIC for Richness and Abundance)


#first i need to add some new variables - a re-categorisation of land use & log species richness & abundance

data$LandUse <- relevel(data$LandUse, ref = 'Primary Vegetation')

data <- data %>%
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
  LogAbund = log(Total_abundance+1),
  LogRichness = log(Species_richness+1)
  )

# The maximal model:
#SS is source/study
#SSBS is source study site block
m1 <- lmer(LogRichness ~ LandUse*RB_tnc*CommonTaxon 
           + (1|SS) + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))

install.packages('performance')
library(performance)
install.packages('see')
install.packages('effectsize')
install.packages('emmeans')
performance::check_model(model1)

install.packages("remotes")
remotes::install_github("easystats/parameters")

# Random-effects structure ------------------------------------------------


re1<- m1
  # lmer(LogRichness ~ LandUse*RB_tnc*CommonTaxon
          # + (1|SS) + (1|SSB), 
          # data = data, REML = F, 
          # control = lmerControl(optimizer = 'bobyqa', 
          #                       optCtrl = list(maxfun = 2e5)))

#remove SSB
re2 <- lmer(LogRichness ~ LandUse*RB_tnc*CommonTaxon
           + (1|SS), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                optCtrl = list(maxfun = 2e5)))

#remove SS 
re3 <- lmer(LogRichness ~ LandUse*RB_tnc*CommonTaxon 
           + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))

#remove all
re4 <- lm(LogRichness ~ LandUse*RB_tnc*CommonTaxon,
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))

AIC(re1,re2,re3,re4)
# re1 has the lowest AIC, so all Random effects will be kept in.

# Fixed-effects structure -------------------------------------------------

## re1 is the best random-effects structure
fe1 <- lmer(LogRichness ~ LandUse*Biome*Realm 
           + (1|SS) + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))
#Only land use
fe2 <- lmer(LogRichness ~ LandUse
           + (1|SS) + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))

# Land Use + Biome
fe3 <- lmer(LogRichness ~ LandUse + Biome 
           + (1|SS) + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))

fe4 <- lmer(LogRichness ~ LandUse + Realm 
           + (1|SS) + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))

fe5 <- lmer(LogRichness ~ LandUse + Realm + Biome
           + (1|SS) + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))

fe5.1 <- lmer(LogRichness ~ LandUse + RB_tnc
             + (1|SS) + (1|SSB), 
             data = data, REML = F, 
             control = lmerControl(optimizer = 'bobyqa', 
                                   optCtrl = list(maxfun = 2e5)))


AIC(fe1, fe2, fe3, fe4, fe5, fe5.1)
# #Biome, but not realm help explain variation in the model.


# now look at their interctions
fe6 <- lmer(LogRichness ~ LandUse*Biome
           + (1|SS) + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))

fe7 <- lmer(LogRichness ~ LandUse*Realm
           + (1|SS) + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))
#LandUse*Biome*Realm = f1
fe8 <- lmer(LogRichness ~ LandUse + LandUse:Biome + LandUse:Realm
           + (1|SS) + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))

fe8.1 <- lmer(LogRichness ~ LandUse*RB_tnc
             + (1|SS) + (1|SSB), 
             data = data, REML = F, 
             control = lmerControl(optimizer = 'bobyqa', 
                                   optCtrl = list(maxfun = 2e5)))


AIC(fe1, fe2, fe6, fe7,fe8, fe8.1)

#conclusion so far: 
#realm does not explain much variation in species richness change. 
#however, when the interactions of land use with biome and realm are
# considered, we find the best fitting model.


#does adding taxa improve anything?
#best model is curerently fe8.1

fe9.1 <- lmer(LogRichness ~ LandUse + CommonTaxon
             + (1|SS) + (1|SSB), 
             data = data, REML = F, 
             control = lmerControl(optimizer = 'bobyqa', 
                                   optCtrl = list(maxfun = 2e5)))

fe9.2 <- lmer(LogRichness ~ LandUse*CommonTaxon
             + (1|SS) + (1|SSB), 
             data = data, REML = F, 
             control = lmerControl(optimizer = 'bobyqa', 
                                   optCtrl = list(maxfun = 2e5)))

fe9.3 <- lmer(LogRichness ~ LandUse + CommonTaxon + Biome
                     + (1|SS) + (1|SSB), 
                     data = data, REML = F, 
                     control = lmerControl(optimizer = 'bobyqa', 
                                           optCtrl = list(maxfun = 2e5)))

fe9.4 <- lmer(LogRichness ~ LandUse*CommonTaxon*Biome
                     + (1|SS) + (1|SSB), 
                     data = data, REML = F, 
                     control = lmerControl(optimizer = 'bobyqa', 
                                           optCtrl = list(maxfun = 2e5)))


fe9.5 <- lmer(LogRichness ~ LandUse*CommonTaxon + LandUse*Biome
                     + (1|SS) + (1|SSB), 
                     data = data, REML = F, 
                     control = lmerControl(optimizer = 'bobyqa', 
                                           optCtrl = list(maxfun = 2e5)))


fe9 <- lmer(LogRichness ~ LandUse + LandUse*Biome + LandUse*Realm + LandUse*CommonTaxon
           + (1|SS) + (1|SSB), 
           data = data, REML = F, 
           control = lmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 2e5)))

AIC(fe8.1, fe9.1, fe9.2, fe9.3, fe9.4, fe9.5, fe9)

# f10 - looking at RB/taxa interaction

fe10 <- lmer(LogRichness ~ LandUse*RB_tnc + LandUse*CommonTaxon 
            + (1|SS) + (1|SSB), 
            data = data, REML = F, 
            control = lmerControl(optimizer = 'bobyqa', 
                                  optCtrl = list(maxfun = 2e5)))

fe10.1 <- lmer(LogRichness ~ LandUse*RB_tnc*CommonTaxon
              + (1|SS) + (1|SSB), 
              data = data, REML = F, 
              control = lmerControl(optimizer = 'bobyqa', 
                                    optCtrl = list(maxfun = 2e5)))

fe10.2 <- lmer(LogRichness ~ LandUse*RB_tnc + LandUse*CommonTaxon + RB_tnc*CommonTaxon
                
              + (1|SS) + (1|SSB), 
              data = data, REML = F, 
              control = lmerControl(optimizer = 'bobyqa', 
                                    optCtrl = list(maxfun = 2e5)))

fe10.3 <- lmer(LogRichness ~ LandUse*RB_tnc + LandUse*CommonTaxon + RB_tnc:CommonTaxon:LandUse
              
              + (1|SS) + (1|SSB), 
              data = data, REML = F, 
              control = lmerControl(optimizer = 'bobyqa', 
                                    optCtrl = list(maxfun = 2e5)))



AIC(fe9.4, fe10, fe10.1, fe10.2, fe10.3)

f11 <- lmer(LogRichness ~ LandUse + RB_tnc:CommonTaxon:LandUse
              + (1|SS) + (1|SSB), 
              data = data, REML = F, 
              control = lmerControl(optimizer = 'bobyqa', 
                                    optCtrl = list(maxfun = 2e5)))

AIC(fe10.1, f11)

# L
# RB
# T
# L:RB
# L:T
# RB:T
# L:RB:T
# L*RB + L*T + RB*T + L:RB:T
#remove RB:T, and the three way interaction 


# root mean squares errors - doesn't take into account 
# mean squared error - 
# AIC
# likelihood ratio test - number of additional degrees of freedom 
# check overfitting


install.packages('AICcmodavg')
library(AICcmodavg)
## S3 method for class 'AIClmerMod'
cand.list <- list('Land Use' = fe2, 
                  'Land Use + Biome' = fe3, 
                  'Land Use + Realm' = fe4, 
                  'Land Use + Realm + Biome' = fe5, 
                  'Land Use + RB' = fe5.1,
                  'LandUse * Biome' = fe6, 
                  'LandUse * Realm' = fe7, 
                  'Land Use + Land Use:Biome + Land Use:Realm' = fe8, 
                  'LandUse*RB'= fe8.1,
                  'LandUse + Taxa' = fe9.1,
                  'LandUse*Taxa' = fe9.2, 
                  'LandUse + Taxa + Biome' = fe9.3,
                  'LandUse*Taxa*Biome' = fe9.4,
                  'LandUse*Taxa + LandUse*Biome' = fe9.5, 
                  'LandUse + LandUse*Biome + LandUse*Realm + LandUse*Taxa' = fe9,
                  'LandUse*RB + LandUse*Taxa' = fe10,
                  'LandUse*RB*Taxa' = fe10.1,
                  'LandUse*RB + LandUse*Taxa + RB*Taxa' = fe10.2,
                  'LandUse*RB + LandUse*Taxa + RB:Taxa:LandUse' = fe10.3,
                  'LandUse + RB:Taxa:LandUse' = f11)

modselect <- aictab(cand.set = cand.list, modnames = NULL,
       second.ord = TRUE, nobs = NULL, sort = TRUE)

fl1 <- flextable(modselect[,1:5])
write.csv(modselect, 'FinalScriptsAndData/Figs/05_modelselection.csv')
