# 04 - GLobal Models and Figs 

#Order of code
# 1. Packages
# 2. Functions
    #runmodels()
# 3. Data
  #Load predicts database & regional biome summary
# 4. Script
  #Selecting regional biomes 
    #(removing data deficient units based on sample size thresholds)
  #Select model structure
    #Selecting appropriate fixed effects and random effects - save to output for use in Table S3. 
  #Impact of sample size
    #changing data deficiency threshold - save to output as Global_CompareSampleSize.csv. saving as Table S4. 
  #Hold-out model
    #Run hold-out models (takes approx 1 hour)
    #Plot & save results
  #Summary stats

# Packages ----------------------------------------------------------------

library(lme4) #glmer()
library(sjPlot)
library(data.table) #rbindlist()
library(ggplot2)
library(performance)
library(devtools)
library(dplyr)
library(flextable) #saving output
# load PREDICTGLMER function 
source('Code/PredictGLMERfunction.R')

install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions")
library(predictsFunctions)

install_github("timnewbold/StatisticalModels")
library(StatisticalModels)

if !dir.exists('Figs') {
  dir.create('Figs')
}
if !dir.exists('output') {
  dir.create('output')
}


# Functions ---------------------------------------------------------------

runmodels1 <- function(data, responseVar, LandUseVar = 'LandUse') {
    '''
    Takes in biodiversity dataset from predicts and runs series of glmers.
    Outputs table of AIC values and R2 values
    responseVar = SpeciesRichness, LogRichness or LogAbund
    LandUseVar = LandUse, LandUse1, LandUse2 etc.
    '''
    m <- NULL
  m[[1]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste(LandUseVar), 
                                     randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
  
  m[[2]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste0(LandUseVar, "*Biome"), 
                                     randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
  
  m[[3]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste0(LandUseVar, "*Realm"), 
                                     randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
  
  m[[4]] <- StatisticalModels::GLMER(modelData = data, responseVar = responseVar,
                                     fitFamily = 'gaussian', fixedStruct = paste0(LandUseVar, "*RB_tnc"), 
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
  
  modelresults <- data.frame('Dataset' = deparse(substitute(data)),
                             "Response" = responseVar,
                             "Fixef" = c("LandUse", "LandUse:Biome", "LandUse:Realm", "LandUse:RegionalBiome"),
                             "LandUseVar" = LandUseVar,
                             "AIC" = aic$AIC,
                             "R2Marginal" = R2$marginal, 
                             "R2Conditional" = R2$conditional,
                             "n" = nrow,
                             "n_RB" = n_distinct(data$RB_tnc))
  modelresults$deltaAIC = modelresults$AIC - min(modelresults$AIC)
  modelresults = arrange(modelresults, AIC)
  return(modelresults)
}

# Data --------------------------------------------------------------------

data <- readRDS('Data/03_PREDICTSModelData_taxa.rds')
dim(data)
ecoregs <- read.csv("Data/04_RBsummary_taxa.csv", stringsAsFactors = T)

# Code --------------------------------------------------------------------
#prep data
#remove data deficient biomes
data <- subset(data, Biome != "Tundra")
data <- subset(data, Biome != "Flooded Grasslands & Savannas")
data <- subset(data, Biome != "Mangroves")
data <- droplevels(data)

# Selecting regional biomes ----------------------------------------------------

ecoregs <- ecoregs %>%
  subset(Biome != "Tundra") %>%
  subset(Biome != "Flooded Grasslands & Savannas") %>%
  subset(Biome != "Mangroves") %>%
  subset(n_sites > 1)

#thresholds 

#entire dataset, with mangroves, tundra, flooded grasslands removed. = data

#regional biomes that have >1 observations in PV and Ag
RB_LUth.1 <- ecoregs[ecoregs$PV.th1 == 1 & ecoregs$Ag.th1 == 1,]
#regional biomes that have >5 observations in PV and Ag
RB_LUth.5 <- ecoregs[ecoregs$PV.th5 == 1 & ecoregs$Ag.th5 == 1,]
#regional biomes that have >25 observations in PV and Ag
RB_LUth.25 <- ecoregs[ecoregs$PV.th25 == 1 & ecoregs$Ag.th25 == 1,]
#regional biomes that have >50 observations in PV and Ag
RB_LUth.50 <- ecoregs[ecoregs$PV.th50 == 1 & ecoregs$Ag.th50 == 1,]

#to create a subset in data:

data_LUth.1 <- filter(data, RB_tnc %in% RB_LUth.1$RB_tnc)
data_LUth.5 <- filter(data, RB_tnc %in% RB_LUth.5$RB_tnc)
data_LUth.25 <- filter(data, RB_tnc %in% RB_LUth.25$RB_tnc)
data_LUth.50 <- filter(data, RB_tnc %in% RB_LUth.50$RB_tnc)

names(ecoregs)
# Select model structure --------------------------------------------------------------

#check random effects structure
r0 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
r1 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|my_taxa)", REML = F)
r2 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
r3 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SSB)+(1|my_taxa)", REML = F)
r4 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SSB)", REML = F)
r5 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)", REML = F)
r6 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|my_taxa)", REML = F)
r7 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc'
                             , randomStruct = 0, REML = F)

AIC(r0$model, r1$model, r2$model, r3$model,r4$model,r5$model,r6$model)
#SS & SSB & taxa is best

#check best land use structure

L0 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l1 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse2*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l2 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse3*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l3 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse4*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
l4 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogRichness",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse5*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
  

aic <- AIC(L0$model, l1$model, l2$model, l3$model, l4$model)
R2 <- NULL
R2[[1]] <- R2GLMER(L0$model)
R2[[2]] <- R2GLMER(l1$model)
R2[[3]] <- R2GLMER(l2$model)
R2[[4]] <- R2GLMER(l3$model)
R2[[5]] <- R2GLMER(l4$model)
R2 <- rbindlist(R2)

#bind model selection results
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
GlobalLUsel$deltaAIC <- GlobalLUsel$AIC - min(GlobalLUsel$AIC)

#test land use selection for abundance 
L0 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogAbund",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l1 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogAbund",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse2*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l2 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogAbund",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse3*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l3 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogAbund",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse4*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
l4 <- StatisticalModels::GLMER(modelData = data_LUth.1, responseVar = "LogAbund",
                               fitFamily = 'gaussian', fixedStruct = 'LandUse5*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

aic <- AIC(L0$model, l1$model, l2$model, l3$model, l4$model)
R2 <- NULL
R2[[1]] <- R2GLMER(L0$model)
R2[[2]] <- R2GLMER(l1$model)
R2[[3]] <- R2GLMER(l2$model)
R2[[4]] <- R2GLMER(l3$model)
R2[[5]] <- R2GLMER(l4$model)
R2 <- rbindlist(R2)

#bind model selection results
GlobalLUsel_a <- data.frame(Dataset = "Global_LUth1", 
                          Response = 'TotalAbundance', 
                          Fixef = 'LandUse*RB', 
                          LandUseVar = c("LandUse", "LandUse2", "LandUse3", "LandUse4", "LandUse5"), 
                          AIC = aic$AIC, 
                          R2Marginal = R2$marginal, 
                          R2Conditional = R2$conditional, 
                          n = nrow(data_LUth.1[(!is.na(data_LUth.1$LogAbund)),]),
                          n_RB = n_distinct(data_LUth.1$RB_tnc)
)
GlobalLUsel_a$deltaAIC <- GlobalLUsel_a$AIC - min(GlobalLUsel_a$AIC)

GlobalLUsel = rbind(GlobalLUsel, GlobalLUsel_a)

#create supp table 3
#small_border = fp_border(color="black", width = 2)

GlobalLUsel = arrange(GlobalLUsel, Response, AIC)
GlobalLUsel$R2Marginal= round(GlobalLUsel$R2Marginal, 2)
GlobalLUsel$AIC= round(GlobalLUsel$AIC, 0)
GlobalLUsel$deltaAIC= round(GlobalLUsel$deltaAIC, 2)
TS3 = flextable(GlobalLUsel[,c(2:4,8,9,6,5,10)])
TS3 <- theme_vanilla(TS3)
TS3 = set_header_labels(TS3, "Fixef" = "Fixed effects", "R2Marginal" = "Marginal R2")
TS3 <- merge_v(TS3, j = ~ Response)
TS3 <- fix_border_issues(TS3)
TS3

save_as_image(TS3, 'output/TableS3_LandUseVar.png')
save_as_docx(TS3, path = 'output/TableS3_LandUseVar.docx')

##Adding biome and realm fixed effects
m0 <- StatisticalModels::GLMER(modelData = data, responseVar = "LogRichness",
                                   fitFamily = 'gaussian', fixedStruct = 'LandUse', 
                                   randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

m1 <- StatisticalModels::GLMER(modelData = data, responseVar = "LogRichness",
                                       fitFamily = 'gaussian', fixedStruct = 'LandUse*Biome', 
                                       randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

m2 <- StatisticalModels::GLMER(modelData = data, responseVar = "LogRichness",
                                   fitFamily = 'gaussian', fixedStruct = 'LandUse*Realm', 
                                   randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

m3 <- StatisticalModels::GLMER(modelData = data, responseVar = "LogRichness",
                                   fitFamily = 'gaussian', fixedStruct = 'LandUse*RB_tnc', 
                                   randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
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
modelresults <- data.frame("Response" = c(names(m0$data[1]), names(m1$data[1]), names(m2$data[1]), names(m3$data[1])),
                          'Dataset' = 'data',
                          "AIC" = aic$AIC,
                          "R2Marginal" = R2$marginal, 
                          "R2Conditional" = R2$conditional)

#are biome and regional biome a significant variables?

m1s <- GLMERSelect(modelData = data, responseVar = "LogRichness",
                   fitFamily = 'gaussian', fixedFactors = c('LandUse', 'Biome'), fixedInteractions = ("LandUse:Biome"), 
                   randomStruct = "(1|SS)+(1|SSB)")
m1s$stats #all fixed effects are significant

m3s <- GLMERSelect(modelData = data, responseVar = "LogRichness",
                   fitFamily = 'gaussian', fixedFactors = c('LandUse', 'RB_tnc'), fixedInteractions = ("LandUse:RB_tnc"), 
                   randomStruct = "(1|SS)+(1|SSB)")

m3s$stats

#model selection suggests that Regional Biome is a stronger predictor than Biome, as it causes a stronger change in AIC.
{
  # how i did the model in lme4
  # model3<- lmer(LogRichness ~ LandUse*RB_tnc
  #               + (1|SS) + (1|SSB), 
  #               data = data, REML = F, 
  #               control = lmerControl(optimizer = 'bobyqa', 
  #                                     optCtrl = list(maxfun = 2e5)))
}
#are biome and regional biome a significant variables? - Abundance
m1as <- GLMERSelect(modelData = data, responseVar = "LogAbund",
                    fitFamily = 'gaussian', fixedFactors = c('LandUse', 'Biome'), fixedInteractions = ("LandUse:Biome"), 
                    randomStruct = "(1|SS)+(1|SSB)")
m1as$stats #all fixed effects are significant

m3as <- GLMERSelect(modelData = data, responseVar = "LogAbund",
                    fitFamily = 'gaussian', fixedFactors = c('LandUse', 'RB_tnc'), fixedInteractions = ("LandUse:RB_tnc"), 
                    randomStruct = "(1|SS)+(1|SSB)")

m3as$stats #all fixed effects significant, the interaction between LU and RB is the strongest



# Impact of sample size ---------------------------------------------------
#Using runmodels() to test the impact of changing data deficiency threshold

alldataM<- runmodels1(data = data, responseVar = 'LogRichness', LandUseVar = 'LandUse')
alldataMA <- runmodels1(data = data, responseVar = 'LogAbund', LandUseVar = 'LandUse')

#this is the dataset that should be used for everything
RB_LUth.1_S <- runmodels1(data = data_LUth.1, responseVar = 'LogRichness', LandUseVar = 'LandUse') 
RB_LUth.1_A <- runmodels1(data = data_LUth.1, responseVar = 'LogAbund', LandUseVar = 'LandUse')

#for upgrade appendix:
#Globalmodsel <- rbind(GlobalLUsel, RB_LUth.1_Sresults, RB_LUth.1_Aresults)
#write.csv(Globalmodsel, "Figs/GlobalModSel.csv", row.names = F)

RB_LUth.5_Sresults <- runmodels1(data = data_LUth.5, responseVar = 'LogRichness') 
RB_LUth.5_Aresults <- runmodels1(data = data_LUth.5, responseVar = 'LogAbund')

RB_LUth.25_Sresults <- runmodels1(data = data_LUth.25, responseVar = 'LogRichness') 
RB_LUth.25_Aresults <- runmodels1(data = data_LUth.25, responseVar = 'LogAbund')

RB_LUth.50_Sresults <- runmodels1(data = data_LUth.50, responseVar = 'LogRichness') 
RB_LUth.50_Aresults <- runmodels1(data = data_LUth.50, responseVar = 'LogAbund')

globalmods <- rbind(alldataM, alldataMA, RB_LUth.1_S, RB_LUth.1_A, RB_LUth.5_Sresults, RB_LUth.5_Aresults, RB_LUth.25_Sresults, RB_LUth.25_Aresults, RB_LUth.50_Sresults, RB_LUth.50_Aresults )
write.csv(globalmods, 'output/Global_compareSampleSize.csv', row.names = F)
#in all cases, changing sample size threshold does not change the output. 
#keep data deficiency to a limit of 25 as a good midpoint between keeping data points and reducing uncertainty. 

# Plot Table S4 ----------------------------------------------------------
globalmods$Fixef = factor(globalmods$Fixef, levels = c("LandUse", "LandUse:Realm", "LandUse:Biome", "LandUse:RegionalBiome"))

#just going to present in a table rather than figure.
globalmods = globalmods %>%
  mutate(minSampleSize = recode_factor(Dataset,
                                   'data'=0,
                                   "data_LUth.1" = 1,
                                   "data_LUth.5" = 5,
                                   "data_LUth.25" = 25,
                                   "data_LUth.50" = 50)
  )

globalmods$R2Marginal = round(globalmods$R2Marginal, 2)
globalmods$AIC = round(globalmods$AIC, 0)
globalmods$deltaAIC = round(globalmods$deltaAIC, 2)
small_border = fp_border_default(color="black", width = 2)
tiny_border = fp_border_default(color="black", width = 1.5)
TS4 = flextable(globalmods[,c(11,2,3,8,9,6,5,10)])
TS4 <- merge_v(TS4, j = ~ minSampleSize)
TS4 <- theme_vanilla(TS4)
TS4 <- fix_border_issues(TS4)
TS4 <- hline(TS4, border = small_border, i = c(8,16,24,32))
TS4 <- hline(TS4, border = tiny_border, i = c(4,12,20,28, 36))
TS4 = set_header_labels(TS4, "minSampleSize" = 'Min sample size', "Fixef"= "Fixed effects", 'R2Marginal' = 'Marginal R2')
TS4
save_as_image(TS4, 'output/TableS4_SampleSize.png')
save_as_docx(TS4, path = 'output/TableS4_SampleSize.docx')

# Hold-out Model ----------------------------------------------------------------
#Running a hold-out model. Run the model 100 times, each time removed 5% of studies.
#data set to run this on : data_LUth.25

set.seed(20210430); # Set the randon number generator seed to ensure replicatable results

N = 100
sample_size = floor(nrow(data_LUth.25)*0.90) # Make our sample size 90% of the original data
studies = unique(data_LUth.25$SSB) #1944
SSB_samplesize = floor(length(studies)*0.90) # 90% of studies 

sample_results = list() # Somewhere to put results for each sample
for (i in 1:N) {
  print(paste0('Iteration ', i))
  # For each attempt, randomly sample the data
  #data_sample = sample_n(data, sample_size)
  #alternate way of sampling where you remove 90% of studies: 
  study_sample = sample(studies, SSB_samplesize)
  data_sample = subset(data_LUth.25, SSB %in% study_sample)
  #might have to run droplevels()
  # Run models using this random sample
  print('Running species richness models')
  sample_dataM <- runmodels1(data = data_sample, responseVar = 'LogRichness', LandUseVar = "LandUse")
  sample_dataM$sample = i # Store sample number with these results
  print('Running abundance models')
  sample_dataMA <- runmodels1(data = data_sample, responseVar = 'LogAbund',LandUseVar = "LandUse")
  sample_dataMA$sample = i # Store sample number with these results
  
  sample_results[[i]] = rbind(sample_dataM, sample_dataMA) # Store in the list
}

# Convert results list to data frame
sample_results_df = rbindlist(sample_results)

# Save results
write.csv(sample_results_df, file = "output/holdout_results_LUth25_taxa.csv", row.names = F)

#need to give each model in each iteration a rank, based on deltaAIC. 
#so rank 1 will be the model with the most support

head(sample_results_df)
sample_results_df = sample_results_df %>%
  arrange(AIC) %>%
  group_by(sample, Response) %>%
  mutate(ranking = row_number())

#trying a different AIC method - best fitting model = 0
sample_results_df = sample_results_df %>%
  arrange(AIC) %>%
  group_by(sample, Response) %>%
  mutate(dAIC_2 = AIC - (min(AIC)))

head(sample_results_df)

# Plot Fig.3 --------------------------------------------------------------
sample_results_df = read.csv('output/holdout_results_LUth25_taxa.csv')
sample_results_df = sample_results_df %>%
  arrange(AIC) %>%
  group_by(sample, Response) %>%
  mutate(dAIC_2 = AIC - (min(AIC)))
# Fixing order of factors
sample_results_df$Fixef = factor(sample_results_df$Fixef, levels = c("LandUse:RegionalBiome","LandUse:Biome", "LandUse:Realm", "LandUse"))
# Boxplot of deltaAIC
sample_results_df$Response = factor(sample_results_df$Response, levels = c('LogRichness', 'LogAbund'))

labels = c('(a) Species richness', '(b) Total abundance')
names(labels) <- c("LogRichness", "LogAbund")

ggplot(sample_results_df, aes(x = Fixef, y = dAIC_2, group = Fixef)) +
  #geom_point(aes(colour = ranking), position = 'jitter') +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Response, scales = "free_y", labeller = labeller(Response = labels), ncol = 1) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 10, hjust = 0),
        strip.background = element_blank(),
        text = element_text(size = 10, colour = 'black'),
        axis.title.x = element_text(vjust = -0.8, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.text = element_text(size = 7, colour = 'black',margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(vjust = -0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        )+
  scale_x_discrete(name = 'Fixed Effects', labels = c('LU:RB', 'LU:Biome', 'LU:Realm', 'LU')) +
  ylab(expression(paste("\n",Delta, "AIC"))) +
  scale_y_continuous(limits = c(0,1200))

ggsave("Figs/Fig3_90pcstudies_25_taxa.png", dpi = 320, width = 14, height = 8, unit = 'cm')
#saving for ecography
ggsave("output/Eco_response/Fig2_90pcstudies_25_taxa.pdf", dpi = 320, width = 6.5, height = 13, unit = 'cm')

  # Summary stats -----------------------------------------------------------
#With sample size threshold 25:
#number of studies
n_studies <- n_distinct(data_LUth.25$Study_name)
n_sources <- n_distinct(data_LUth.25$Source_ID)
n_sites <- n_distinct(data_LUth.25$SSBS)
#number of sources
#number of sites

#With sample size threshold 1:
#number of studies
n_studies1 <- n_distinct(data_LUth.1$Study_name)
n_sources1 <- n_distinct(data_LUth.1$Source_ID)
n_sites1 <- n_distinct(data_LUth.1$SSBS)
n_rb1 <- n_distinct(data_LUth.1$RB_tnc)
#number of sources
#number of sites
