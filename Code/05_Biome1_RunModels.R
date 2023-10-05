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
library(officer)
library(ggplot2)
library(devtools)
library(StatisticalModels)
library(data.table) #rbindlist() #function for getting predicted values from all models
library(cowplot) #get_legend, grid_plot()
library(patchwork) #join plots together
library(emmeans) #post-hoc testing
source('code/BiomeModelFunctions.R')

# 2. DATA -----------------------------------------------------------------
data <- readRDS('Data/03_PREDICTSModelData_taxa.rds')
 dim(data)
 

# 4. SCRIPT ---------------------------------------------------------------
# 4.1 Prep and subset data ------------------------------------------------

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
table(Biome1_abund$Realm, Biome1_abund$LandUse, useNA = 'always')

#remove rows where Use Intensity is NA to keep dataset the same across models
Biome1 = Biome1[!is.na(Biome1$Use_intensity),]
# 4.2 Selecting Model Structure -------------------------------------------

#test random effect structure 
#check random effects structure

r0 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
r1 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|my_taxa)", REML = F)
r2 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
r3 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SSB)+(1|my_taxa)", REML = F)
r4 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SSB)", REML = F)
r5 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)", REML = F)
r6 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|my_taxa)", REML = F)
r7 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc'
                               , randomStruct = 0, REML = F)
AIC(r0$model, r1$model, r2$model, r3$model,r4$model,r5$model,r6$model)
#SS & SSB & my_taxa is best

#check best land use structure

L10 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l11 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse2*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l12 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse3*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l13 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse4*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
l14 <- StatisticalModels::GLMER(modelData = Biome1, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse5*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)



aic <- AIC(L10$model, l11$model, l12$model, l13$model, l14$model)
R2 <- NULL
R2[[1]] <- R2GLMER(L10$model)
R2[[2]] <- R2GLMER(l11$model)
R2[[3]] <- R2GLMER(l12$model)
R2[[4]] <- R2GLMER(l13$model)
R2[[5]] <- R2GLMER(l14$model)
R2 <- rbindlist(R2)

B1LUsel <- data.frame(Dataset = "Biome1", 
                      Response = 'Species_richness', 
                      Fixef = 'LandUse*Realm', 
                      LandUseVar = c("LandUse", "LandUse2", "LandUse3", "LandUse4", "LandUse5"), 
                      AIC = aic$AIC, 
                      R2Marginal = R2$marginal, 
                      R2Conditional = R2$conditional, 
                      n = nrow(Biome1[!is.na(Biome1$Use_intensity),]),
                      n_RB = n_distinct(Biome1$RB_tnc)
                      )
B1LUsel$deltaAIC <-  B1LUsel$AIC - min(B1LUsel$AIC)

B1LUsel = arrange(B1LUsel, AIC)

#lowest AIC & highest R2 values is landuse 1.

# 4.3 Run Models --------------------------------------------------

B1sr<- Biomemodels1(data = Biome1, responseVar = 'Species_richness', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3', fitfamily = 'poisson')
B1A <- Biomemodels1(data = Biome1, responseVar = 'LogAbund', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3', fitfamily = 'gaussian')

B1 <- rbind(B1sr, B1A)
B1$Dataset = 'Biome1'
B1 <- rbind(B1LUsel, B1)
write.csv(B1, "output/B1ModSel.csv", row.names = F)

names(B1)[3] = 'Fixed Effects'
names(B1)[6] = 'Marginal R2'
names(B1)[10] = 'dAIC'

B1$`Marginal R2`= round(B1$`Marginal R2`, 2)
B1$AIC = round(B1$AIC, 0)
B1$dAIC = round(B1$dAIC,2)

TS5 = flextable(B1[,c(2:4,8,9,6,5,10)])
small_border = fp_border(color="black", width = 2)
tiny_border = fp_border(color="black", width = 1.5)
TS5 <- merge_v(TS5, j = ~ Response)
TS5 <- theme_vanilla(TS5)
TS5 <- fix_border_issues(TS5)
TS5 <- hline(TS5, border = small_border, i = c(5,9))
TS5
save_as_image(TS5, 'Output/TableS5_Biome1Models.png')
save_as_docx(TS5, path = 'Output/TableS5_Biome1Models.docx')

# 4.4 Plot models ----------------------------------------------------------

### Overall Biome  ----------------------------------------------------------

#factors and order of factors in this model
LandUseFact <- c("Primary Vegetation", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use")

# overall biome - Land Use
#FOR CHECKING - FIGURE NOT USED IN MANUSCRIPT
#model:
b0 <- StatisticalModels::GLMER(modelData = Biome1, 
                               responseVar = "Species_richness",
                               fitFamily = 'poisson', 
                               fixedStruct = 'LandUse', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

nd <- data.frame(LandUse = factor(levels(b0$data$LandUse),
                                  levels = levels(b0$data$LandUse)),
                 Species_richness = 0)

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
                               responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LU_UI_3', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

nd <- data.frame(LU_UI_3 = factor(levels(b2$data$LU_UI_3),
                                  levels = levels(b2$data$LU_UI_3)),
                 Species_richness = 0)
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
Biome1save = Biome1
Biome1 = filter(Biome1, my_taxa != 'Fungi')

b1 <- StatisticalModels::GLMER(modelData = Biome1, 
                               responseVar = "Species_richness",
                               fitFamily = 'poisson', 
                               fixedStruct = 'LandUse*Realm',
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
#n = 5536
#list realms in this biome
RealmB1 <- data.frame(Realm = unique(Biome1$Realm))

#check its the right realms
test <- sapply(RealmB1, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})


B1_LU<- apply(RealmB1, 1, FUN = realmPredsLU, model = b1$model, data = Biome1)
B1_LU <- data.table::rbindlist(B1_LU)
for ( i in 1:nrow(B1_LU)) {
  B1_LU$n[i] <- nrow(subset(Biome1, Realm == B1_LU$Realm[i] & LandUse == B1_LU$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")


figB1.a <- ggplot(B1_LU[B1_LU$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), linewidth = 1) +
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
ggsave("Figs/attempt2/Biome1_RealmLU_taxa.png")

#colours = Afrotropic = '#DD8D29', Indo-Malay = '#E2D200', Neotropic = '#46ACC8'

### Fig4.B ------------------------------------------------------------------

b3 <- StatisticalModels::GLMER(modelData = Biome1,
                               responseVar = "Species_richness",
                               fitFamily = 'poisson', 
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


figB1.b <- ggplot(B1_LUUI3[B1_LUUI3$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.5), size = 5) + 
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.5), linewidth = 1) +
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
                                 randomStruct = "(1|SS)+(1|SSB)", REML = F)

B1_LU_a <- apply(RealmB1, 1, FUN = realmPredsLU_a, model = b1a1$model, data = Biome1_abund[!is.na(Biome1_abund$Use_intensity),])
B1_LU_a <- data.table::rbindlist(B1_LU_a)
for ( i in 1:nrow(B1_LU_a)) {
  B1_LU_a$n[i] <- nrow(subset(Biome1_abund[!is.na(Biome1_abund$Use_intensity),], Realm == B1_LU_a$Realm[i] & LandUse == B1_LU_a$LU[i]))
}



figB1.c <- ggplot(B1_LU_a[B1_LU_a$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ),  data = B1_LU_a[B1_LU_a$n > 25 & B1_LU_a$upper > 110,], colour = 'black', position = position_dodge(width = 0.6) ) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr')) +
  scale_y_continuous(name = 'Change in Total Abundance (%)') +
  theme_bw()+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black')
  ) +
  ggtitle('a')

figB1.c


### Fig 4.d -----------------------------------------------------------------

b1a3 <- StatisticalModels::GLMER(modelData = Biome1_abund[!is.na(Biome1_abund$Use_intensity),], responseVar = "LogAbund",
                                 fitFamily = 'gaussian', fixedStruct = 'LU_UI_3*Realm', 
                                 randomStruct = "(1|SS)+(1|SSB)", REML = F)

B1_LUUI3_a <- apply(RealmB1, 1, FUN = realmPredsLU_UI_a, model = b1a3$model, data = Biome1_abund[!is.na(Biome1_abund$Use_intensity),])

B1_LUUI3_a <- data.table::rbindlist(B1_LUUI3_a)
for ( i in 1:nrow(B1_LUUI3_a)) {
  B1_LUUI3_a$n[i] <- nrow(subset(Biome1_abund[!is.na(Biome1_abund$Use_intensity),], Realm == B1_LUUI3_a$Realm[i] & LU_UI_3 == B1_LUUI3_a$LU[i])) 
}


B1_LUUI3_a <- cbind(B1_LUUI3_a, LandUse, UseIntensity)

figB1.d <- ggplot(B1_LUUI3_a[B1_LUUI3_a$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  #geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =115, ),  data = B1_LUUI3_a[B1_LUUI3_a$n > 25 & B1_LUUI3_a$upper > 100,], colour = 'black', position = position_dodge(width = 0.6) ) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use - Use Intensity', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = '') +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    #legend.position = 'top',
    text = element_text(size = 20, colour = 'black',),
    legend.title = element_blank())+
  ggtitle('b')

figB1.d



# Fig4 --------------------------------------------------------------------

l <- cowplot::get_legend(figB1.a + theme(legend.position = "bottom"))

#plotting species richness only, abundance for supplementary:
(cowplot::plot_grid(l) + plot_spacer())/(figB1.a + figB1.b) + 
  plot_layout(heights = unit(c(1,10), "cm")) & theme(legend.position = "none")
#plot_annotation(title = "main title") 
ggsave("Figs/Fig4_SR.png", width = 13.7, height = 12)
#make sure graphics box is tall & wide enough to stop text overlapping.
#Saving 13.7 x 12 in image

#abundance plots for supplementary:
(cowplot::plot_grid(l) + plot_spacer())/(figB1.c + figB1.d) + 
  plot_layout(heights = unit(c(1,10), "cm")) & theme(legend.position = "none")
#plot_annotation(title = "main title") 
ggsave("Figs/FigS1_A.png", width = 13.7, height = 12)

#All 4 plots:
# cowplot::plot_grid(l)/(figB1.a + figB1.b + figB1.c + figB1.d) + 
#   plot_layout(heights = unit(c(1,20), "cm")) & theme(legend.position = "none")
# #plot_annotation(title = "main title") 
# ggsave("Figs/Fig4.png", width = 13.7, height = 12)
# #make sure graphics box is tall & wide enough to stop text overlapping.
# #Saving 13.7 x 12 in image

write.csv(Biome1, "Figs/Maps/Biome1.csv", row.names = F)


# Taxon case study --------------------------------------------------------

unique(Biome1$my_taxa)
table(Biome1$LandUse, Biome1$my_taxa)

Biome1 %>%
  group_by(my_taxa) %>%
  summarise(min = min(Species_richness), max = max(Species_richness))

#just vertebrates
Biome1_v = filter(Biome1, my_taxa == 'Vertebrate')

b1v <- StatisticalModels::GLMER(modelData = Biome1_v, 
                               responseVar = "Species_richness",
                               fitFamily = 'poisson', 
                               fixedStruct = 'LandUse*Realm',
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
#list realms in this biome
RealmB1v <- data.frame(Realm = unique(Biome1_v$Realm))

#check its the right realms
test <- sapply(RealmB1v, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})


B1_LUv<- apply(RealmB1v, 1, FUN = realmPredsLU, model = b1v$model, data = Biome1_v)
B1_LUv <- data.table::rbindlist(B1_LUv)
for ( i in 1:nrow(B1_LUv)) {
  B1_LUv$n[i] <- nrow(subset(Biome1_v, Realm == B1_LUv$Realm[i] & LandUse == B1_LUv$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")


figB1v.a <- ggplot(B1_LUv[B1_LUv$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.7), linewidth = 0.7) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), size = 1.5) +
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

figB1v.a

#just invertebrates
Biome1_i = filter(Biome1, my_taxa == 'Invertebrate')

b1i <- StatisticalModels::GLMER(modelData = Biome1_i, 
                                responseVar = "Species_richness",
                                fitFamily = 'poisson', 
                                fixedStruct = 'LandUse*Realm',
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
#list realms in this biome
RealmB1i <- data.frame(Realm = unique(Biome1_i$Realm))

#check its the right realms
test <- sapply(RealmB1i, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})

B1_LUi<- apply(RealmB1i, 1, FUN = realmPredsLU, model = b1i$model, data = Biome1_i)
B1_LUi <- data.table::rbindlist(B1_LUi)
for ( i in 1:nrow(B1_LUi)) {
  B1_LUi$n[i] <- nrow(subset(Biome1_i, Realm == B1_LUi$Realm[i] & LandUse == B1_LUi$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")


figB1i.a <- ggplot(B1_LUi[B1_LUi$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.7), linewidth = 0.7) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), size = 1.5) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr')) +
  scale_y_continuous(name = '') +
  #geom_text(aes(y = lower ), position="identity", color = 'black', size = 3, angle = 45) +       
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    #legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black'))+
  ggtitle("b")

figB1i.a
figB1v.a

#just plants
Biome1_p = filter(Biome1, my_taxa == 'Plant')

b1p <- StatisticalModels::GLMER(modelData = Biome1_p, 
                                responseVar = "Species_richness",
                                fitFamily = 'poisson', 
                                fixedStruct = 'LandUse*Realm',
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
#n = 5536
#list realms in this biome
RealmB1p <- data.frame(Realm = unique(Biome1_p$Realm))

#check its the right realms
test <- sapply(RealmB1p, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})


B1_LUp<- apply(RealmB1p, 1, FUN = realmPredsLU, model = b1p$model, data = Biome1_p)
B1_LUp <- data.table::rbindlist(B1_LUp)
for ( i in 1:nrow(B1_LUp)) {
  B1_LUp$n[i] <- nrow(subset(Biome1_p, Realm == B1_LUp$Realm[i] & LandUse == B1_LUp$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")


figb1p.a <- ggplot(B1_LUp[B1_LUp$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.7), linewidth = 0.7) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), linewidth = 1.5) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr')) +
  scale_y_continuous(name = '') +
  #geom_text(aes(y = lower ), position="identity", color = 'black', size = 3, angle = 45) +       
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    #legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black'))+
  ggtitle("c")

figb1p.a
figB1i.a
figB1v.a

l <- cowplot::get_legend(figB1v.a + theme(legend.position = "bottom"))

(cowplot::plot_grid(l)+plot_spacer())/(figB1v.a + figB1i.a + figb1p.a) + 
  plot_layout(heights = unit(c(1,10), "cm")) & theme(legend.position = "none")
#plot_annotation(title = "main title") 

ggsave("Figs/taxa/Biome1_SR_poisson.png", width = 13.7, height = 6)