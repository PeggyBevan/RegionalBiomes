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
source('code/BiomeModelFunctions.R')

# 3. DATA -----------------------------------------------------------------
data <- readRDS('Data/03_PREDICTSModelData_taxa.rds')
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
r0 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
r1 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|my_taxa)", REML = F)
r2 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
r3 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SSB)+(1|my_taxa)", REML = F)
r4 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SSB)", REML = F)
r5 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)", REML = F)
r6 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|my_taxa)", REML = F)
r7 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc'
                               , randomStruct = 0, REML = F)
AIC(r0$model, r1$model, r2$model, r3$model,r4$model,r5$model,r6$model)
#SS & SSB & my_taxa is best

#check best land use structure

L10 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                                fitFamily = 'poisson', fixedStruct = 'LandUse*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l11 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                                fitFamily = 'poisson', fixedStruct = 'LandUse2*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l12 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                                fitFamily = 'poisson', fixedStruct = 'LandUse3*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l13 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                                fitFamily = 'poisson', fixedStruct = 'LandUse4*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
l14 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
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
B4LUsel$deltaAIC <- B4LUsel$AIC - min(B4LUsel$AIC)
B4LUsel = arrange(B4LUsel, AIC)


# 4.3: Run Models ----------------------------------------------------------
#remove where use intensity = NA
Biome4 = Biome4[!is.na(Biome4$Use_intensity),]
Biome4 = filter(Biome4, my_taxa != 'Fungi')
Biome4 = filter(Biome4, Realm != 'Australasia')

#model selection
#using biome models function from biome1 script
B4sr<- Biomemodels1(data = Biome4, responseVar = 'Species_richness', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3', fitfamily = 'poisson')
B4A <- Biomemodels1(data = Biome4, responseVar = 'LogAbund', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3', fitfamily = 'gaussian')

B4 <- rbind(B4sr, B4A)
B4$Dataset = 'Biome4'
B4 <- rbind(B4LUsel, B4)
write.csv(B4, "output/B4ModSel.csv", row.names = F)

names(B4)[3] = 'Fixed Effects'
names(B4)[6] = 'Marginal R2'
names(B4)[10] = 'dAIC'

B4$`Marginal R2`= round(B4$`Marginal R2`, 2)
B4$AIC = round(B4$AIC, 0)
B4$dAIC = round(B4$dAIC,2)

TS6 = flextable(B4[,c(2:4,8,9,6,5,10)])
small_border = fp_border(color="black", width = 2)
tiny_border = fp_border(color="black", width = 1.5)
TS6 <- merge_v(TS6, j = ~ Response)
TS6 <- theme_vanilla(TS6)
TS6 <- fix_border_issues(TS6)
TS6 <- hline(TS6, border = small_border, i = c(5,9))
TS6
save_as_image(TS6, 'Output/TableS6_Biome4Models.png')
save_as_docx(TS6, path = 'Output/TableS6_Biome4Models.docx')

#run models to be used in figures
b0 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

b1 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)", REML = F)

b2 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LU_UI_3', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)

b3 <- StatisticalModels::GLMER(modelData = Biome4, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LU_UI_3*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)


AIC(b0$model, b1$model, b2$model, b3$model)
#adding realm improves model, land use/ use intensity is best combo

R2GLMER(b0$model)
R2GLMER(b1$model)
R2GLMER(b2$model)
R2GLMER(b3$model)

Biome4_abund = Biome4[!is.na(Biome4$LogAbund),]
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
  B4_LU$n[i] <- nrow(subset(Biome4, Realm == B4_LU$Realm[i] & LandUse == B4_LU$LU[i]))
}

level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")

figB4.a <- ggplot(B4_LU[B4_LU$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#A6D854", '#8DA0CB', "#FFD92F")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), linewidth = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), linewidth = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100), hjust = 1.1,  data = B4_LU[B4_LU$n > 25 & B4_LU$upper > 100,], colour = 'black') +
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

B4_LUUI3 <- apply(RealmB4, 1, FUN = realmPredsLU_UI, model = b3$model, data = Biome4[!is.na(Biome4$Use_intensity),])
B4_LUUI3 <- data.table::rbindlist(B4_LUUI3)
for ( i in 1:nrow(B4_LUUI3)) {
  B4_LUUI3$n[i] <- nrow(subset(Biome4[!is.na(Biome4$Use_intensity),], Realm == B4_LUUI3$Realm[i] & LU_UI_3 == B4_LUUI3$LU[i])) 
}
LandUse = c("Primary Vegetation", 'Secondary Vegetation', 'Secondary Vegetation', 'Agriculture', 'Agriculture')
UseIntensity <- c('PV', 'Minimal Use', 'Instense Use', 'Minimal Use', 'Instense Use')
B4_LUUI3 <- cbind(B4_LUUI3, LandUse, UseIntensity)

level_order <- c("Primary Vegetation", "Secondary Vegetation_Minimal use", "Secondary Vegetation_Intense use", "Agriculture_Minimal use", "Agriculture_Intense use")

figB4.b <- ggplot(B4_LUUI3[B4_LUUI3$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  scale_colour_manual(values = c("#A6D854", '#8DA0CB', "#FFD92F")) +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), linewidth = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), linewidth = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1,  data = B4_LUUI3[B4_LUUI3$n > 25 & B4_LUUI3$upper > 100,], colour = 'black') +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
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
  B4_LU_a$n[i] <- nrow(subset(Biome4_abund[!is.na(Biome4_abund$Use_intensity),], Realm == B4_LU_a$Realm[i] & LandUse == B4_LU_a$LU[i])) 
}


figB4.c <- ggplot(B4_LU_a[B4_LU_a$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#A6D854", '#8DA0CB', "#FFD92F")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), linewidth = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), linewidth = 2) +
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
  ggtitle('a')

figB4.c
#ggsave("FinalScriptsAndData/Figs/attempt2/Biome4_LU_aRealm.png")


### Fig 5.d -----------------------------------------------------------------

B4_LUUI3_a <- apply(RealmB4, 1, FUN = realmPredsLU_UI_a, model = ba3$model, data = Biome4_abund[!is.na(Biome4_abund$Use_intensity),])
B4_LUUI3_a <- data.table::rbindlist(B4_LUUI3_a)
for ( i in 1:nrow(B4_LUUI3_a)) {
  B4_LUUI3_a$n[i] <- nrow(subset(Biome4_abund[!is.na(Biome4_abund$Use_intensity),], Realm == B4_LUUI3_a$Realm[i] & LU_UI_3 == B4_LUUI3_a$LU[i])) 
}

B4_LUUI3_a <- cbind(B4_LUUI3_a, LandUse, UseIntensity)

figB4.d <- ggplot(B4_LUUI3_a[B4_LUUI3_a$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) +
  scale_colour_manual(values = c("#A6D854", '#8DA0CB', "#FFD92F")) + 
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), size = 2) +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1,  angle = 45, data = B4_LUUI3_a[B4_LUUI3_a$n > 25 & B4_LUUI3_a$upper > 100,], colour = 'black', position = position_dodge(width = 0.6)) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = 'Land Use - Use Intensity', labels = c('PV', 'SV-M', 'SV-I' ,'Agr-M', 'Agr-I')) +
  scale_y_continuous(name = '') +
  theme_bw()+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black')) +
  ggtitle('b')

figB4.d


# Fig 5 -------------------------------------------------------------------


l <- cowplot::get_legend(figB4.a + theme(legend.position = "bottom"))

#plotting species richness only, abundance for supplementary:
(cowplot::plot_grid(l) + plot_spacer())/(figB4.a + figB4.b) + 
  plot_layout(heights = unit(c(1,10), "cm")) & theme(legend.position = "none")
#plot_annotation(title = "main title") 
ggsave("Figs/Fig5_SR.png", width = 13.7, height = 12)
#make sure graphics box is tall & wide enough to stop text overlapping.
#Saving 13.7 x 12 in image

#abundance plots for supplementary:
(cowplot::plot_grid(l) + plot_spacer())/(figB4.c + figB4.d) + 
  plot_layout(heights = unit(c(1,10), "cm")) & theme(legend.position = "none")
#plot_annotation(title = "main title") 
ggsave("Figs/FigS2_A.png", width = 13.7, height = 12)

#All 4 plots:
# cowplot::plot_grid(l)/(figB1.a + figB1.b + figB1.c + figB1.d) + 
#   plot_layout(heights = unit(c(1,20), "cm")) & theme(legend.position = "none")
# #plot_annotation(title = "main title") 
# ggsave("Figs/Fig4.png", width = 13.7, height = 12)
# #make sure graphics box is tall & wide enough to stop text overlapping.
# #Saving 13.7 x 12 in imagewrite.csv(Biome4, "Figs/Maps/Biome4.csv", row.names = F)

write.csv(Biome4, "Figs/Maps/Biome4.csv", row.names = F)

# Taxon case study --------------------------------------------------------

#I think the ranges of individual species might be varying a lot, so some kind of scaling is needed

unique(Biome4$my_taxa)
table(Biome4$LandUse, Biome4$my_taxa)

Biome4 %>%
  group_by(my_taxa) %>%
  summarise(min = min(Species_richness), max = max(Species_richness))

#just vertebrates
Biome4_v = filter(Biome4, my_taxa == 'Vertebrate')

b4v <- StatisticalModels::GLMER(modelData = Biome4_v, 
                                responseVar = "Species_richness",
                                fitFamily = 'poisson', 
                                fixedStruct = 'LandUse*Realm',
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
#list realms in this biome
RealmB4v <- data.frame(Realm = unique(Biome4_v$Realm))

#check its the right realms
test <- sapply(RealmB4v, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})


B4_LUv<- apply(RealmB4v, 1, FUN = realmPredsLU, model = b4v$model, data = Biome4_v)
B4_LUv <- data.table::rbindlist(B4_LUv)
for ( i in 1:nrow(B4_LUv)) {
  B4_LUv$n[i] <- nrow(subset(Biome4_v, Realm == B4_LUv$Realm[i] & LandUse == B4_LUv$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")


figB4v.a <- ggplot(B4_LUv[B4_LUv$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  scale_colour_manual(values = c("#A6D854", '#8DA0CB', "#FFD92F")) + 
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.7), size = 0.7) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), size = 1.5) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr'), drop = F) +
  scale_y_continuous(name = 'Change in Species Richness (%)') +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1,  angle = 45, data = B4_LUv[B4_LUv$n > 25 & B4_LUv$upper > 100,], colour = 'black', position = position_dodge(width = 0.6)) +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    #legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black')) +
  ggtitle("a")

figB4v.a

#just invertebrates
Biome4_i = filter(Biome4, my_taxa == 'Invertebrate')

B4i <- StatisticalModels::GLMER(modelData = Biome4_i, 
                                responseVar = "Species_richness",
                                fitFamily = 'poisson', 
                                fixedStruct = 'LandUse*Realm',
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
#list realms in this biome
RealmB4i <- data.frame(Realm = unique(Biome4_i$Realm))

#check its the right realms
test <- sapply(RealmB4i, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})

B4_LUi<- apply(RealmB4i, 1, FUN = realmPredsLU, model = B4i$model, data = Biome4_i)
B4_LUi <- data.table::rbindlist(B4_LUi)
for ( i in 1:nrow(B4_LUi)) {
  B4_LUi$n[i] <- nrow(subset(Biome4_i, Realm == B4_LUi$Realm[i] & LandUse == B4_LUi$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")
B4_LUi

#the sample size is slightly higher here because neotropic pasture has infinite values, so presenting the resultso f this are not helpful. 
figB4i.a <- ggplot(B4_LUi[B4_LUi$n > 36,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  scale_colour_manual(values = c("#A6D854", '#8DA0CB', "#FFD92F")) + 
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.7), linewidth = 0.7) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), linewidth = 1.5) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr'), drop = F) +
  scale_y_continuous(name = '') +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1,  angle = 45, data = B4_LUi[B4_LUi$n > 36 & B4_LUi$upper > 100,], colour = 'black', position = position_dodge(width = 0.6)) +
  #geom_text(aes(y = lower ), position="identity", color = 'black', size = 3, angle = 45) +       
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    #legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black'))+
  ggtitle("b")

figB4i.a
figB4v.a

#just plants
Biome4_p = filter(Biome4, my_taxa == 'Plant')

B4p <- StatisticalModels::GLMER(modelData = Biome4_p, 
                                responseVar = "Species_richness",
                                fitFamily = 'poisson', 
                                fixedStruct = 'LandUse*Realm',
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
#n = 5536
#list realms in this biome
RealmB4p <- data.frame(Realm = unique(Biome4_p$Realm))

#check its the right realms
test <- sapply(RealmB4p, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})


B4_LUp<- apply(RealmB4p, 1, FUN = realmPredsLU, model = B4p$model, data = Biome4_p)
B4_LUp <- data.table::rbindlist(B4_LUp)
for ( i in 1:nrow(B4_LUp)) {
  B4_LUp$n[i] <- nrow(subset(Biome4_p, Realm == B4_LUp$Realm[i] & LandUse == B4_LUp$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")


figB4p.a <- ggplot(B4_LUp[B4_LUp$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  scale_colour_manual(values = c('#8DA0CB', "#FFD92F")) + geom_errorbar(width = 0.2, position = position_dodge(width = 0.7), linewidth = 0.7) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), linewidth = 0.7) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), linewidth = 1.5) +
  coord_cartesian(ylim = c(-100,100)) +
  scale_x_discrete(name = '', labels = c('PV', 'SV', 'PF', 'Pa' ,'Cr'), drop = F) +
  scale_y_continuous(name = '') +
  geom_text(aes(x = LU, group = Realm, label = signif(upper,3), y =100, ), hjust = 1,  angle = 45, data = B4_LUp[B4_LUp$n > 25 & B4_LUp$upper > 100,], colour = 'black', position = position_dodge(width = 0.6)) +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    #legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 20, colour = 'black'))+
  ggtitle("c")

figB4p.a
figB4i.a
figB4v.a

l <- cowplot::get_legend(figB4i.a + theme(legend.position = "bottom"))

(cowplot::plot_grid(l)+plot_spacer())/(figB4v.a + figB4i.a + figB4p.a) + 
  plot_layout(heights = unit(c(1,10), "cm")) & theme(legend.position = "none")
#plot_annotation(title = "main title") 

ggsave("Figs/taxa/Biome4_SR_poisson.png", width = 13.7, height = 6)




# Plotting whole biome ----------------------------------------------------
#Not in manuscript, just for visualisation & data checking:
#whole biome, just land use 

nd <- data.frame(LandUse = factor(levels(b0$data$LandUse),
                                  levels = levels(b0$data$LandUse)),
                 Species_richness = 0)

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
                 Species_richness = 0)

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


