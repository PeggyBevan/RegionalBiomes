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
source('code/BiomeModelFunctions.R')

# 2. DATA -----------------------------------------------------------------
data <- readRDS('Data/03_PREDICTSModelData.rds')
dim(data)

# 3. SCRIPT ---------------------------------------------------------------
## 3.1 Prep and subset data ------------------------------------------------
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

Biome7 = Biome7[!is.na(Biome7$Use_intensity),]
## 4.2 Selecting model structure -------------------------------------------
#check random effects & land use structure

#test random effects
r0 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
r1 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|my_taxa)", REML = F)
r2 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)+(1|SSB)", REML = F)
r3 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SSB)+(1|my_taxa)", REML = F)
r4 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SSB)", REML = F)
r5 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|SS)", REML = F)
r6 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc', 
                               randomStruct = "(1|my_taxa)", REML = F)
r7 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*RB_tnc'
                               , randomStruct = 0, REML = F)
AIC(r0$model, r1$model, r2$model, r3$model,r4$model,r5$model,r6$model)
#SS & SSB is best

#check best land use structure

L70 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l71 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse2*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l72 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse3*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)

l73 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse4*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)
l74 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "Species_richness",
                               fitFamily = 'poisson', fixedStruct = 'LandUse5*Realm', 
                               randomStruct = "(1|SS)+(1|SSB)+(1|my_taxa)", REML = F)



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
B7LUsel$deltaAIC <- B7LUsel$AIC - min(B7LUsel$AIC)
B7LUsel = arrange(B7LUsel, AIC)


## 4.3 Run Models ----------------------------------------------------------

#model selection
#using biome models function from biome1 script
B7sr<- Biomemodels1(data = Biome7, responseVar = 'Species_richness', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3', fitfamily = 'poisson')
B7A <- Biomemodels1(data = Biome7, responseVar = 'LogAbund', LandUseVar = "LandUse", UseIntensityVar = 'LU_UI_3', fitfamily = 'gaussian')

B7 <- rbind(B7sr, B7A)
B7$Dataset = 'Biome7'
B7 <- rbind(B7LUsel, B7)
write.csv(B7, "output/B7ModSel.csv", row.names = F)

names(B7)[3] = 'Fixed Effects'
names(B7)[6] = 'Marginal R2'
names(B7)[10] = 'dAIC'
B7$`Marginal R2`= round(B7$`Marginal R2`, 2)
B7$AIC = round(B7$AIC, 0)
B7$dAIC = round(B7$dAIC,2)


TS7 = flextable(B7[,c(2:4,8,9,6,5,10)])
small_border = fp_border(color="black", width = 2)
tiny_border = fp_border(color="black", width = 1.5)
TS7 <- merge_v(TS7, j = ~ Response)
TS7 <- theme_vanilla(TS7)
TS7 <- fix_border_issues(TS7)
TS7 <- hline(TS7, border = small_border, i = c(5,9))
TS7
save_as_image(TS7, 'Output/TableS7_Biome7Models.png')
save_as_docx(TS7, path = 'Output/TableS7_Biome7Models.docx')

b70 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "Species_richness",
                                fitFamily = 'poisson', fixedStruct = 'LandUse', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

b71 <- StatisticalModels::GLMER(modelData = Biome7[!is.na(Biome7$Use_intensity),], responseVar = "Species_richness",
                                fitFamily = 'poisson', fixedStruct = 'LandUse*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

b72 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "Species_richness",
                                fitFamily = 'poisson', fixedStruct = 'LU_UI_3', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)

b73 <- StatisticalModels::GLMER(modelData = Biome7, responseVar = "Species_richness",
                                fitFamily = 'poisson', fixedStruct = 'LU_UI_3*Realm', 
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)


AIC(b70$model, b71$model, b72$model, b73$model)

#Abundance
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

## 4.3 Plot Models ----------------------------------------------------------

RealmB7 <- data.frame(Realm = unique(Biome7$Realm))

test <- sapply(RealmB7, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})

# Fig 6.a -----------------------------------------------------------------

B7_LU<- apply(RealmB7, 1, FUN = realmPredsLU, model = b71$model, data = Biome7)
B7_LU <- data.table::rbindlist(B7_LU)
for ( i in 1:nrow(B7_LU)) {
  B7_LU$n[i] <- nrow(subset(Biome7, Realm == B7_LU$Realm[i] & LandUse == B7_LU$LU[i]))
}

level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")

figB7.a <- ggplot(B7_LU[B7_LU$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), linewidth = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), linewidth = 2) +
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


B7_LUUI3 <- apply(RealmB7, 1, FUN = realmPredsLU_UI, model = b73$model, data = Biome7)
B7_LUUI3 <- data.table::rbindlist(B7_LUUI3)
for ( i in 1:nrow(B7_LUUI3)) {
  B7_LUUI3$n[i] <- nrow(subset(Biome7, Realm == B7_LUUI3$Realm[i] & LU_UI_3 == B7_LUUI3$LU[i])) 
}
B7_LUUI3 <- cbind(B7_LUUI3, LandUse, UseIntensity)


figB7.b <- ggplot(B7_LUUI3[B7_LUUI3$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), linewidth = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), linewidth = 2) +
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
  B7_LU_a$n[i] <- nrow(subset(Biome7_abund[!is.na(Biome7_abund$Use_intensity),], Realm == B7_LU_a$Realm[i] & LandUse == B7_LU_a$LU[i])) 
}

figB7.c <- ggplot(B7_LU_a[B7_LU_a$n > 25,], aes(x = factor(LU, levels = level_order_LU), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) + 
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), linewidth = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), linewidth = 2) +
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
  B7_LUUI3_a$n[i] <- nrow(subset(Biome7_abund[!is.na(Biome7_abund$Use_intensity),], Realm == B7_LUUI3_a$Realm[i] & LU_UI_3 == B7_LUUI3_a$LU[i])) 
}

B7_LUUI3_a <- cbind(B7_LUUI3_a, LandUse, UseIntensity)

figB7.d <- ggplot(B7_LUUI3_a[B7_LUUI3_a$n > 25,], aes(x = factor(LU, levels = LandUseFact), y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.6), size = 5) +
  scale_colour_manual(values = c("#66C2A5", "#E78AC3", "#8DA0CB")) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.6), linewidth = 1) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.6), linewidth = 2) +
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
ggsave("Figs/Fig6.png", width = 13.7, height = 12)
#make sure graphics box is tall & wide enough to stop text overlapping.
#Saving 13.7 x 12 in image


write.csv(Biome7, "Figs/Maps/Biome7.csv")

# Taxon case study --------------------------------------------------------
#I think the ranges of individual species might be varying a lot, so some kind of scaling is needed

unique(Biome7$my_taxa)
table(Biome7$LandUse, Biome7$my_taxa)

Biome7 %>%
  group_by(my_taxa) %>%
  summarise(min = min(Species_richness), max = max(Species_richness))

#just vertebrates
Biome7_v = filter(Biome7, my_taxa == 'Vertebrate')

b7v <- StatisticalModels::GLMER(modelData = Biome7_v, 
                                responseVar = "Species_richness",
                                fitFamily = 'poisson', 
                                fixedStruct = 'LandUse*Realm',
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
#list realms in this biome
Realmb7v <- data.frame(Realm = unique(Biome7_v$Realm))

#check its the right realms
test <- sapply(Realmb7v, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})


b7_LUv<- apply(Realmb7v, 1, FUN = realmPredsLU, model = b7v$model, data = Biome7_v)
b7_LUv <- data.table::rbindlist(b7_LUv)
for ( i in 1:nrow(b7_LUv)) {
  b7_LUv$n[i] <- nrow(subset(Biome7_v, Realm == b7_LUv$Realm[i] & LandUse == b7_LUv$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")


figb7v.a <- ggplot(b7_LUv[b7_LUv$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  scale_colour_manual(values = c("#A6D854", '#8DA0CB', "#FFD92F")) + 
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.7), size = 0.5) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), size = 1) +
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
  ggtitle("Vertebrates")

figb7v.a

#just invertebrates
Biome7_i = filter(Biome7, my_taxa == 'Invertebrate')

b7i <- StatisticalModels::GLMER(modelData = Biome7_i, 
                                responseVar = "Species_richness",
                                fitFamily = 'poisson', 
                                fixedStruct = 'LandUse*Realm',
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
#list realms in this biome
Realmb7i <- data.frame(Realm = unique(Biome7_i$Realm))

#check its the right realms
test <- sapply(Realmb7i, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})

b7_LUi<- apply(Realmb7i, 1, FUN = realmPredsLU, model = b7i$model, data = Biome7_i)
b7_LUi <- data.table::rbindlist(b7_LUi)
for ( i in 1:nrow(b7_LUi)) {
  b7_LUi$n[i] <- nrow(subset(Biome7_i, Realm == b7_LUi$Realm[i] & LandUse == b7_LUi$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")


figb7i.a <- ggplot(b7_LUi[b7_LUi$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  scale_colour_manual(values = c("#E78AC3", "#A6D854", '#8DA0CB', "#FFD92F")) + 
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.7), size = 0.5) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), size = 1) +
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
  ggtitle("Invertebrates")

figb7i.a
figb7v.a

#just plants
Biome7_p = filter(Biome7, my_taxa == 'Plant')

b7p <- StatisticalModels::GLMER(modelData = Biome7_p, 
                                responseVar = "Species_richness",
                                fitFamily = 'poisson', 
                                fixedStruct = 'LandUse*Realm',
                                randomStruct = "(1|SS)+(1|SSB)", REML = F)
#n = 5536
#list realms in this biome
Realmb7p <- data.frame(Realm = unique(Biome7_p$Realm))

#check its the right realms
test <- sapply(Realmb7p, FUN = function(realm) {
  cat(paste0(realm,'\n'))
})


b7_LUp<- apply(Realmb7p, 1, FUN = realmPredsLU, model = b7p$model, data = Biome7_p)
b7_LUp <- data.table::rbindlist(b7_LUp)
for ( i in 1:nrow(b7_LUp)) {
  b7_LUp$n[i] <- nrow(subset(Biome7_p, Realm == b7_LUp$Realm[i] & LandUse == b7_LUp$LU[i]))
}
level_order_LU <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Pasture","Cropland")


figb7p.a <- ggplot(b7_LUp[b7_LUp$n > 25,], aes(x = factor(LU, levels = level_order_LU),y = y, ymax = upper, ymin = lower, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  scale_colour_manual(values = c('#8DA0CB', "#FFD92F")) + geom_errorbar(width = 0.2, position = position_dodge(width = 0.7), size = 0.5) +
  geom_errorbar(width = 0, aes(ymin = lower75, ymax = upper75), position = position_dodge(width = 0.7), size = 1) +
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
  ggtitle("Plants")

figb7p.a
figb7i.a
figb7v.a

l <- cowplot::get_legend(figb7i.a + theme(legend.position = "bottom"))

cowplot::plot_grid(l)/(figb7v.a + figb7i.a + figb7p.a) + 
  plot_layout(heights = unit(c(1,10), "cm")) & theme(legend.position = "none")
#plot_annotation(title = "main title") 

ggsave("Figs/taxa/Biome7_SR_poisson.png", width = 13.7, height = 6)



