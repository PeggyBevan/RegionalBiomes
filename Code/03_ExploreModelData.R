# 03_ExploreData
#Order of code:
# 1. Packages
# 2. Load Data
# 3. Script
  #Regional biomes
#  - creating ecoregs - list of all existing RBs & how many studies are in BD for each one
#  - save 01_RBSamplingEffort.png - the table that shows representation across regional biomes in the predicts data base - before anything has been removed. Table S1
  
  #Land Use
# - create a table to show land-use sampling effort, including 'human-dominated' - the sum of pasture, cropland and plantation forest.
# - saved as 05_RBLandUseSummary.csv
  #Taxon Representation
# - exploring sampling effort by taxon, in each biome, realm and regional biome
# - create table '03_TaxaSamplingEffort.png'
  #Summary table of RBs & Save
# - create summary table of the final regional biomes left - because i have not removed much this is basically a summary of RBs. '04_FinalList.png'
  #Regional Biome csv
# - save a csv with columns biome, biome number, realm, RB code, n_studies, n_sources, N_sites, n_taxa, n_ all land use types, total ecoregions, proportion of ecoregions covered. and number of observations in each taxa. Saved to data/04_RBsummary.csv
# - used later on to determine best sample size threshold for regional biome inclusion criteria. 


# Packages ----------------------------------------------------------------
library(rgdal)  # readOGR() spTransform()
library(dplyr) # n_distinct(), count()
library(flextable) #flextable()
library(tidyr) #spread()
library(officer) #prop_section, for saving tables as .docx

# Data --------------------------------------------------------------------

BD <- readRDS('Data/03_PREDICTSModelData_taxa.rds')

# The TNC egoregion map has the same ecoregions as those in the PREDICTS databse. I got it from here: http://maps.tnc.org/gis_data.html
tnc_ecoregions <- readOGR(dsn = 'Data/terr-ecoregions-TNC', layer = 'tnc_terr_ecoregions')
# 814 features, 16 fields


# Script ------------------------------------------------------------------

n_rb <- n_distinct(BD$RB_tnc)
 #There are 46 regional biomes represented by the Predicts database.

#Spatial Representation of regional biomes within the PREDICTS database
#I want to look at: 

# * number of studies in each RB  (n_studies)
# * number of ecoregions covered within each RB (n_ecoregions) as a proportion of the total available (t_ecoregions); (p_ecoregions) 
# * proportion of observations at each land cover type  (n_pv; n_sv; n_pf; n_crl; n_pas; n_urb)
# * number of taxa covered  (n_taxa)
# 
# In this shapefile, regional biomes are given a code with their geographic realm and biome number. For example, Tropical and Subtropical Moist Broadleaf forest (biome 01) in the indo-malay realm has the code **IM1**.

#list all RBs
rb_tncCodes<-unique(tnc_ecoregions@data$RealmMHT)

# make dataframe with realm code and biome
ecoregs <- data.frame('RB_tnc' = unique(tnc_ecoregions@data$RealmMHT)) %>%
  mutate(Realm_code = substr(RB_tnc, 1,2),
         Biome_num = ifelse(nchar(as.character(RB_tnc))<=3, 
                            substr(RB_tnc, 3,3),
                            substr(RB_tnc,3,4)
         ),
         Realm = recode_factor(Realm_code,
                               'IM' = 'Indo-Malay',
                               'PA' = 'Palearctic',
                               'AA' = 'Australasia',
                               'NA' = 'Nearctic',
                               'NT' = 'Neotropic',
                               'AT' = 'Afrotropic',
                               'OC' = 'Oceania',
                               'AN' = 'Antarctica' 
         ),
         Biome = recode_factor(Biome_num,
                               "1" = 'Tropical & Subtropical Moist Broadleaf Forests',
                               "4" = 'Temperate Broadleaf & Mixed Forests',
                               "12" = 'Mediterranean Forests, Woodlands & Scrub',
                               "7" = 'Tropical & Subtropical Grasslands, Savannas & Shrublands',
                               "5" = 'Temperate Conifer Forests',
                               "8" = 'Temperate Grasslands, Savannas & Shrublands',
                               "13" = 'Deserts & Xeric Shrublands',
                               "14" = 'Mangroves',
                               "2" = 'Tropical & Subtropical Dry Broadleaf Forests',
                               "10" = 'Montane Grasslands & Shrublands',
                               "3" = 'Tropical & Subtropical Coniferous Forests',
                               "6" = 'Boreal Forests/Taiga',
                               "9" = 'Flooded Grasslands & Savannas',
                               "11" = 'Tundra',
                               "98" = "Inland Water", 
                               "99" = "Rock and Ice"
         )
  )


ecoregs$Biome_num<-as.numeric(ecoregs$Biome_num)
ecoregs<-arrange(ecoregs, Biome_num, Realm)
#remove rock & ice and inland water
ecoregs <- filter(ecoregs, Biome_num<=14)

# Regional Biomes ------------------------------------------------------
# calculate number of studies etc. for each RB
for (i in (1:nrow(ecoregs))) {
  rb = ecoregs$RB_tnc[i]
  tnc_code = paste(ecoregs$RB_tnc[i])
  ecoregs$n_studies[i] = n_distinct(BD$Study_name[BD$RB_tnc==rb])
  ecoregs$n_sources[i] = n_distinct(BD$Source_ID[BD$RB_tnc==rb])
  ecoregs$n_sites[i] = n_distinct(BD$SSBS[BD$RB_tnc==rb]) 
  ecoregs$n_ecoregions[i] = n_distinct(BD$Ecoregion[BD$RB_tnc==rb])
  ecoregs$n_taxa[i] = n_distinct(BD$my_taxa[BD$RB_tnc==rb])
  ecoregs$n_country[i] =n_distinct(BD$Country[BD$RB_tnc==rb])
  ecoregs$t_ecoregions[i] = n_distinct(tnc_ecoregions@data$ECO_NAME[tnc_ecoregions@data$RealmMHT==tnc_code])
  ecoregs$p_ecoregions[i] = round(ecoregs$n_ecoregions[i]/ecoregs$t_ecoregions[i], 2)
}

write.csv(ecoregs, 'Data/04_RBsummary.csv', row.names = F)

#That ends up looking like this:

#taxa is not included in this as it is just geographic representation
###TABLE 1
#if starting code from here:
ecoregs <- read.csv("Data/04_RBsummary.csv", stringsAsFactors = T)
ecoregs$Biome_num<-as.numeric(ecoregs$Biome_num)
ecoregs<-arrange(ecoregs, Biome_num, Realm)
f1 <- flextable(ecoregs[,c(5,4,1,6:13)]) 

typology <- data.frame(
  col_keys = f1$col_keys,
  type = c("Regional Biome", "Regional Biome",
           "Regional Biome",
           'Sampling Effort', 'Sampling Effort', 
           'Sampling Effort', 'Sampling Effort', 'Sampling Effort', 'Sampling Effort',
           'Proportion of Ecoregions sampled', 'Proportion of Ecoregions sampled'),
  what = c("Biome", "Realm", "Code", "Studies", "Sources", 
           'Sites', 'Ecoregions', 'Taxa', "Countries", 'Total ecoregions',
           'Proportion of Ecoregions sampled'),
  stringsAsFactors = FALSE )
f1
formatflext <- function(ftable) {
  ftable <- set_header_df(ftable, mapping = typology, key = 'col_keys')
  ftable <- merge_h(ftable, part = "header")
  ftable <- merge_v(ftable, part = "header")
  ftable <- align(ftable, part = 'header', align = 'left')
  ftable <- theme_vanilla(ftable)
}

f1 <- formatflext(f1)
f1 = align(f1, part = 'header', align = 'left')
f1 <- merge_v(f1, j = ~ Biome)
f1 <- fix_border_issues(f1)
#conditional formatting
f1 <- bg(f1, bg = '#cccccc', i = ~ n_studies < 1, j = 2:11)
f1
#save_as_docx(f1, path = 'testtable.docx')
#Create fig directory if it doesn't exist
if dir.exists('Figs') == FALSE {
  dir.create('Figs')
}
save_as_image(f1, path = 'Figs/01_RBSamplingEffort.png')
save_as_docx(f1, pr_section = prop_section(page_size(orient = 'landscape')), path = 'Figs/01_RBSamplingEffort.docx')

#split this into two tables, so it can go over two pages
part1 = ecoregs[1:45,]
part2 = ecoregs[46:64,]

TS1a = flextable(part1[,c(5,4,1,6:13)]) 
TS1b = flextable(part2[,c(5,4,1,6:13)]) 

TS1a <- formatflext(TS1a)
TS1a = align(TS1a, part = 'header', align = 'left')
TS1a <- merge_v(TS1a, j = ~ Biome)
TS1a <- fix_border_issues(TS1a)
#conditional formatting
TS1a <- bg(TS1a, bg = '#cccccc', i = ~ n_studies < 1, j = 2:11)
TS1a
#add cont. footer
TS1a = add_footer_row(TS1a, values = 'cont.', colwidths = 11)
TS1a = align(TS1a, part = 'footer', align = 'right')
#save_as_docx(f1, path = 'testtable.docx')

if dir.exists('Output') == FALSE {
  dir.create('Output')
}
save_as_image(TS1a, path = 'Output/TS1a_RegionalBiomeList.png')
save_as_docx(TS1a, path = 'Output/TS1a_RegionalBiomeList.docx')

TS1b <- formatflext(TS1b)
TS1b = align(TS1b, part = 'header', align = 'left')
TS1b <- merge_v(TS1b, j = ~ Biome)
TS1b <- fix_border_issues(TS1b)
#conditional formatting
TS1b <- bg(TS1b, bg = '#cccccc', i = ~ n_studies < 1, j = 2:11)
TS1b
#save_as_docx(f1, path = 'testtable.docx')
save_as_image(TS1b, path = 'Output/TS1b_RegionalBiomeList.png')
save_as_docx(TS1b, path = 'Output/TS1b_RegionalBiomeList.docx')



# Land Use ----------------------------------------------------------------

LU <- BD %>% count(RB_tnc, LandUse) %>% 
  spread(LandUse, n)
ecoregs<- merge(ecoregs, LU, by.x = 'RB_tnc', all.x = TRUE)

#create agriculture count
ecoregs$Agriculture <- rowSums(ecoregs[,16:18], na.rm = T)

# Taxon Sampling Effort ---------------------------------------------------

# representation by taxa
levels(BD$my_taxa)
n_distinct(BD$my_taxa)

table(BD$my_taxa)

##How are taxa spread across regional biomes?
# Across biomes
d <- as.data.frame(table(BD$Biome, BD$my_taxa))
d <- spread(d, Var2,Freq)
names(d)[1]<-'Biome'
flextable(d)
# fungi and herptiles mainly underrepresented. Mammals are less studied in grassland biomes. Birds and plants are very well sampled.

# Across realms
e <- as.data.frame(table(BD$Realm, BD$my_taxa))
e <- spread(e, Var2, Freq)
names(e)[1]<-'Realm'
flextable(e)

#fungi are not well studied anywhere tropical. 
# no mammal studies in nearctic realm. 
#everywhere else pretty well spread out.
table(BD$)

Taxa <- BD %>% count(RB_tnc, my_taxa) %>% spread(my_taxa, n)
ecoregs<- merge(ecoregs, Taxa, by.x = 'RB_tnc', all.x = TRUE)

f4 <- ecoregs[,c(5,4,1,3,20:23)]
f4 <- arrange(f4, Biome_num)
f4[is.na(f4)] = 0 
f4 <- flextable(f4[,c(1,2,5:8)])
f4 <- merge_v(f4, j = ~ Biome)

typology <- data.frame(
  col_keys = f4$col_keys,
  type = c("Regional Biome", "Regional Biome",
           'Number of samples of each taxa', 
           'Number of samples of each taxa', 
           'Number of samples of each taxa', 
           'Number of samples of each taxa'
           ),
  what = c("Biome", "Realm", 
           'Plant', "Fungi",  
           'Invertebrate', 
           'Vertebrate'),
  stringsAsFactors = FALSE )

f4 <- set_header_df(f4, mapping = typology, key = 'col_keys')
f4 <- merge_h(f4, part = "header")
f4 <- merge_v(f4, part = "header")
f4 <- align(f4, i = 1:2, part = 'header', align = 'left')
f4 <- theme_vanilla(f4)
f4 <- fix_border_issues(f4)
f4
save_as_image(f4, 'Figs/03_TaxaSamplingEffort.png')
save_as_docx(f4, path = 'Figs/03_TaxaSamplingEffort.docx')

#Some regional biomes are only represented by one or two taxa groups, meaning they probably aren't representative of the RB as a whole. 
# --we have to assume that the representation of taxa are proportional to the number of taxa in that RB.

#supplementary fig - taxon per regional biome, only studied biomes
ecoregs_dd <- ecoregs %>%
  subset(Biome != "Tundra") %>%
  subset(Biome != "Flooded Grasslands & Savannas") %>%
  subset(Biome != "Mangroves") %>%
  subset(n_sites > 1)

ecoregs_dd = subset(ecoregs_dd, Biome_num == 1 | Biome_num == 4 | Biome_num == 7)

ft <- ecoregs_dd[,c(5,4,3,14:17)]
ft <- arrange(ft, Biome_num)
ft[is.na(ft)] = 0 
ft <- flextable(ft[,c(1,2,4:7)])
ft <- merge_v(ft, j = ~ Biome)

typology <- data.frame(
  col_keys = ft$col_keys,
  type = c("Regional Biome", "Regional Biome",
           'Number of sites', 
           'Number of sites', 
           'Number of sites', 
           'Number of sites'
  ),
  what = c("Biome", "Realm", 
           'Fungi',  
           'Invertebrate', 'Plant', 
           'Vertebrate'),
  stringsAsFactors = FALSE )

ft <- set_header_df(ft, mapping = typology, key = 'col_keys')
ft <- merge_h(ft, part = "header")
ft <- merge_v(ft, part = "header")
ft <- align(ft, i = 1:2, part = 'header', align = 'left')
ft <- theme_vanilla(ft)
ft <- fix_border_issues(ft)
ft
save_as_image(ft, 'Figs/03_TaxaSamplingEffort.png')
save_as_docx(ft, path = 'Figs/03_TaxaSamplingEffort.docx')

# 7. Summary table of RBs & Save -------------------------------------------------


BD$Biome_num <- as.numeric(BD$Biome_num)
#samples per regional biome
g <- BD %>%
  group_by(Biome_num, Biome, Realm, RB_tnc) %>%
  count()

#number of studies
n_distinct(BD$Study_name)
#number of independent sources
n_distinct(BD$Source_ID)

names(g)[5] <- 'N_sites'

f5 <- flextable(g[2:5])
f5 <- theme_vanilla(f5)
f5 <- merge_v(f5, j = ~ Biome)
f5 <- fix_border_issues(f5)
save_as_image(f5, 'Figs/04_FinalList.png')

# Save summary csv --------------------------------------------------------
##create threshold variables 

ecoregs <- ecoregs %>%
  mutate(sites.th1 = if_else(n_sites > 1, 1, 0),
         sites.th25 = if_else(n_sites > 25, 1, 0),
         sites.th50 = if_else(n_sites > 50, 1, 0),
         PV.th1 = if_else(`Primary Vegetation`> 1, 1, 0, missing = 0),
         PV.th5 = if_else(`Primary Vegetation`> 5, 1, 0, missing = 0),
         PV.th25 =  if_else(`Primary Vegetation`> 25, 1, 0, missing = 0),
         PV.th50 = if_else(`Primary Vegetation`> 50, 1, 0, missing = 0),
         Ag.th1 = if_else(Agriculture > 1, 1, 0, missing = 0),
         Ag.th5 = if_else(Agriculture > 5, 1, 0, missing = 0),
         Ag.th25 = if_else(Agriculture > 25, 1, 0, missing = 0),
         Ag.th50 = if_else(Agriculture > 50, 1, 0, missing = 0)
         )
#save ecoregs
write.csv(ecoregs, "Data/04_RBsummary_taxa.csv", row.names = F)

#for supplementary material

sp.ecos <- ecoregs[,c(1:6,8,9,11,12,13)]
write.csv(sp.ecos, "Figs/RBcoverage_taxa.csv", row.names = F)


# Table S2 - LandUse Explanation Table -----------------------------------------

TS2df = data.frame('Predominant habitat' = 
                     c('Primary vegetation',
                       'Mature secondary vegetation',
                       'Intermediate secondary vegetation',
                       'Young secondary vegetation',
                       'Plantation forest',
                       'Pasture',
                       'Cropland',
                       'Urban'),
                   'Description' = c('Native vegetation that is not known or inferred to have ever been completely destroyed by human actions or by extreme natural events that do not normally play a role in ecosystem dynamics',
                                     "Regeneration following complete removal of primary vegetation; architectural structure approaching that of primary vegetation, corresponding to a completed succession",
                                   "Regeneration following complete removal of primary vegetation; mixed architecture showing a mid-successional stage",
                                   "Regeneration following complete removal of primary vegetation; simple architecture representing an early successional stage",
                                   "Previously cleared areas that people have planted with crop trees or crop shrubs for commercial or subsistence harvesting of wood and/or fruit",
                                   "Land where livestock is known to be grazed regularly or permanently. The plant species may be predominantly native (as in rangelands) or strongly associated with humans (as in European-style pastures)",
                                   "Land that people have planted with herbaceous crops, even if these crops will be fed to livestock once harvested",
                                   "Areas with human habitation and/or buildings, where the primary vegetation has been removed, and where such vegetation as is present is predominantly managed for civic or personal amenity"),
                 'LandUse' = c('Primary Vegetation','Secondary Vegetation','Secondary Vegetation','Secondary Vegetation', 'Plantation forest', 'Pasture', 'Cropland', 'Not included'),
                 'LandUse2' = c('Natural Vegetation','Natural Vegetation','Natural Vegetation','Natural Vegetation', 'Agriculture', 'Agriculture', 'Agriculture', 'Not included'),
                 "LandUse3" = c('Primary Vegetation','Secondary Vegetation','Secondary Vegetation','Secondary Vegetation', 'Agriculture', 'Agriculture', 'Agriculture', 'Not included'),
                 "LandUse4" = c('Primary Vegetation', 'Secondary Vegetation', 'Secondary Vegetation','Secondary Vegetation', 'Harvested', 'Pasture', 'Harvested', 'Not included'),
                 "LandUse5" = c('Primary Vegetation','Secondary Vegetation','Secondary Vegetation','Secondary Vegetation', 'Plantation forest', 'Agriculture', 'Agriculture', 'Not included')
)

TS2 = flextable(TS2df)
#merge some columns
typology <- data.frame(
  col_keys = TS2$col_keys,
  higher = c("", "",
           "New classification",
           'New classification', 'New classification', 
           'New classification', 'New classification'),
  what = TS2$col_keys,
  stringsAsFactors = FALSE)
TS2
formatflext <- function(ftable) {
  ftable <- set_header_df(ftable, mapping = typology, key = 'col_keys')
  ftable <- merge_h(ftable, part = "header")
  ftable <- merge_v(ftable, part = "header")
  ftable <- align(ftable, part = 'header', align = 'left')
  ftable <- theme_vanilla(ftable)
}

TS2 <- formatflext(TS2)
TS2 = set_header_labels(TS2, "Predominant.habitat" = 'Predominant habitat')
TS2 = set_table_properties(TS2, layout = "autofit")

TS2 <- merge_v(TS2, j = ~ LandUse + LandUse2 + LandUse3 + LandUse4 + LandUse5)
TS2 <- fix_border_issues(TS2)
save_as_image(TS2, 'output/TS2_LandUseDescriptors.png')
save_as_docx(TS2, path = 'output/TS2_LandUseDescriptors.docx')
