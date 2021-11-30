# 03_ExploreData

#4. Regional biomes
#  - creating ecoregs - list of all existing RBs & how many studies are in BD for each one
#  - save 01_RBSamplingEffort.png - the table that shows representation across regional biomes in the predicts data base - before anything has been removed.
# 5. Land Use
# - create a table to show land use sampling effort, including 'human-dominated' - the sum of pasture, cropland and plantation forest.
# - this has been saved as 05_RBLandUseSummary.csv
# 6. Taxon Representation
# - exploring sampling effort by taxon, in each biome, realm and regional biome
# - create table '03_TaxaSamplingEffort.png'
# 7. Summary table of RBs & Save
# - create summary table of the final regional biomes left - because i have not removed much this is basically a summary of RBs. '04_FinalList.png'

# i want a csv with columns biome, biome number, realm, RB code, n_studies, n_sources, N_sites, n_taxa, n_ all land use types, total ecoregions, proportion of ecoregions covered.


# Packages ----------------------------------------------------------------

library(rgdal)  # readOGR() spTransform()
library(dplyr) # n_distinct(), count()
library(flextable) #flextable()
library(tidyr) #spread()


# Data --------------------------------------------------------------------

BD <- readRDS('Data/03_PREDICTSModelData.rds')

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

# 4. Regional Biomes ------------------------------------------------------

# calculate number of studies etc. for each RB
for (i in (1:nrow(ecoregs))) {
  rb = ecoregs$RB_tnc[i]
  tnc_code = paste(ecoregs$RB_tnc[i])
  ecoregs$n_studies[i] = n_distinct(BD$Study_name[BD$RB_tnc==rb])
  ecoregs$n_sources[i] = n_distinct(BD$Source_ID[BD$RB_tnc==rb])
  ecoregs$n_sites[i] = n_distinct(BD$SSBS[BD$RB_tnc==rb]) 
  ecoregs$n_ecoregions[i] = n_distinct(BD$Ecoregion[BD$RB_tnc==rb])
  ecoregs$n_taxa[i] = n_distinct(BD$Taxon[BD$RB_tnc==rb])
  ecoregs$n_country[i] =n_distinct(BD$Country[BD$RB_tnc==rb])
  ecoregs$t_ecoregions[i] = n_distinct(tnc_ecoregions@data$ECO_NAME[tnc_ecoregions@data$RealmMHT==tnc_code])
  ecoregs$p_ecoregions[i] = round(ecoregs$n_ecoregions[i]/ecoregs$t_ecoregions[i], 2)
}

#write.csv(ecoregs, 'FinalScriptsAndData/Data/04_RBsummary.csv', row.names = F)

#That ends up looking like this:

# k1 <- kbl(ecoregs[,c(5,4,1,6:10)], escape = F) %>%
#   kable_classic(full_width= F) %>%
#   kable_styling(fixed_thead = T) %>%
#   column_spec(2:9, background = ifelse(ecoregs$n_studies == 0, "#666", 'white')) %>%
#   collapse_rows(columns = 1, valign = "top")
# 
# save_kable(k1,'FinalScriptsAndData/testkable.pdf')

#taxa is not included in this as it is just geographic representation
###TABLE 1
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

formatflext <- function(ftable) {
  ftable <- set_header_df(ftable, mapping = typology, key = 'col_keys')
  ftable <- merge_h(ftable, part = "header")
  ftable <- merge_v(ftable, part = "header")
  ftable <- align(ftable, i = 1:2, part = 'header', align = 'left')
  ftable <- theme_vanilla(ftable)
}

f1 <- formatflext(f1)

f1 <- merge_v(f1, j = ~ Biome)
f1 <- fix_border_issues(f1)
#conditional formatting
f1 <- bg(f1, bg = '#cccccc', i = ~ n_studies < 1, j = 2:8)
#save_as_docx(f1, path = 'FinalScriptsAndData/testtable.docx')
save_as_image(f1, path = 'FinalScriptsAndData/Figs/01_RBSamplingEffort.png')


# 5. Land Use ----------------------------------------------------------------

LU <- BD %>% count(RB_tnc, LandUse) %>% 
  spread(LandUse, n)

ecoregs<- merge(ecoregs, LU, by.x = 'RB_tnc', all.x = TRUE)

#create agriculture count

ecoregs$Agriculture <- rowSums(ecoregs[,16:18], na.rm = T)


# Taxon Sampling Effort ---------------------------------------------------

# representation by taxa
levels(BD$Study_common_taxon)
n_distinct(BD$Study_common_taxon)
n_distinct(BD$Taxon)

table(BD$Rank_of_study_common_taxon)

# there are a lot of rows that have no rank - these also have no common taxon - these all have species names though - they are mainly plants (SORTED USING COALSESCE FUNCTION AT BEGINNING OF SCRIPT)
table(BD$CommonTaxon)
table(BD$CommonTaxon_Phylum)

##How are taxa spread across regional biomes?
# Across biomes
d <- as.data.frame(table(BD$Biome, BD$CommonTaxon))
d <- spread(d, Var2,Freq)
names(d)[1]<-'Biome'
flextable(d)
# fungi and herptiles mainly underrepresented. Mammals are less studied in grassland biomes. Birds and plants are very well sampled.

# Across realms
e <- as.data.frame(table(BD$Realm, BD$CommonTaxon))
e <- spread(e, Var2, Freq)
names(e)[1]<-'Realm'
flextable(e)

#fungi are not well studied anywhere tropical. 
# no mammal studies in nearctic realm. 
#everywhere else pretty well spread out.

Taxa <- BD %>% count(RB_tnc, CommonTaxon) %>% spread(CommonTaxon, n)
ecoregs<- merge(ecoregs, Taxa, by.x = 'RB_tnc', all.x = TRUE)
ecoregs$Vertebrate <- rowSums(ecoregs[,c(22,24,25)], na.rm = T)

f4 <- flextable(ecoregs[,])

f4 <- f4[, c(4,2,3,1, 5:10)]
f <- arrange(f, Biome_num)
f4 <- flextable(f[,2:10])
f4 <- merge_v(f4, j = ~ Biome)

f.2 <- tidyr::gather(f, 'Common Taxa', 'n', 5:10)


typology <- data.frame(
  col_keys = f4$col_keys,
  type = c("Regional Biome", "Regional Biome",
           "Regional Biome", 'Number of samples of each taxa', 
           'Number of samples of each taxa', 
           'Number of samples of each taxa', 
           'Number of samples of each taxa',
           'Number of samples of each taxa', 
           'Number of samples of each taxa'),
  what = c("Biome", "Realm", "Code", 
           'Plant', "Fungi",  
           'Herptile', 'Invertebrate', 
           'Bird', 'Mammal'),
  stringsAsFactors = FALSE )

f4 <- set_header_df(f4, mapping = typology, key = 'col_keys')
f4 <- merge_h(f4, part = "header")
f4 <- merge_v(f4, part = "header")
f4 <- align(f4, i = 1:2, part = 'header', align = 'left')
f4 <- theme_vanilla(f4)
f4 <- fix_border_issues(f4)
f4
save_as_image(f4, 'FinalScriptsAndData/Figs/03_TaxaSamplingEffort.png')

#Some regional biomes are only represented by one or two taxa groups, meaning they probably aren't representative of the RB as a whole. 
# --we have to assume that the representation of taxa are proportional to the number of taxa in that RB. 

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
save_as_image(f5, 'FinalScriptsAndData/Figs/04_FinalList.png')



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
write.csv(ecoregs, "FinalScriptsAndData/Data/04_RBsummary.csv", row.names = F)

#for appendix 

ap.ecos <- ecoregs[,c(1:6,8,9,11,12,13)]
write.csv(ap.ecos, "FinalScriptsAndData/Figs/RBcoverage.csv", row.names = F)

ap.ecos[ap.ecos$n_sites == 0,]



