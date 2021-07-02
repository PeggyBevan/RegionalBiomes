# Functions that return TraitDatasets

.TrophicLevels<-function(){
  return(c("Herbivore","Carnivore","Fungivore","Detritivore","Omnivore","Non-feeding"))
}

.ReadCategoricalTrait <- function(raw.data, column, unit='Category', 
                                  ranks='Species', levels=NULL, ...) {
  # Helper that processes categorical values in a source data.frame

  # Not all taxa will have a value in the column
  raw.data <- droplevels(raw.data[''!=raw.data[,column],])

  # Some species appear more than once - ensure values in column are unique
  values <- tapply(as.character(Capitalize(StripWhitespace(raw.data[,column]))), 
                   raw.data$Taxon, unique)
  ok <- 1==sapply(values, length)

  if(!all(ok)) {
    stop(paste('Multiple values of ', column, 'for',
               paste(names(ok)[!ok], collapse=', ')))
  } else {
    if(!is.null(levels)) {
      values <- EnsureFactorLevels(values, levels)
    }
    return (TraitDataset(taxa=names(values),
                         values=unname(values), 
                         ranks=ranks, 
                         unit=unit,
                         ...))
    return (trait)
  }
}

.InvertebrateLengthToMass <- function(lengths){
  
  # Currently uses a single allometry for all invertebrates from
  # Sabo et al. (2002) Journal of the North American Benthological Society
  # 21: 336-343.
  # This paper used data for 302 terrestrial invertebrate taxa, including
  # the orders Araneae, Coleoptera, Homoptera, Hymenoptera, Orthoptera,
  # Microcoryphia and Diptera.
  # We use here the overall relationship for all terrestrial invertebrates.
  # This relationship had an R2 of 0.81, and likely masks variation within
  # invertebrate clades.
  
  # Use the published allometry to estimate dry mass
  dry.mass <- 0.032 * lengths^2.63
  
  # Data from the same paper on average proportion of mass that is water
  av.water.content <- 0.925
  
  # Use this value to calculate wet mass
  wet.mass <- dry.mass / (1 - av.water.content)
  
  # Convert from mg to g
  wet.mass <- wet.mass / 1000
  
  return(wet.mass)
  
}

.AmphibianLengthToMass <- function(lengths){
  
  data(svlmass)
  
  model <- lm(log(Mass_g)~log(SVL_mm),data=svlmass)
  
  m0 <- exp(model$coefficients['(Intercept)'])
  b <- model$coefficients['log(SVL_mm)']
  
  masses <- m0 * lengths^b
  
  return(masses)
  
}

.ReptileLengthToMass <- function(lengths_cm,families){
  
  data(reptilesvlmass)
  
  mod <- lm(log(adult_body_mass_g)~log(adult_svl_cm)+family+
              log(adult_svl_cm):family,data=reptilesvlmass)
  
  
  nd <- data.frame(adult_svl_cm=lengths_cm,family=families)
  
  nd$good_family <- nd$family %in% reptilesvlmass$family
  nd$mass <- NA
  
  nd2 <- nd[nd$good_family,]
  
  nd$mass[nd$good_family] <- suppressWarnings(predict(mod,newdata=nd2))
  
  return(exp(nd$mass))
  
}

ButchartAvianBodyMass <- function(path, ...) {
  raw.data <- read.csv(path, ...)
  return (TraitDataset(reference='Butchart',
                       trait='Adult wet mass',
                       kingdom='Animalia',
                       unit='log10 (g)',
                       taxa=raw.data$Sci.name, 
                       ranks='Species',
                       values=log10(raw.data$Overall.mean.mass)))
}

ButchartGenerationLength <- function(path, ...){
  raw.data <- read.csv(path, ...)
  
  raw.data <- raw.data[(raw.data$Binomial!=""),]
  
  return(TraitDataset(reference = 'Butchart',
                      trait = 'Generation length',
                      kingdoms = 'Animalia',
                      unit = 'log10 (years)',
                      taxa = raw.data$Binomial,
                      ranks = 'Species',
                      values = log10(raw.data$Generation_length)))
}

PantheriaMammalianBodyMass <- function(path, sep='\t', header=TRUE, ...) {
  warning("Try the new function - Pantheria - it returns more than just body mass")
  return (Pantheria(path, sep, header, ...)$Pantheria_mass)
}

MyhrvoldAmniotes <- function(path, ...){
  
  raw.data <- read.csv(path, ...)
  
  raw.data <- as.data.frame(lapply(raw.data,function(x){
    if ((class(x)=="numeric") | (class(x)=="integer")){
      x[x <= -999] <- NA
    }
    
    return(x)
  }))
  
  raw.data$adult_body_mass_g[is.na(raw.data$adult_body_mass_g)] <- apply(
    X = raw.data[is.na(raw.data$adult_body_mass_g),c(
      'female_body_mass_g','male_body_mass_g')],
    MARGIN = 1,FUN = mean,na.rm = TRUE)
  raw.data$adult_body_mass_g[is.na(raw.data$adult_body_mass_g)] <- raw.data$no_sex_body_mass_g[
    is.na(raw.data$adult_body_mass_g)]
  
  raw.data$adult_svl_cm[is.na(raw.data$adult_svl_cm)] <- apply(
    X = raw.data[is.na(raw.data$adult_svl_cm),c(
      'female_svl_cm','male_svl_cm')],
    MARGIN = 1,FUN = mean,na.rm = TRUE)
  raw.data$adult_svl_cm[is.na(raw.data$adult_svl_cm)] <- raw.data$no_sex_svl_cm[
    is.na(raw.data$adult_svl_cm)]
  
  raw.data$adult_body_mass_g[(is.na(raw.data$adult_body_mass_g)) & 
                               (raw.data$class=="Reptilia")] <- 
    .ReptileLengthToMass(lengths_cm = raw.data$adult_svl_cm[(
      is.na(raw.data$adult_body_mass_g)) & 
        (raw.data$class=="Reptilia")],
    families = droplevels(raw.data$family[(is.na(raw.data$adult_body_mass_g)) & 
                                 (raw.data$class=="Reptilia")]))
  
  maturity.days <- apply(X = raw.data[,c('female_maturity_d','male_maturity_d')],
                         MARGIN = 1,FUN = mean,na.rm = TRUE)
  maturity.days[is.na(maturity.days)] <- raw.data$no_sex_maturity_d[is.na(maturity.days)]
  
  return(list(Adult_wet_mass=TraitDataset(reference='Myhrvold',
                                          trait='Adult wet mass',
                                          kingdom='Animalia',
                                          unit='log10 (g)',
                                          taxa=paste(raw.data$genus,raw.data$species),
                                          ranks='Species',
                                          values=log10(raw.data$adult_body_mass_g)),
              Maturity_time=TraitDataset(reference = 'Myhrvold',
                                         trait = 'Maturity time',
                                         kingdoms = 'Animalia',
                                         unit = 'log10 (days)',
                                         taxa = paste(raw.data$genus,raw.data$species),
                                         ranks = 'Species',
                                         values = log10(maturity.days)),
              Maximum_longevity=TraitDataset(reference = 'Myhrvold',
                                             trait='Maximum longevity',
                                             kingdoms = 'Animalia',
                                             unit = 'log10 (years)',
                                             taxa = paste(raw.data$genus,raw.data$species),
                                             ranks = 'Species',
                                             values = log10(raw.data$maximum_longevity_y))))
  
}

PantheriaMammals <- function(path, sep='\t', ...) {
  raw.data <- read.csv(path, sep=sep, ...)

  # -999 indicates no data
  mass <- raw.data[,c('MSW05_Binomial', 'X5.1_AdultBodyMass_g')]
  mass <- mass[mass$X5.1_AdultBodyMass_g > 0,]

  habitat <- raw.data[,c('MSW05_Binomial', 'X12.1_HabitatBreadth')]
  habitat <- habitat[habitat$X12.1_HabitatBreadth > 0,]

  diet <- raw.data[,c('MSW05_Binomial', 'X6.1_DietBreadth')]
  diet <- diet[diet$X6.1_DietBreadth > 0,]

  tl <- raw.data[,c('MSW05_Binomial', 'X6.2_TrophicLevel')]
  tl <- tl[tl$X6.2_TrophicLevel > 0,]
  tl$Trophic_level <- ''
  tl$Trophic_level[1 == tl$X6.2_TrophicLevel] <- 'Herbivore'
  tl$Trophic_level[2 == tl$X6.2_TrophicLevel] <- 'Omnivore'
  tl$Trophic_level[3 == tl$X6.2_TrophicLevel] <- 'Carnivore'
  levels <- .TrophicLevels()
  
  tl$Trophic_level <- EnsureFactorLevels(tl$Trophic_level, levels=levels)

  return (list(Adult_wet_mass=TraitDataset(reference='Pantheria',
                                           trait='Adult wet mass',
                                           kingdom='Animalia',
                                           unit='log10 (g)',
                                           taxa=mass$MSW05_Binomial,
                                           ranks='Species',
                                           values=log10(mass$X5.1_AdultBodyMass_g)),
               Habitat_breadth=TraitDataset(reference='Pantheria',
                                            trait='Habitat breadth',
                                            kingdom='Animalia',
                                            unit='Count',
                                            taxa=habitat$MSW05_Binomial,
                                            ranks='Species',
                                            values=habitat$X12.1_HabitatBreadth),
               Diet_breadth=TraitDataset(reference='Pantheria',
                                            trait='Diet breadth',
                                            kingdom='Animalia',
                                            unit='Count',
                                            taxa=diet$MSW05_Binomial,
                                            ranks='Species',
                                            values=diet$X6.1_DietBreadth),
               Trophic_level=TraitDataset(reference='Pantheria',
                                            trait='Trophic level',
                                            kingdom='Animalia',
                                            unit='Category',
                                            taxa=tl$MSW05_Binomial,
                                            ranks='Species',
                                            values=tl$Trophic_level)))
}

TRY <- function(path, sep='\t', fileEncoding='latin1', ...) {
  # TODO col.classes?
  .Log('Reading TRY raw data\n')
  raw.data <- read.csv(path, sep=sep, fileEncoding=fileEncoding, ...)

  # Subset TRY data by data which relates to plant.traits, excluding unused 
  # levels/factors and any accessory data (e.g. latitude, humidity)
  raw.data <- raw.data['x'==raw.data$TraitVsNonTrait,]

  # Just the traits that we are interested in
  traits <- c('Seed mass', 'Plant height vegetative', 'Plant height generative')
  raw.data <- raw.data[raw.data$TraitName %in% traits,]

  # Remove duplicates and NA

  # We combine several data release files so might have duplicated rows.
  bad <- duplicated(raw.data)

  # From the TRY release notes:
  # The data may contain duplicates, if the same data have been contributed to 
  # TRY from different contributors. If we have identified an entry as 
  # duplicate you will find the ID of the original entry in OrigObsDataID.
  bad <- bad | raw.data$OrigObsDataID %in% raw.data$ObsDataID

  # Remove values that are NA, Inf or zero
  bad <- bad | is.na(raw.data$StdValue) | is.infinite(raw.data$StdValue) | 0==raw.data$StdValue
  if(any(bad)) {
    .Log('Removing', sum(bad), 'duplicated, NA, Inf or 0 TRY measurements\n')
    raw.data <- raw.data[!bad,]
  }

  # Just the columns we are interested in
  raw.data <- raw.data[,c('TraitName', 'StdValue', 'AccSpeciesName')]

  .Log('Dropping levels\n')
  raw.data <- droplevels(raw.data)

  # Compute mean of log10
  raw.data$Log10StdValue <- log10(raw.data$StdValue)

  # Sanity check
  stopifnot(!any(is.na(raw.data$Log10StdValue) | is.infinite(raw.data$Log10StdValue)))

  F <- function(trait.name) {
    # Compute species-level mean log10 trait value

    # Just the trait that we are interested in
    trait <- droplevels(raw.data[trait.name==raw.data$TraitName,])

    # tapply gives us an array - elements are species-level mean trait - names 
    # are species names. The steps taken above guarantee that all values of 
    # Log10StdValue will be non-NA.
    return (tapply(trait$Log10StdValue, trait$AccSpeciesName, mean))
  }

  .Log('Computing species-level means\n')
  seed.mass <- F('Seed mass')
  veg.height <- F('Plant height vegetative')
  gen.height <- F('Plant height generative')

  return (list(Seed_mass=TraitDataset(reference='TRY',
                                      trait='Seed mass',
                                      kingdom='Plantae',
                                      unit='log10 (g)',
                                      taxa=names(seed.mass), 
                                      ranks='Species',
                                      values=as.vector(seed.mass)), 
                Vegetative_height=TraitDataset(reference='TRY',
                                       trait='Vegetative height',
                                       kingdom='Plantae',
                                       unit='log10 (m)',
                                       taxa=names(veg.height),
                                       ranks='Species',
                                       values=as.vector(veg.height)), 
               Generative_height=TraitDataset(reference='TRY',
                                       trait='Generative height',
                                       kingdom='Plantae',
                                       unit='log10 (m)',
                                       taxa=names(gen.height), 
                                       ranks='Species',
                                       values=as.vector(gen.height))))
}

CooperHerptileSVL <- function(path, ...) {
  raw.data <- read.csv(path, ...)
  mass.log10 <- log10(.AmphibianLengthToMass(raw.data$svl))
  return (list(LengthDerivedVolume=TraitDataset(reference='Cooper et al 2008',
                                                trait='Length-derived volume',
                                                kingdom='Animalia',
                                                unit='3*log10 (mm)',
                                                taxa=raw.data$binomial, 
                                                ranks='Species',
                                                values=3*log10(raw.data$svl)),
               Adult_wet_mass_inf=TraitDataset(reference='Cooper et al 2008',
                                 trait='Adult wet mass (inferred)',
                                 kingdom='Animalia',
                                 unit='log10 (g)',
                                 taxa=raw.data$binomial,
                                 ranks='Species',
                                 values=mass.log10)))
}

SeniorHerptileSVL <- function(path, ...) {
  raw.data <- read.csv(path, ...)

  # Define ordering of best to worst SVL estimates based on the following logic:
  #
  # 1). Calculated > given because the former is calculated from given male AND 
  #     female data, whilst the latter is simply an overall species estimate as
  #     stated in literature (often without the original data from which it 
  #     was calculated)
  # 2). Midrange > mean because it is less sensitive to extreme values
  # 3). Calculated/given > female/male because it is more representative of the 
  #     species overall (females are generally much larger than males) 
  # 4). Mean > min/max because it is more representative
  # 5). Female values are not > males and min is not > max - this is an arbitrary 
  #     ordering based on the assumption that if you had both values you would 
  #     be calculating the midrange instead

  priority <- c('SVL_calc_midrange','SVL_calc_mean','SVL_given_midrange',
                'SVL_given_mean','SVL_given_min','SVL_given_max','SVL_female_midrange', 
                'SVL_male_midrange', 'SVL_female_mean','SVL_male_mean','SVL_female_min',
                'SVL_female_max','SVL_male_min','SVL_male_max')

  raw.data$Best <- NA
  for (c in priority) {
    raw.data$Best[is.na(raw.data$Best)] <- raw.data[,c][is.na(raw.data$Best)]
  }

  raw.data <- na.omit(raw.data[,c('Taxon','Rank','Best')])

  mass.log10 <- log10(.AmphibianLengthToMass(raw.data$Best))
  
  # SeniorData currently has data for 96 species not covered by Cooper et al. 
  # (57 still lacking SVL, 12 of which have a congener in Cooper et al.)

  return (list(LengthDerivedVolume=TraitDataset(reference='Senior',
                                                trait='Length-derived volume',
                                                kingdom='Animalia',
                                                unit='3*log10 (mm)',
                                                taxa=raw.data$Taxon, 
                                                ranks=raw.data$Rank,
                                                values=3*log10(raw.data$Best)),
               Adult_wet_mass_inf=TraitDataset(reference='Senior',
                                 trait='Adult wet mass (inferred)',
                                 kingdom='Animalia',
                                 unit='log10 (g)',
                                 taxa=raw.data$Taxon,
                                 ranks=raw.data$Rank,
                                 values=mass.log10)))
}

BickfordAmphibians <- function(path, ...){
  raw.data <- read.csv(path, ...)
  
  mass.log10 <- log10(.AmphibianLengthToMass(raw.data$size))
  
  return(list(Adult_wet_mass_ord_inf=TraitDataset(reference='Bickford',
                                                  trait = 'Adult wet mass (ordinal, inferred)',
                                                  unit = 'log10 (g)',
                                                  taxa = raw.data$SCI_NAME,
                                                  ranks = 'Species',
                                                  kingdoms = 'Animalia',
                                                  values = mass.log10),
              LengthDerivedVolume=TraitDataset(reference='Bickford',
                                               trait='Length-derived volume',
                                               kingdom='Animalia',
                                               unit='3*log10 (mm)',
                                               taxa=raw.data$SCI_NAME, 
                                               ranks='Species',
                                               values=3*log10(raw.data$size))))
  
}

GBIFGeographicRange <- function(dir) {
  # The number of 1 degree grid cells in which a species occurs in GBIF

  # The GBIF occurrence CSV files that we downloaded earlier
  .Log('Listing GBIF occurrence files\n')
  files <- list.files(dir, pattern=glob2rx('*csv'), recursive=TRUE, 
                      full.names=TRUE)

  .Log('Reading', length(files), 'GBIF files\n')
  records <- lapply(files, function(path) read.csv(path))

  stopifnot(all(sapply(records, nrow)>0))

  # Extract species names from file paths
  kingdom <- basename(dirname(dirname(files)))
  species <- basename(files)
  species <- gsub('.csv', '', species) # Remove .csv extension

  # log10(total area of occurrence in km^2)
  log.range <- sapply(records, function(rows) {
    log10(sum(DegreeCellAreaKM((rows$minLatitude+rows$maxLatitude)/2,1,1)))
  })

  return (TraitDataset(reference='GBIF',
                       trait='Geographic range',
                       kingdom=kingdom,
                       unit='log(10) (square km)',
                       taxa=species, 
                       ranks='Species', 
                       values=log.range))
}

GBIFGeographicRangeO <- function(dir,gridSize=110000) {
  # The number of grid cells of specified size in which a species occurs in GBIF
  
  files<-.GBIFFiles(dir)
  
  records<-.GBIFRecords(files)
  
  # Extract species names from file paths
  kingdom <- basename(dirname(files))
  species <- basename(files)
  species <- gsub('.csv', '', species) # Remove .csv extension
  
  grids<-.GridGBIF(records,gridSize)
  
  log.range<-sapply(grids,function(x){
    
    return(log10((nrow(x) * gridSize^2)/1000/1000))
    
  })
  
  return (TraitDataset(reference='GBIF',
                       trait=paste('Geographic range (occupancy - ',gridSize/1000,' kilometres)',sep=''),
                       kingdom=kingdom,
                       unit='log(10) (square km)',
                       taxa=species, 
                       ranks='Species', 
                       values=log.range))
}

GBIFGeographicRangeE <- function(dir,gridSize=110000) {
  # The number of grid cells of specified size in which a species occurs in GBIF
  
  ls<-.LandSeaMask(gridSize)
  
  files<-.GBIFFiles(dir)
  
  records<-.GBIFRecords(files)
  
  # Extract species names from file paths
  kingdom <- basename(dirname(files))
  species <- basename(files)
  species <- gsub('.csv', '', species) # Remove .csv extension
  
  grids<-.GridGBIF(records,gridSize)
  
  log.range<-sapply(grids,function(x){
    
    x$minLongitudeShift<-x$minLongitude+17367530
    x$minLongitudeShift[(x$minLongitudeShift>17367530)]<-x$minLongitudeShift[
      (x$minLongitudeShift>17367530)]-34735060
    
    # Order by latitude
    x<-x[order(x$minLatitude),]
    # Find the latitude containing the 2.5th and 97.5th percentile
    lower.lat<-x$minLatitude[head(which(cumsum(x$count)>sum(x$count)*0.025),1)]
    upper.lat<-x$minLatitude[tail(which(cumsum(x$count)<sum(x$count)*0.975),1)]
    
    # Clip the latitude percentiles out of the land-sea mask
    ls<-ls[(ls$y>=lower.lat),]
    ls<-ls[(ls$y<=upper.lat),]
    
    if (diff(range(x$minLongitude)) <= diff(range(x$minLongitudeShift))){
      # if(TRUE){
      # Order by longitude
      x<-x[order(x$minLongitude),]
      # Find the longitude containing the 2.5th and 97.5th percentile
      lower.long<-x$minLongitude[tail(which(!(cumsum(x$count)>ceiling(
        sum(x$count)*0.025))),1)]
      upper.long<-x$minLongitude[head(which(!(cumsum(x$count)<=floor(
        sum(x$count)*0.975))),1)]
      
      # Clip the longitude percentiles out of the land-sea mask
      ls<-ls[(ls$x>=lower.long),]
      ls<-ls[(ls$x<=upper.long),]
      
      
    } else {
      # Order by shifted longitude
      x<-x[order(x$minLongitudeShift),]
      # Find the longitude containing the 2.5th and 97.5th percentile
      lower.long<-x$minLongitudeShift[tail(which(!(cumsum(x$count)>ceiling(
        sum(x$count)*0.025))),1)]
      upper.long<-x$minLongitudeShift[head(which(!(cumsum(x$count)<=floor(
        sum(x$count)*0.975))),1)]
      
      # Create shifted x coordinates in the land-sea mask
      ls$xShift<-ls$x+17367530
      ls$xShift[(ls$xShift>17367530)]<-ls$xShift[(ls$xShift>17367530)]-34735060
      
      # Clip the longitude percentiles out of the land-sea mask
      ls<-ls[(ls$xShift>=lower.long),]
      ls<-ls[(ls$xShift<=upper.long),]
      
    }
    
    # Remove sea cells
    ls<-ls[(ls$mask==1),]
    
    return(log10((nrow(ls) * gridSize^2)/1000/1000))
    
  })
  
  return (TraitDataset(reference='GBIF',
                       trait=paste('Geographic range (extent - ',gridSize/1000,' kilometres)',sep=''),
                       kingdom=kingdom,
                       unit='log(10) (square km)',
                       taxa=species, 
                       ranks='Species', 
                       values=log.range))
}

.GBIFFiles<-function(dir){
  .Log('Listing GBIF occurrence files\n')
  files <- list.files(dir, pattern=glob2rx('*csv'), recursive=TRUE, 
                      full.names=TRUE)
  return(files)
}

.GBIFRecords<-function(files){
  .Log('Reading', length(files), 'GBIF files\n')
  records <- lapply(files, function(path) read.csv(path))
  
  stopifnot(all(sapply(records, nrow)>0))
  
  return(records)
}

.GridGBIF<-function(records,gridSize){
  
  wgsCRS <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
  
  all.longs <- seq(from = -17367530,to = 17367530,by = gridSize)
  all.lats <- seq(from = -7342230,to = 7342230,by = gridSize)
  
  grids<-lapply(records,function(x){
    points <- SpatialPointsDataFrame(
      coords = data.frame(x=x$decimallongitude,
                          y=x$decimallatitude),
      data=data.frame(z=rep(0,nrow(x))),proj4string = wgsCRS)
    pointsProj <- spTransform(x = points,CRSobj = behrCRS)
    
    
    minLongitude<-sapply(pointsProj@coords[,'x'],function(l) max(all.longs[(l>=all.longs)]))
    minLongitude[17367530==minLongitude] <- 17367530 - gridSize
    minLatitude<-sapply(pointsProj@coords[,'y'],function(l) max(all.lats[(l>=all.lats)]))
    minLatitude[7342230==minLatitude] <- 7342230 - gridSize
    
    res <- table(minLongitude,minLatitude)
    res <- as.data.frame(res, responseName='count')
    
    # Drop zero counts
    res <- res[0<res$count,]
    rownames(res) <- NULL
    
    # Compute grid cells
    res$minLatitude <- as.numeric(as.character(res$minLatitude))
    res$minLongitude <- as.numeric(as.character(res$minLongitude))
    res$maxLatitude <- gridSize+res$minLatitude
    res$maxLongitude <- gridSize+res$minLongitude
    # stopifnot(all(res$minLatitude >= -7342230))
    # stopifnot(all(res$maxLatitude <= 7342230))
    # stopifnot(all(res$minLongitude >= -17367530))
    # stopifnot(all(res$maxLatitude <= 17367530))
    
    res <- res[order(res$minLatitude, res$maxLatitude, res$minLongitude, res$maxLongitude),]
    res <- res[,c('minLatitude', 'maxLatitude', 'minLongitude', 'maxLongitude', 'count')]
    
    return(res)
  })
  
  return(grids)
}

.LandSeaMask<-function(gridSize){
  .Log('Generating land-sea mask\n')
  behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

  data(LandSea)
  
  LandSea$mask=1
  
  LandSea <- spTransform(x = LandSea,CRSobj = behrCRS)
  
  r <- SpatialGridDataFrame(grid = GridTopology(cellcentre.offset = c(-17367530+(gridSize/2),
                                                                      -7342230+(gridSize/2)),
                                                cellsize = c(gridSize,gridSize),
                                                cells.dim = c(ceil(34735060/gridSize),
                                                              ceil(14684460/gridSize))),
                            data = data.frame(z=rep(0,ceil(34735060/gridSize)*ceil(14684460/gridSize))))
  
  
  r<-raster(r)
  
  rp<-rasterize(LandSea,r,'mask')
  
  ls.df<-data.frame(coordinates(rp)-0.5)
  ls.df$mask<-rp@data@values
  ls.df$mask[is.na(ls.df$mask)]<-0  
  
  return(ls.df)
  
}

GBIFRecordDensity <- function(dir){
  # The number of records in GBIF divided by the number of 1 degree grid cells in which a species occurs
  
  # The GBIF occurrence CSV files that we downloaded earlier
  .Log('Listing GBIF occurrence files\n')
  files <- list.files(dir, pattern=glob2rx('*csv'), recursive=TRUE, 
                      full.names=TRUE)
  
  .Log('Reading', length(files), 'GBIF files\n')
  records <- lapply(files, function(path) read.csv(path))
  
  stopifnot(all(sapply(records, nrow)>0))
  
  # Extract species names from file paths
  kingdom <- basename(dirname(dirname(files)))
  species <- basename(files)
  species <- gsub('.csv', '', species) # Remove .csv extension
  
  # log10(total area of occurrence in km^2)
  num.records <- sapply(records, function(rows) {
    sum(rows$count)
  })
  
  num.cells <- sapply(records, function(rows) {
    length(rows$count)
  })
  
  return (TraitDataset(reference='GBIF',
                       trait='GBIF record density',
                       kingdom=kingdom,
                       unit='log(10) (records per cell)',
                       taxa=species, 
                       ranks='Species', 
                       values=log10(num.records/num.cells)))
}

BirdlifeRangeArea <- function(path){
  
  raw.data<-read.table(path,header=TRUE,sep="\t")
  
  return(TraitDataset(reference='Birdlife International',
                      trait='Range area (EOO)',
                      unit='log10 (square km)',
                      taxa=raw.data$Species,
                      ranks='Species',
                      kingdoms='Animalia',
                      values=log10(raw.data$range_size_m2/(1000*1000))))
  
}

IUCNRangeArea <- function(path){
  
  raw.data<-read.csv(path)
  
  return(TraitDataset(reference='IUCN',
                      trait='Range area (EOO)',
                      unit='log10 (square km)',
                      taxa=raw.data$Species,
                      ranks='Species',
                      kingdoms='Animalia',
                      values=log10(raw.data$range_size_m2/(1000*1000))))
  
}

GBIFOccurrences <- function(dir) {
  warning("This function renamed to GBIFGeographicRange\n")
  return (GBIFGeographicRange(dir))
}

SantiniDensity <- function(path, ...) {
  
  raw.data <- read.csv(path, ...)
  
  return(TraitDataset(reference='Santini',
                      trait='Average local density',
                      unit='log10 (individuals per km2)',
                      taxa=raw.data$Species,
                      ranks='Species',
                      kingdoms='Animalia',
                      values=log10(raw.data$Density_km2)))
  
}

EdgarBeetleLength <- function(path, ...) {
  # Carabid beetle body length
  # TODO Categorical data:
  # stopifnot(all(raw.data$FLIGHT %in% c('yes','no')))
  # flight <- ifelse('yes'==raw.data$FLIGHT, TRUE, FALSE)
  # larval.diet <- StripWhitespace(Capitalize(as.character(raw.data$D.LARV)))
  # adult.diet <- StripWhitespace(Capitalize(as.character(raw.data$D.ADULT)))
  # larval.habitat <- StripWhitespace(Capitalize(as.character(raw.data$H.LARV)))
  # adult.habitat <- StripWhitespace(Capitalize(as.character(raw.data$H.ADULT)))

  raw.data <- read.csv(path, ...)

  taxon <- paste(raw.data$GENUS, raw.data$SPECIES)

  l <- MeanOfRangeOrValue(raw.data$B.LENGTH)

  Adult_tl<-data.frame(taxon=taxon,D.ADULT=raw.data[,'D.ADULT'])
  Adult_tl$adult_TL<-NA
  Adult_tl$adult_TL['carnivorous' == Adult_tl$D.ADULT] <- 'Carnivore'
  Adult_tl$adult_TL['fungivorous' == Adult_tl$D.ADULT] <- 'Fungivore'
  Adult_tl$adult_TL['herbivorous' == Adult_tl$D.ADULT] <- 'Herbivore'
  Adult_tl$adult_TL['non-feeding' == Adult_tl$D.ADULT] <- 'Non-feeding'
  Adult_tl$adult_TL['omnivorous' == Adult_tl$D.ADULT] <- 'Omnivore'
  Adult_tl$adult_TL['saprophagous' == Adult_tl$D.ADULT] <- 'Detritivore'
  Adult_tl$adult_TL['varied' == Adult_tl$D.ADULT] <- NA
  Adult_tl$adult_TL<-factor(Adult_tl$adult_TL)
  Adult_tl$adult_TL[(Adult_tl$adult_TL=='NA')]<-NA
  Adult_tl$adult_TL<-droplevels(Adult_tl$adult_TL)
  
  
  
  Adult_tl <- droplevels(Adult_tl[which(!is.na(Adult_tl$adult_TL)),])
  
  levels <- .TrophicLevels()
  Adult_tl$adult_TL <- EnsureFactorLevels(Adult_tl$adult_TL, levels=levels)
  
  mass <- .InvertebrateLengthToMass(l)
  
  return (list(Length_derived_Volume=TraitDataset(reference='Edgar',
                       trait='Length-derived volume',
                       kingdom='Animalia',
                       unit='3*log10 (mm)',
                       taxa=taxon, 
                       ranks='Species',
                       values=3*log10(l)),
               Length_derived_mass=TraitDataset(reference='Edgar',
                                                trait='Length-derived mass',
                                                kingdom='Animalia',
                                                unit='g',
                                                taxa=taxon,
                                                ranks='Species',
                                                values=mass),
               Trophic_level=TraitDataset(reference = 'Edgar',
                                          trait='Adult Trophic level',
                                          kingdom='Animalia',
                                          unit='Category',
                                          taxa=Adult_tl$taxon, 
                                          ranks='Species',
                                          values=Adult_tl$adult_TL)))
}

IUCNStatus <- function(dir) {
  # IUCN Red List status

  # The CSV files that we downloaded earlier
  .Log('Listing IUCN files\n')
  files <- list.files(dir, pattern=glob2rx('*csv'), recursive=TRUE, 
                      full.names=TRUE)

  .Log('Reading', length(files), 'IUCN files\n')
  records <- lapply(files, function(path) read.csv(path))

  stopifnot(all(sapply(records, nrow)>0))

  # Take the most recent primary row that is not a subspecies
  records <- lapply(records, function(rows) {
      rows <- rows[is.na(rows$Infrarank) & rows$Modified.Year==max(rows$Modified.Year) & 'true'==rows$Primary,]
      return (rows[,c('Scientific.Name','Category')])
    })

  # Might leave us with some empty records
  records <- records[0<sapply(records, nrow)]
  stopifnot(all(1==sapply(records, nrow)))

  records <- do.call('rbind', records)

  # Map from species name to kingdom
  kingdom <- basename(dirname(dirname(files)))
  names(kingdom) <- gsub('.csv', '', basename(files)) # Remove .csv extension

  stopifnot(all(!is.na(match(records$Scientific.Name, names(kingdom)))))
  records$Kingdom <- kingdom[match(records$Scientific.Name, names(kingdom))]

  # Order factor levels alphabetically
  
  records$Category <- EnsureFactorLevels(records$Category, IUCNStatuses())
  return(TraitDataset(reference='IUCN',
                      trait='IUCN Red List',
                      kingdom=records$Kingdom,
                      unit='Status',
                      taxa=records$Scientific.Name, 
                      ranks='Species', 
                      values=records$Category))
}

CITESStatus <- function(path, ...) {
  cites <- read.csv(path, ...)

  cites <- droplevels(cites[cites$Rank %in% c('SPECIES','SUBSPECIES'),])
  levels(cites$Rank) <- c('Species','Infraspecies')
  cites <- cites[,c('Kingdom','Scientific.Name','Listing','Rank')]

  cites <- cites[!duplicated(cites),]

  stopifnot(!any(duplicated(cites$Scientific.Name)))

  cites$Listing <- EnsureFactorLevels(cites$Listing, c('I','I/II','II','III'))
  return (TraitDataset(reference='CITES',
                       trait='CITES status',
                       kingdom=cites$Kingdom,
                       unit='Appendix',
                       taxa=cites$Scientific.Name, 
                       ranks=cites$Rank, 
                       values=cites$Listing))
}

MeiriAndFeldmanReptilianBodyMass <- function(path, ...) {
  raw.data <- read.csv(path, ...)

  raw.data <- raw.data[,c('reptile.database.name','maximum.size..log.g.')]
  # Remove duplicates
  raw.data <- unique(raw.data)

  # Compute ranks from length of name
  ranks <- strsplit(as.character(raw.data$reptile.database.name), ' ')
  ranks <- sapply(ranks, length)
  stopifnot(all(ranks %in% 1:3))
  ranks <- ifelse(1==ranks, 'Genus', ifelse(2==ranks, 'Species', 'Infraspecies'))

  # Values are ln(g) - convert to log10(g) for consistency with other traits
  values <- log10(exp(raw.data$maximum.size..log.g.))
  return (TraitDataset(reference='Meiri and Feldman',
                       trait='Maximum wet mass',
                       kingdom='Animalia',
                       unit='log10 (g)',
                       taxa=raw.data$reptile.database.name, 
                       ranks=ranks,
                       values=values))
}

GilbertHoverflyThoraxVolume <- function(path, ...) {
  raw.data <- read.csv(path, ...)
  return (TraitDataset(reference='Gilbert',
                       trait='Thorax volume',
                       kingdom='Animalia',
                       unit='log10 (mm^3)',
                       taxa=gsub('_', ' ', (raw.data$Binomial)),
                       ranks='Species',
                       values=log10(raw.data$THVOL)))
}

LeitchPlantPloidyAdjustedCValue <- function(path, ...) {
  # Data from Ilia Leitch
  raw.data <- read.csv(path, ...)

  # Ploidy not know for some records
  raw.data <- raw.data[!is.na(raw.data$C_value_pg_adjusted_for_ploidy),]

  # Some records can not be used - see OK_to_use_notes for reasons
  raw.data <- droplevels(raw.data[raw.data$OK_to_use,])

  # Some species appear more than once - compute species-level means
  values <- log10(raw.data$C_value_pg_adjusted_for_ploidy)
  values <- tapply(values, raw.data$Parsed_name, mean)
  return (TraitDataset(reference='Leitch',
                       trait='C-value adjusted for ploidy',
                       kingdom='Plantae',
                       unit='log10 (pg)',
                       taxa=names(values),
                       ranks='Species',
                       values=as.vector(values)))
}

LeitchPlantCValue <- function(path, ...) {
  # Data from Ilia Leitch
  raw.data <- read.csv(path, ...)

  # Some records can not be used - see OK_to_use_notes for reasons
  raw.data <- droplevels(raw.data[raw.data$OK_to_use,])

  # Some species appear more than once - compute species-level means
  values <- log10(raw.data$C_value_pg)
  values <- tapply(values, raw.data$Parsed_name, mean)

  return (TraitDataset(reference='Leitch',
                       trait='C-value',
                       kingdom='Plantae',
                       unit='log10 (pg)',
                       taxa=names(values),
                       ranks='Species',
                       values=as.vector(values)))
}

.HortonArachnidBodyLength <- function(raw.data) {
  # TODO Make this into a more general method
  # Helper that processes Horton body length data
  raw.data <- cbind(raw.data, Male_body_length=MeanOfRangeOrValue(raw.data$M.BODY.LENGTH))
  raw.data <- cbind(raw.data, Female_body_length=MeanOfRangeOrValue(raw.data$F.BODY.LENGTH))
  raw.data <- cbind(raw.data, Body_length=MeanOfRangeOrValue(raw.data$A.BODY.LENGTH))

  # Where Body_length is empty, compute as the mean of Male_body_length and
  # Female_body_length
  missing <- is.na(raw.data$Body_length)
  raw.data$Body_length[missing] <- apply(raw.data[missing,c('Male_body_length', 'Female_body_length')], 1, mean, na.rm=TRUE)

  # Only rows for which all three of M.BODY.LENGTH, F.BODY.LENGTH and
  # A.BODY.LENGTH should have a Body_length of NA
  stopifnot(all(''==raw.data[is.na(raw.data$Body_length),
                c('M.BODY.LENGTH','F.BODY.LENGTH','A.BODY.LENGTH')]))

  # Not all taxa will have a value in the column
  raw.data <- droplevels(raw.data[!is.na(raw.data$Body_length),])

  # Some species appear more than once - compute species-level means
  values <- tapply(3*log10(raw.data$Body_length), 
                   raw.data$Taxon, mean)

  return (TraitDataset(taxa=names(values),
                       ranks='Species',
                       values=as.vector(values), 
                       reference='Horton',
                       trait='Length-derived volume',
                       kingdom='Animalia',
                       unit='3*log10 (mm)'))
}

HortonArachnid <- function(path, ...) {
  # Returns a list of arachnid body length, hunting method and diet collated
  # by Matt Horton
  raw.data <- read.csv(path, ...)

  # Not all records are resolved to species level
  raw.data <- droplevels(raw.data[''!=raw.data$Genus & ''!=raw.data$Species,])

  raw.data$Taxon <- paste(StripWhitespace(raw.data$Genus), 
                          StripWhitespace(raw.data$Species))

  bl <- .HortonArachnidBodyLength(raw.data)
  hm <- .ReadCategoricalTrait(raw.data, column='HUNTING.METHOD',
                              trait='Hunting method',
                              reference='Horton',
                              kingdom='Animalia')
  diet <- .ReadCategoricalTrait(raw.data, column='DIET', trait='Diet type', 
                                reference='Horton',
                                kingdom='Animalia')
  es <- .ReadCategoricalTrait(raw.data, column='GENERALIST.SPECIALIST', 
                              trait='Ecological specialism',
                              reference='Horton',
                              kingdom='Animalia')
  bm <- TraitDataset(reference = 'Horton',trait = 'Length-derived mass',
                     unit = 'g',
                     taxa = bl$Values$Taxon,
                     ranks = 'Species',
                     kingdom='Animalia',
                     values = .InvertebrateLengthToMass(10^(bl$Values$Value/3)))
  
  
  # Check diet and hunting method make sense
  stopifnot(all.equal(hm$Values$Taxon, diet$Values$Taxon))
  check <- table(hm$Values$Value, diet$Values$Value)
  stopifnot(0==check['Hunter','Herbivore'])
  stopifnot(0==check['Phytophagous','Predator'])
  stopifnot(0==check['Web weaver','Herbivore'])
  return (list(Diet=diet,
               Hunting_method=hm,
               Body_length=bl,
               Length_derived_mass=bm,
               Ecological_specialist=es))
}

.SuFormacidaeBodyLength <- function(raw.data) {
  # Helper that reads ant body length data
  raw.data <- cbind(raw.data, Body_length=MeanOfRangeOrValue(raw.data$Body.Length..mm...worker.))

  raw.data <- droplevels(raw.data[!is.na(raw.data$Body_length),])

  # Some species appear more than once - compute species-level means
  values <- tapply(3*log10(raw.data$Body_length), 
                   raw.data$Taxon, mean)

  return (TraitDataset(taxa=names(values),
                       ranks='Species',
                       values=as.vector(values), 
                       reference='Su',
                       trait='Length-derived volume',
                       kingdom='Animalia',
                       unit='3*log10 (mm)'))
}

.SuHeadDimension <- function(raw.data, min.col, max.col, trait) {
  # Helper that reads ant head dimension data
  raw.data <- droplevels(raw.data[!is.na(raw.data[,min.col]) & 
                                  !is.na(raw.data[,max.col]),])

  values <- apply(raw.data[,c(min.col,max.col)], 1, mean)

  # Some species appear more than once - compute species-level means
  values <- tapply(log10(values), raw.data$Taxon, mean)

  return (TraitDataset(taxa=names(values),
                       ranks='Species',
                       values=as.vector(values), 
                       reference='Su',
                       trait=trait,
                       kingdom='Animalia',
                       unit='log10 (mm)'))
}

SuFormicidae <- function(path, ...) {
  raw.data <- read.csv(path, ...)

  # Remove records not at species level
  raw.data <- droplevels(raw.data[!raw.data$Species %in% c('Spp.','Genus'),])

  raw.data$Taxon <- paste(StripWhitespace(raw.data$Genus), 
                          StripWhitespace(raw.data$Species))

  body.length <- .SuFormacidaeBodyLength(raw.data)
  head.width <- .SuHeadDimension(raw.data, 'Head.W.Min', 'Head.W.Max', 'Head width')
  head.length <- .SuHeadDimension(raw.data, 'Head.L.Min', 'Head.L.Max', 'Head length')
  habitat <- .ReadCategoricalTrait(raw.data, 'Habitat',
                                   trait='Habitat', reference='Su', 
                                   kingdom='Animalia')
  #   diet <- .ReadCategoricalTrait(raw.data, 'Diet',
  #                                 trait='Trophic level', reference='Su', 
  #                                 kingdom='Animalia')
  #   
  Diet1<-data.frame(Taxon=paste(raw.data$Genus,raw.data$Species),Diet=raw.data$Diet)
  Diet1<-droplevels(Diet1[(Diet1$Diet != ""),])
  stopifnot(all(tapply(Diet1$Diet,Diet1$Taxon,function(x) return(length(unique(x))))<=1))

  Diet2<-data.frame(Taxon=paste(raw.data$Genus,raw.data$Species),Diet=raw.data$Diet..Generalised.)
  Diet2<-droplevels(Diet2[(Diet2$Diet != ""),])
  Diet2<-unique(Diet2)
  Diet2.sp<-split(Diet2$Diet,Diet2$Taxon)
  Diet2.sp<-unlist(lapply(Diet2.sp,paste,collapse=" & "))
  

  Trophic_level<-data.frame(Taxon=unique(Diet1$Taxon))
  Trophic_level$Diet<-Diet2$Diet[match(Trophic_level$Taxon,Diet2$Taxon)]
  Trophic_level$Diet[(Trophic_level$Diet=="Unknown")]<-NA
  Trophic_level<-droplevels(Trophic_level)
  
  Trophic_level$Trophic_level<-NA
  
  Trophic_level$Trophic_level['Carnivorous' == Trophic_level$Diet] <- 'Carnivore'
  Trophic_level$Trophic_level['Fungus' == Trophic_level$Diet] <- 'Fungivore'
  Trophic_level$Trophic_level['Honeydew' == Trophic_level$Diet] <- 'Herbivore'
  Trophic_level$Trophic_level['Honeydew & Predatory' == Trophic_level$Diet] <- 'Omnivore'
  Trophic_level$Trophic_level['Honeydew & Scavenger' == Trophic_level$Diet] <- 'Omnivore'
  Trophic_level$Trophic_level['Nectar' == Trophic_level$Diet] <- 'Herbivore'
  Trophic_level$Trophic_level['Nectar//Food Bodies' == Trophic_level$Diet] <- 'Herbivore'
  Trophic_level$Trophic_level['Nectar/Food Bodies' == Trophic_level$Diet] <- 'Herbivore'
  Trophic_level$Trophic_level['Nitrogenous Sources' == Trophic_level$Diet] <- 'Detritivore'
  Trophic_level$Trophic_level['Pollen' == Trophic_level$Diet] <- 'Herbivore'
  Trophic_level$Trophic_level['Predator+honeydew' == Trophic_level$Diet] <- 'Omnivore'
  Trophic_level$Trophic_level['Predatory' == Trophic_level$Diet] <- 'Carnivore'
  Trophic_level$Trophic_level['Predatory & Scavenger' == Trophic_level$Diet] <- 'Omnivore'
  Trophic_level$Trophic_level['Predatory + Honeydew' == Trophic_level$Diet] <- 'Omnivore'
  Trophic_level$Trophic_level['Regurgitation' == Trophic_level$Diet] <- 'Omnivore'
  Trophic_level$Trophic_level['Scavenger' == Trophic_level$Diet] <- 'Detritivore'
  Trophic_level$Trophic_level['Scavengers & Seed Gathering' == Trophic_level$Diet] <- 'Omnivore'
  Trophic_level$Trophic_level['Seed Gathering' == Trophic_level$Diet] <- 'Herbivore'
  Trophic_level$Trophic_level['Seeds+Nectar+honeydew' == Trophic_level$Diet] <- 'Herbivore'  
  
  Trophic_level$Trophic_level[(is.na(Trophic_level$Trophic_level))]<-paste(Diet1$Diet[match(
    Trophic_level$Taxon[(is.na(Trophic_level$Trophic_level))],Diet1$Taxon)])
  
  Trophic_level$Trophic_level['Carnivorous' == Trophic_level$Trophic_level] <- 'Carnivore'  
  Trophic_level$Trophic_level['Herbivorous' == Trophic_level$Trophic_level] <- 'Herbivore'  
  Trophic_level$Trophic_level['Omnivorous' == Trophic_level$Trophic_level] <- 'Omnivore'  
  
  Trophic_level$Trophic_level<-factor(Trophic_level$Trophic_level)
  
  levels <- .TrophicLevels()
  Trophic_level$Trophic_level <- EnsureFactorLevels(Trophic_level$Trophic_level, levels=levels)

  bm <- TraitDataset(reference = 'Su',trait = 'Length-derived mass',
                     unit = 'g',
                     taxa = body.length$Values$Taxon,
                     ranks = 'Species',
                     kingdom='Animalia',
                     values = .InvertebrateLengthToMass(10^(body.length$Values$Value/3)))
  
  return (list(Body_length=body.length,
               Length_derived_mass=bm,
               Head_width=head.width, 
               Head_length=head.length, 
               Habitat=habitat, 
               Trophic_level=TraitDataset(reference='SuFormicidae',
                                          trait='Trophic level',
                                          kingdom='Animalia',
                                          unit='Category',
                                          taxa=Trophic_level$Taxon,
                                          ranks='Species',
                                          values=Trophic_level$Trophic_level)))

}

MiddletonWellingHymenoptera <- function(path, ...) {
  raw.data <- read.csv(path, ...)

  # Remove records not at species level
  raw.data <- droplevels(raw.data[''!=raw.data$Species,])

  raw.data$Taxon <- paste(StripWhitespace(raw.data$Genus), 
                          StripWhitespace(raw.data$Species))

  stopifnot(all(!duplicated(raw.data$Taxon)))

  # Forewing length: M/F/range vs average. (for only a few species per family)
  wl <- cbind(MeanOfRangeOrValue(raw.data$Forewing.Length.Male),
              MeanOfRangeOrValue(raw.data$Forewing.Length.Female))
  raw.data$Forewing.length <- apply(wl, 1, mean, na.rm=TRUE)
  stopifnot(all(raw.data$Forewing.length>0 | is.na(raw.data$Forewing.length)))

  # Wingspan: range or average. (for ~85% of species!)
  ws <- cbind(MeanOfRangeOrValue(raw.data$Wingspan.male),
              MeanOfRangeOrValue(raw.data$Wingspan.female))
  raw.data$Wingspan <- apply(ws, 1, mean, na.rm=TRUE)
  stopifnot(all(raw.data$Wingspan>0 | is.na(raw.data$Wingspan)))

  # No of larval foodplant species: count. Lawrence, this will need to be logged for working out CWM.
  summary(raw.data$Approx..Number.of.Larval.Food.plants)
  stopifnot(all(raw.data$Approx..Number.of.Larval.Food.plants>0 |
                is.na(raw.data$Approx..Number.of.Larval.Food.plants)))

  # Specificity: One species, One genus, One family, One order, > 1 order.
  # Lawrence will write code to produce a column that recodes this as 0-4 scale,
  # so an average can be taken. You need to tell Lawrence what number each level
  # of the factor should get (e.g., 0 = >1 order, 1 = 1 order, and so on)
  raw.data$Larval.Specificity <- Capitalize(StripWhitespace(as.character(raw.data$Larval.Specificity)))
  specificity <- c('More than one order',
                   'One order',
                   'One family',
                   'One genus',
                   'One species')
  stopifnot(all(raw.data$Larval.Specificity %in% c('',specificity)))
  # Empty values do not appear in levels and so will be assigned a value of NA
  raw.data$Larval.Specificity <- factor(raw.data$Larval.Specificity, levels=specificity)

  # Voltinism: Univoltine, Bivoltine, Multivoltine, Highly variable. 
  # Lots of univoltine, lots of multivoltine. Univoltine vs rest is one 
  # natural way to score this; that can be done in the R adapter code, so 
  # Lawrence will work out the binary variable if you tell him which of your 
  # categories need to be combined into each of the new binary categories.
  raw.data$Voltinism <- as.factor(StripWhitespace(tolower(raw.data$Voltinism)))
  summary(raw.data$Voltinism)
  raw.data$Is.univoltine <- 'univoltine'==raw.data$Voltinism
  # Missing values
  raw.data$Is.univoltine[''==raw.data$Voltinism] <- NA
  summary(raw.data$Is.univoltine)

  # Flight months: as a count of months or as a range of counts of months 
  # (so, e.g., "3-7" means "from 3 to 7 months", not "from month 3 to month 7").
  raw.data$Flight.duration <- DiffOfRangeOrValue(raw.data$Flight.Months..number.)
  stopifnot(all((raw.data$Flight.duration>=1 & raw.data$Flight.duration<=12) | 
                 is.na(raw.data$Flight.duration)))

  # Habitat: comma-separated integers corresponding to IUCN habitats. 
  # One useful summary might be a count of the natural major habitats. 
  # another might be whether it is recorded from artificial habitats (14 and 15). 
  # Have these data for maybe as many as 500 species, mostly for butterflies. 
  # Probably the most sensible solution is to have some yarg code for handling 
  # IUCN habitat classifications. Lawrence will write this.
  # http://www.iucnredlist.org/technical-documents/classification-schemes/habitats-classification-scheme-ver3
  raw.data$Artificial <- sapply(strsplit(as.character(raw.data$Habitat), ','),
                                function(v) any(substr(v, 1, 2) %in% c('14','15')))
  raw.data$Artificial[0==nchar(as.character(raw.data$Habitat))] <- NA

  return(list(Flight_duration=TraitDataset(taxa=raw.data$Taxon,
                                           ranks='Species',
                                           values=raw.data$Flight.duration,
                                           reference='Middleton',
                                           trait='Flight duration',
                                           kingdom='Animalia',
                                           unit='months'),
              Wingspan=TraitDataset(taxa=raw.data$Taxon,
                                    ranks='Species',
                                    values=log10(raw.data$Wingspan),
                                    reference='Middleton',
                                    trait='Wingspan',
                                    kingdom='Animalia',
                                    unit='log10 (mm)'),
              Forewing_length=TraitDataset(taxa=raw.data$Taxon,
                                           ranks='Species',
                                           values=log10(raw.data$Forewing.length),
                                           reference='Middleton',
                                           trait='Forewing length',
                                           kingdom='Animalia',
                                           unit='log10 (mm)'),
              Host_diversity=TraitDataset(taxa=raw.data$Taxon,
                                           ranks='Species',
                                           values=log10(raw.data$Approx..Number.of.Larval.Food.plants),
                                           reference='Middleton',
                                           trait='Host diversity',
                                           kingdom='Animalia',
                                           unit='log10 (N hosts)'),
              Host_specificity=TraitDataset(taxa=raw.data$Taxon,
                                           ranks='Species',
                                           values=raw.data$Larval.Specificity,
                                           reference='Middleton',
                                           trait='Host specificity',
                                           kingdom='Animalia',
                                           unit='Category'),
              Univoltine=TraitDataset(taxa=raw.data$Taxon,
                                      ranks='Species',
                                      values=raw.data$Is.univoltine,
                                      reference='Middleton',
                                      trait='Univoltine',
                                      kingdom='Animalia',
                                      unit='Logical'),
              Found_in_artificial_habitat=TraitDataset(taxa=raw.data$Taxon,
                                                       ranks='Species',
                                                       values=raw.data$Artificial,
                                                       reference='Middleton',
                                                       trait='Found in artificial habitat',
                                                       kingdom='Animalia',
                                                       unit='Logical')))
}

KisslingMammals <- function(path, sep='\t', ...) {
  raw.data <- read.csv(path, sep=sep, ...)
  raw.data <- droplevels(raw.data['NotAssigned'!=raw.data$TrophicLevel,])

  levels <- .TrophicLevels()
  raw.data$TrophicLevel <- EnsureFactorLevels(raw.data$TrophicLevel, levels=levels)

  .Log('TODO Other traits from Kissling\n')

  return (list(Trophic_level=TraitDataset(reference='Kissling',
                                           trait='Trophic level',
                                           kingdom='Animalia',
                                           unit='Category',
                                           taxa=paste(raw.data$Genus, raw.data$Species),
                                           ranks='Species',
                                           values=raw.data$TrophicLevel)))
}

AmphibBIO <- function(path, ...) {
  
  raw.data <- read.csv(path, ...)
  
  mass.log10 <- log10(.AmphibianLengthToMass(raw.data$Body_size_mm))
  
  maturity.av <- apply(X = raw.data[,c(
    'Age_at_maturity_min_y','Age_at_maturity_max_y')],
    MARGIN = 1,mean,na.rm=FALSE)
  
  maturity.days <- maturity.av * 365
  
  diet.cols <- c('Arthro'='Diet_invert',
                 'Vert'='Diet_vert_gen',
                 'Fruits'='Diet_fruit',
                 'Seeds'='Diet_seeds',
                 'Flowers'='Diet_flowers',
                 'Leaves'='Diet_plant')
  
  raw.data[,names(diet.cols)][is.na(raw.data[,names(diet.cols)])] <- 0
  
  diet.bad <- 0==rowSums(raw.data[,names(diet.cols)])
  # Three rows have a value of zero in all of cols
  
  diets <- lapply(names(diet.cols), function(col) {
    print(col)
    TraitDataset(reference='AmphiBIO',
                 trait=diet.cols[col],
                 kingdom='Animalia',
                 unit='Yes-No',
                 taxa=raw.data$Species[!diet.bad],
                 ranks='Species',
                 values=raw.data[!diet.bad,col])
  })
  names(diets) <- unname(diet.cols)
  
  # A diet category derived from the ten diet columns
  diet <- rep('', nrow(raw.data))
  plant.cols <- c('Fruits', 'Seeds', 'Flowers','Leaves')
  animal.cols <- c('Vert', 'Arthro')
  all.cols<-c('Fruits', 'Seeds', 'Leaves', 'Flowers','Vert', 'Arthro')
  
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Fruits==1] <- 'Frugivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Seeds==1] <- 'Granivore'
  diet[1<rowSums(raw.data[,plant.cols]) & rowSums(raw.data[,animal.cols])==0] <- 'Herbivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Flowers==1]<-'Herbivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Leaves==1]<-'Herbivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Arthro==1] <- 'Invertivore'
  diet[rowSums(raw.data[,animal.cols])>=1 & rowSums(raw.data[,plant.cols])==0] <- 'Carnivore'
  
  
  # Anything not fitting one of the above definitins is an omnivore
  diet['' == diet] <- 'Omnivore'
  
  stopifnot(all('' != diet))
  diet <- EnsureFactorLevels(diet, levels=c(
    'Frugivore', 'Granivore', 'Herbivore', 'Plant',
    'Omnivore', 'Invertivore', 'Carnivore', 'Animal'))
  
  
  trophic_level <- paste(diet)
  trophic_level[(trophic_level=="Frugivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Granivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Herbivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Omnivore")] <- "Omnivore"
  trophic_level[(trophic_level=="Invertivore")] <- "Carnivore"
  trophic_level[(trophic_level=="Carnivore")] <- "Carnivore"
  trophic_level[(trophic_level=="NA")]<-NA
  trophic_level <- factor(trophic_level)
  
  
  specialisation <- rep('', nrow(raw.data)) 
  specialisation<- rowSums(raw.data[,all.cols]>0)
  specialisation[0==specialisation]<-NA
  
  
  
  
  
  return(list(LengthDerivedVolume=TraitDataset(reference='AmphiBIO',
                                               trait='Length-derived volume',
                                               kingdom='Animalia',
                                               unit='3*log10 (mm)',
                                               taxa=raw.data$Species, 
                                               ranks='Species',
                                               values=3*log10(raw.data$Body_size_mm)),
              Adult_wet_mass_inf=TraitDataset(reference='AmphiBIO',
                                              trait='Adult wet mass (inferred)',
                                              kingdom='Animalia',
                                              unit='log10 (g)',
                                              taxa=raw.data$Species,
                                              ranks='Species',
                                              values=mass.log10),
              Maturity_time=TraitDataset(reference = 'AmphiBIO',
                                         trait = 'Maturity time',
                                         kingdoms = 'Animalia',
                                         unit = 'log10 (days)',
                                         taxa = raw.data$Species,
                                         ranks = 'Species',
                                         values = log10(maturity.days)),
              Maximum_longevity=TraitDataset(reference = 'AmphiBIO',
                                             trait='Maximum longevity',
                                             kingdoms = 'Animalia',
                                             unit = 'log10 (years)',
                                             taxa = raw.data$Species,
                                             ranks = 'Species',
                                             values = log10(raw.data$Longevity_max_y)),
              Diet=TraitDataset(reference='AmphiBIO',
                                trait='Diet',
                                kingdom='Animalia',
                                unit='Category',
                                taxa=raw.data$Species,
                                ranks='Species',
                                values=diet),
              Trophic_Level=TraitDataset(reference='AmphiBIO',
                                         trait='Trophic Level',
                                         kingdom='Animalia',
                                         unit='Category',
                                         taxa=raw.data$Species,
                                         ranks='Species',
                                         values=trophic_level),
              Specialisation=TraitDataset(reference='AmphiBIO',
                                          trait='Specialisation',
                                          kingdom='Animalia',
                                          unit='Number of categories',
                                          taxa=raw.data$Species,
                                          ranks='Species',
                                          values=specialisation)))
  
  
}

MorrillReptileDiet <- function(path, ...){
  
  raw.data <- read.csv(path, ...)
  
  trophic.level <- gsub(" ","",paste(raw.data$TrophicLevel))
  trophic.level[trophic.level=="Carnivorous"] <- "Carnivore"
  trophic.level[trophic.level=="Herbivorous"] <- "Herbivore"
  trophic.level[trophic.level=="Omnivorous"] <- "Omnivore"
  trophic.level[trophic.level=="NA"] <- NA
  
  return(list(Trophic_level=TraitDataset(reference='Morrill',
                                         trait='Trophic level',
                                         kingdom='Animalia',
                                         unit='Category',
                                         taxa=raw.data$Binomial,
                                         ranks='Species',
                                         values=trophic.level)))
  
}

ScharfReptiles <- function(path, ...){
  
  raw.data <- read.csv(path, ...)
  
  maturity.days <- (raw.data$age.at.maturity..months./12)*365
  
  trophic.level <- paste(raw.data$diet)
  trophic.level[trophic.level=="Carnivorous"] <- "Carnivore"
  trophic.level[trophic.level=="Herbivorous"] <- "Herbivore"
  trophic.level[trophic.level=="Omnivorous"] <- "Omnivore"
  
  return(list(Adult_wet_mass=TraitDataset(reference='Scharf',
                                          trait='Adult wet mass',
                                          kingdom='Animalia',
                                          unit='log10 (g)',
                                          taxa=raw.data$species,
                                          ranks='Species',
                                          values=log10(exp(raw.data$log.mass..g.))),
              Maximum_longevity=TraitDataset(reference = 'Scharf',
                                             trait='Maximum longevity',
                                             kingdoms = 'Animalia',
                                             unit = 'log10 (years)',
                                             taxa = raw.data$species,
                                             ranks = 'Species',
                                             values = log10(raw.data$Longevity..years.)),
              Maturity_time=TraitDataset(reference = 'Scharf',
                                         trait = 'Maturity time',
                                         kingdoms = 'Animalia',
                                         unit = 'log10 (days)',
                                         taxa = raw.data$species,
                                         ranks = 'Species',
                                         values = log10(maturity.days)),
              Trophic_level=TraitDataset(reference='Scharf',
                                         trait='Trophic level',
                                         kingdom='Animalia',
                                         unit='Category',
                                         taxa=raw.data$species,
                                         ranks='Species',
                                         values=trophic.level)))
  
  
}

SekerciogluAvianDiet<-function(path,...){
  raw.data<-read.csv(path,...)
  
  diet<-raw.data[,c('Latin','Primary_Diet')]
  
  tl <- raw.data[,c('Latin','Primary_Diet')]
  tl$Trophic_level <- ''
  tl$Trophic_level['FI' == tl$Primary_Diet] <- 'Carnivore'
  tl$Trophic_level['FISC' == tl$Primary_Diet] <- 'Omnivore'
  tl$Trophic_level['FIVE' == tl$Primary_Diet] <- 'Carnivore'
  tl$Trophic_level['FR' == tl$Primary_Diet] <- 'Herbivore'
  tl$Trophic_level['FRIN' == tl$Primary_Diet] <- 'Omnivore'
  tl$Trophic_level['FRNE' == tl$Primary_Diet] <- 'Herbivore'
  tl$Trophic_level['FRSE' == tl$Primary_Diet] <- 'Herbivore'
  tl$Trophic_level['IN' == tl$Primary_Diet] <- 'Carnivore'
  tl$Trophic_level['INNE' == tl$Primary_Diet] <- 'Omnivore'
  tl$Trophic_level['INSC' == tl$Primary_Diet] <- 'Omnivore'
  tl$Trophic_level['INSE' == tl$Primary_Diet] <- 'Omnivore'
  tl$Trophic_level['INVE' == tl$Primary_Diet] <- 'Carnivore'
  tl$Trophic_level['INVG' == tl$Primary_Diet] <- 'Omnivore'
  tl$Trophic_level['NE' == tl$Primary_Diet] <- 'Herbivore'
  tl$Trophic_level['OM' == tl$Primary_Diet] <- 'Omnivore'
  tl$Trophic_level['PL' == tl$Primary_Diet] <- 'Herbivore'
  tl$Trophic_level['SC' == tl$Primary_Diet] <- 'Detritivore'
  tl$Trophic_level['SCVE' == tl$Primary_Diet] <- 'Omnivore'
  tl$Trophic_level['SE' == tl$Primary_Diet] <- 'Herbivore'
  tl$Trophic_level['VE' == tl$Primary_Diet] <- 'Carnivore'
  tl$Trophic_level['VG' == tl$Primary_Diet] <- 'Herbivore'
  
  levels <- .TrophicLevels()
  tl$Trophic_level <- EnsureFactorLevels(tl$Trophic_level, levels=levels)
  
  
  return(list(Diet=TraitDataset(reference='Sekercioglu',
                                trait='Diet',
                                kingdom='Animalia',
                                unit='Category',
                                taxa=diet$Latin,
                                ranks='Species',
                                values=diet$Primary_Diet),
              Trophic_level=TraitDataset(reference='Sekercioglu',
                                         trait='Trophic level',
                                         kingdom='Animalia',
                                         unit='Category',
                                         taxa=tl$Latin,
                                         ranks='Species',
                                         values=tl$Trophic_level)))
  
}

PacificiMammals <- function(path, ...) {
  raw.data <- read.csv(path, ...)
  raw.data <- raw.data[!duplicated(raw.data$Scientific_name),]
  
  raw.data$Max_longevity_d[(raw.data$Max_longevity_d=="no information")] <- NA
  raw.data <- droplevels(raw.data)
  
  raw.data$Max_longevity_d <- paste(raw.data$Max_longevity_d)
  raw.data$Max_longevity_d[(raw.data$Max_longevity_d=="NA")] <- NA
  
  raw.data$Max_longevity_d <- as.numeric(raw.data$Max_longevity_d)
  
  .Log('TODO Other traits from Pacifici\n')
  return (list(Adult_wet_mass=TraitDataset(reference='Pacifici',
                                           trait='Adult wet mass',
                                           kingdom='Animalia',
                                           unit='log10 (g)',
                                           taxa=raw.data$Scientific_name,
                                           ranks='Species',
                                           values=log10(raw.data$AdultBodyMass_g)),
               Generation_length=TraitDataset(reference='Pacifici',
                                           trait='Generation length',
                                           kingdom='Animalia',
                                           unit='log10 (days)',
                                           taxa=raw.data$Scientific_name,
                                           ranks='Species',
                                           values=log10(raw.data$GenerationLength_d)),
               Maximum_longevity=TraitDataset(reference = 'Pacifici',
                                              trait = 'Maximum longevity',
                                              kingdoms = 'Animalia',
                                              unit = 'log10 (days)',
                                              taxa = raw.data$Scientific_name,
                                              ranks = 'Species',
                                              values = log10(raw.data$Max_longevity_d))))
}

AnAge <- function(path, ...) {
  
  raw.data <- read.csv(path, ...)
  
  raw.data[(raw.data == -9999)] <- NA
  
  raw.data$Binomial <- paste(raw.data$Genus,raw.data$Species)
  
  raw.data$Maturity_days <- apply(raw.data[,c(
    'Female.maturity..days.','Male.maturity..days.')],
    1,mean,na.rm=TRUE)
  
  return(list(Maturity_time=TraitDataset(reference = 'AnAge',
                                    trait = 'Maturity time',
                                    kingdoms = 'Animalia',
                                    unit = 'log10 (days)',
                                    taxa = raw.data$Binomial,
                                    ranks = 'Species',
                                    values = log10(raw.data$Maturity_days)),
              Maximum_longevity=TraitDataset(reference = 'AnAge',
                                             trait='Maximum longevity',
                                             kingdoms = 'Animalia',
                                             unit = 'log10 (years)',
                                             taxa = raw.data$Binomial,
                                             ranks = 'Species',
                                             values = log10(raw.data$Maximum.longevity..yrs.))))
  
}

VertebrateClimateTolerance <- function(pathTmax,pathTmin,pathPrecip,...){
  raw.data.tmax <- read.csv(pathTmax)
  raw.data.tmin <- read.csv(pathTmin)
  raw.data.precip <- read.csv(pathPrecip)
  
  return(list(Max_temperature_mean=TraitDataset(reference="Newbold",
                                                trait="Mean maximum temperature",
                                                kingdom="Animalia",
                                                unit="Degrees C times 10",
                                                taxa=raw.data.tmax$Binomial,
                                                ranks="Species",
                                                values=raw.data.tmax$Mean),
              Max_temperature_sd=TraitDataset(reference="Newbold",
                                              trait="SD of maximum temperature",
                                              kingdom="Animalia",
                                              unit="Degrees C times 10",
                                              taxa=raw.data.tmax$Binomial,
                                              ranks="Species",
                                              values=raw.data.tmax$SD),
              Min_temperature_mean=TraitDataset(reference="Newbold",
                                                trait="Mean minimum temperature",
                                                kingdom="Animalia",
                                                unit="Degrees C times 10",
                                                taxa=raw.data.tmin$Binomial,
                                                ranks="Species",
                                                values=raw.data.tmin$Mean),
              Min_temperature_sd=TraitDataset(reference="Newbold",
                                              trait="SD of minimum temperature",
                                              kingdom="Animalia",
                                              unit="Degrees C times 10",
                                              taxa=raw.data.tmin$Binomial,
                                              ranks="Species",
                                              values=raw.data.tmin$SD),
              Precipitation_mean=TraitDataset(reference="Newbold",
                                              trait="Mean precipitation",
                                              taxa=raw.data.precip$Binomial,
                                              kingdom="Animalia",
                                              unit="mm",
                                              ranks="Species",
                                              values=raw.data.precip$Mean),
              Precipitation_sd=TraitDataset(reference="Newbold",
                                            trait="SD of precipitation",
                                            taxa=raw.data.precip$Binomial,
                                            kingdom="Animalia",
                                            unit="mm",
                                            ranks="Species",
                                            values=raw.data.precip$SD)))
}

WilmanMammals <- function(path, sep='\t', ...) {
  raw.data <- read.csv(path, sep=sep, ...)
  
  # Create a dataset for each column
  # A mapping from raw.data column name to a more descriptive name
  diet.cols <- c('Diet.Inv'='Diet_invert',
                 'Diet.Vend'='Diet_vert_endo',
                 'Diet.Vect'='Diet_vert_ecto',
                 'Diet.Vfish'='Diet_vert_fish',
                 'Diet.Vunk'='Diet_vert_unkown',
                 'Diet.Scav'='Diet_scavenging',
                 'Diet.Fruit'='Diet_fruit',
                 'Diet.Nect'='Diet_nectar',
                 'Diet.Seed'='Diet_seeds',
                 'Diet.PlantO'='Diet_plant_other')
  
  diet.bad <- 100!=rowSums(raw.data[,names(diet.cols)])
  # Three rows have a value of zero in all of cols
  stopifnot(all(c('Myzopoda aurita','Mystacina robusta','Mystacina tuberculata') ==
                  raw.data$Scientific[diet.bad]))
  
  diets <- lapply(names(diet.cols), function(col) {
    print(col)
    TraitDataset(reference='Wilman',
                 trait=diet.cols[col],
                 kingdom='Animalia',
                 unit='Percent',
                 taxa=raw.data$Scientific[!diet.bad],
                 ranks='Species',
                 values=raw.data[!diet.bad,col])
  })
  names(diets) <- unname(diet.cols)
  
  # A diet category derived from the ten diet columns
  diet <- rep('', nrow(raw.data))
  diet[100 == raw.data$Diet.Fruit] <- 'Frugivore'
  diet[100 == raw.data$Diet.Nect] <- 'Nectarivore'
  diet[100 == raw.data$Diet.Seed] <- 'Granivore'
  diet[100 == raw.data$Diet.PlantO] <- 'Herbivore'
  plant.cols <- c('Diet.Fruit', 'Diet.Nect', 'Diet.Seed', 'Diet.PlantO')
  diet[''==diet & 100 == rowSums(raw.data[,plant.cols])] <- 'Plant'
  
  diet[100 == raw.data$Diet.Inv] <- 'Invertivore'
  carnivore.cols <- c('Diet.Vend', 'Diet.Vect', 'Diet.Vfish', 'Diet.Vunk')
  diet[100 == rowSums(raw.data[, carnivore.cols])] <- 'Carnivore'
  animal.cols <- c('Diet.Inv', 'Diet.Vend', 'Diet.Vect', 'Diet.Vfish', 'Diet.Vunk')
  diet[''==diet & 100 == rowSums(raw.data[,animal.cols])] <- 'Animal'
  
  # No data provided for these three species
  diet['Myzopoda aurita' == raw.data$Scientific] <- 'Invertivore'
  diet['Mystacina robusta' == raw.data$Scientific] <- 'Omnivore'
  diet['Myzopoda aurita' == raw.data$Scientific] <- 'Omnivore'
  
  # Anything not fitting one of the above definitins is an omnivore
  diet['' == diet] <- 'Omnivore'
  
  stopifnot(all('' != diet))
  diet <- EnsureFactorLevels(diet, levels=c(
    'Frugivore', 'Nectarivore', 'Granivore', 'Herbivore', 'Plant',
    'Omnivore', 'Invertivore', 'Carnivore', 'Animal'))
  
  
  trophic_level <- paste(diet)
  trophic_level[(trophic_level=="Frugivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Nectarivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Granivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Herbivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Plant")] <- "Herbivore"
  trophic_level[(trophic_level=="Omnivore")] <- "Omnivore"
  trophic_level[(trophic_level=="Invertivore")] <- "Carnivore"
  trophic_level[(trophic_level=="Carnivore")] <- "Carnivore"
  trophic_level[(trophic_level=="Animal")] <- "Carnivore"
  trophic_level <- factor(trophic_level)
  
  
  specialisation <- rep('', nrow(raw.data)) 
  all.cols<-c('Diet.Inv', 'Diet.Vend', 'Diet.Vect', 'Diet.Vfish','Diet.Vunk', 'Diet.Fruit', 'Diet.Nect', 'Diet.Seed', 'Diet.PlantO')
  specialisation<- rowSums(raw.data[,all.cols]>0)
  specialisation[0==specialisation]<-NA
  
  # An activity category derived from the Activity.* columns
  activity <- rep('', nrow(raw.data))
  activity[1 == raw.data$Activity.Nocturnal & 0 == raw.data$Activity.Crepuscular & 0 == raw.data$Activity.Diurnal] <- 'ObiglateNoctural'
  activity[1 == raw.data$Activity.Nocturnal & 1 == raw.data$Activity.Crepuscular & 0 == raw.data$Activity.Diurnal] <- 'NocturnalAndCrepuscular'
  activity[0 == raw.data$Activity.Nocturnal & 0 == raw.data$Activity.Crepuscular & 1 == raw.data$Activity.Diurnal] <- 'ObiglateDiurnal'
  activity[0 == raw.data$Activity.Nocturnal & 1 == raw.data$Activity.Crepuscular & 1 == raw.data$Activity.Diurnal] <- 'DiurnalAndCrepuscular'
  activity['' == activity] <- 'Opportunist'
  stopifnot(all('' != activity))
  
  foraging.mapping <- c(M='Marine',
                        G='Ground level, including aquatic',
                        S='Scansorial',
                        Ar='Arboreal',
                        A='Aerial')
  
  return (c(diets,
            list(Diet=TraitDataset(reference='Wilman',
                                   trait='Diet',
                                   kingdom='Animalia',
                                   unit='Category',
                                   taxa=raw.data$Scientific,
                                   ranks='Species',
                                   values=diet),
                 Trophic_Level=TraitDataset(reference='Wilman',
                                            trait='Trophic Level',
                                            kingdom='Animalia',
                                            unit='Category',
                                            taxa=raw.data$Scientific,
                                            ranks='Species',
                                            values=trophic_level),
                 Specialisation=TraitDataset(reference='Wilman',
                                             trait='Specialisation',
                                             kingdom='Animalia',
                                             unit='Number of categories',
                                             taxa=raw.data$Scientific,
                                             ranks='Species',
                                             values=specialisation),
                 Activity_time=TraitDataset(reference='Wilman',
                                            trait='Activity_time',
                                            kingdom='Animalia',
                                            unit='Category',
                                            taxa=raw.data$Scientific,
                                            ranks='Species',
                                            values=activity),
                 Adult_wet_mass=TraitDataset(reference='Wilman',
                                             trait='Adult wet mass',
                                             kingdom='Animalia',
                                             unit='log10 (g)',
                                             taxa=raw.data$Scientific,
                                             ranks='Species',
                                             values=log10(raw.data$BodyMass.Value)),
                 Is_nocturnal=TraitDataset(reference='Wilman',
                                           trait='Is_nocturnal',
                                           kingdom='Animalia',
                                           unit='Logical',
                                           taxa=raw.data$Scientific,
                                           ranks='Species',
                                           values=raw.data$Activity.Nocturnal>0),
                 Is_crepuscular=TraitDataset(reference='Wilman',
                                             trait='Is_crepuscular',
                                             kingdom='Animalia',
                                             unit='Logical',
                                             taxa=raw.data$Scientific,
                                             ranks='Species',
                                             values=raw.data$Activity.Crepuscular>0),
                 Is_diurnal=TraitDataset(reference='Wilman',
                                         trait='Is_diurnal',
                                         kingdom='Animalia',
                                         unit='Logical',
                                         taxa=raw.data$Scientific,
                                         ranks='Species',
                                         values=raw.data$Activity.Diurnal>0),
                 Foraging_strategy=TraitDataset(reference='Wilman',
                                                trait='Foraging strategy',
                                                kingdom='Animalia',
                                                unit='Category',
                                                taxa=raw.data$Scientific,
                                                ranks='Species',
                                                values=foraging.mapping[raw.data$ForStrat.Value]))))
}

WilmanBirds <- function(path, sep='\t', ...) {
  raw.data <- read.csv(path,sep=sep,...)
  
  raw.data <- droplevels(raw.data[(raw.data$Scientific!=""),])
  
  # Create a dataset for each column
  # A mapping from raw.data column name to a more descriptive name
  diet.cols <- c('Diet.Inv'='Diet_invert',
                 'Diet.Vend'='Diet_vert_endo',
                 'Diet.Vect'='Diet_vert_ecto',
                 'Diet.Vfish'='Diet_vert_fish',
                 'Diet.Vunk'='Diet_vert_unkown',
                 'Diet.Scav'='Diet_scavenging',
                 'Diet.Fruit'='Diet_fruit',
                 'Diet.Nect'='Diet_nectar',
                 'Diet.Seed'='Diet_seeds',
                 'Diet.PlantO'='Diet_plant_other')
  
  diets <- lapply(names(diet.cols), function(col) {
    print(col)
    TraitDataset(reference='Wilman',
                 trait=diet.cols[col],
                 kingdom='Animalia',
                 unit='Percent',
                 taxa=raw.data$Scientific,
                 ranks='Species',
                 values=raw.data[,col])
  })
  names(diets) <- unname(diet.cols)
  
  # A diet category derived from the ten diet columns
  diet <- rep('', nrow(raw.data))
  diet[100 == raw.data$Diet.Fruit] <- 'Frugivore'
  diet[100 == raw.data$Diet.Nect] <- 'Nectarivore'
  diet[100 == raw.data$Diet.Seed] <- 'Granivore'
  diet[100 == raw.data$Diet.PlantO] <- 'Herbivore'
  plant.cols <- c('Diet.Fruit', 'Diet.Nect', 'Diet.Seed', 'Diet.PlantO')
  diet[''==diet & 100 == rowSums(raw.data[,plant.cols])] <- 'Plant'
  
  diet[100 == raw.data$Diet.Inv] <- 'Invertivore'
  carnivore.cols <- c('Diet.Vend', 'Diet.Vect', 'Diet.Vfish', 'Diet.Vunk')
  diet[100 == rowSums(raw.data[, carnivore.cols])] <- 'Carnivore'
  animal.cols <- c('Diet.Inv', 'Diet.Vend', 'Diet.Vect', 'Diet.Vfish', 'Diet.Vunk')
  diet[''==diet & 100 == rowSums(raw.data[,animal.cols])] <- 'Animal'
  
  # Anything not fitting one of the above definitins is an omnivore
  diet['' == diet] <- 'Omnivore'
  
  stopifnot(all('' != diet))
  diet <- EnsureFactorLevels(diet, levels=c(
    'Frugivore', 'Nectarivore', 'Granivore', 'Herbivore', 'Plant',
    'Omnivore', 'Invertivore', 'Carnivore', 'Animal'))
  
  
  
  trophic_level <- paste(diet)
  trophic_level[(trophic_level=="Frugivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Nectarivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Granivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Herbivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Plant")] <- "Herbivore"
  trophic_level[(trophic_level=="Omnivore")] <- "Omnivore"
  trophic_level[(trophic_level=="Invertivore")] <- "Carnivore"
  trophic_level[(trophic_level=="Carnivore")] <- "Carnivore"
  trophic_level[(trophic_level=="Animal")] <- "Carnivore"
  trophic_level <- factor(trophic_level)
  
  
  specialisation <- rep('', nrow(raw.data)) 
  all.cols<-c('Diet.Inv', 'Diet.Vend', 'Diet.Vect', 'Diet.Vfish','Diet.Vunk', 'Diet.Fruit', 'Diet.Nect', 'Diet.Seed', 'Diet.PlantO')
  specialisation<- rowSums(raw.data[,all.cols]>0)
  specialisation[0==specialisation]<-NA
  
  # An activity category derived from the Activity.* columns
  activity <- rep('', nrow(raw.data))
  activity[1 == raw.data$Activity.Nocturnal & 0 == raw.data$Activity.Crepuscular & 0 == raw.data$Activity.Diurnal] <- 'ObiglateNoctural'
  activity[1 == raw.data$Activity.Nocturnal & 1 == raw.data$Activity.Crepuscular & 0 == raw.data$Activity.Diurnal] <- 'NocturnalAndCrepuscular'
  activity[0 == raw.data$Activity.Nocturnal & 0 == raw.data$Activity.Crepuscular & 1 == raw.data$Activity.Diurnal] <- 'ObiglateDiurnal'
  activity[0 == raw.data$Activity.Nocturnal & 1 == raw.data$Activity.Crepuscular & 1 == raw.data$Activity.Diurnal] <- 'DiurnalAndCrepuscular'
  activity['' == activity] <- 'Opportunist'
  stopifnot(all('' != activity))
  
  foraging.mapping <- c(M='Marine',
                        G='Ground level, including aquatic',
                        S='Scansorial',
                        Ar='Arboreal',
                        A='Aerial')
  
  return (c(diets,
            list(Diet=TraitDataset(reference='Wilman',
                                   trait='Diet',
                                   kingdom='Animalia',
                                   unit='Category',
                                   taxa=raw.data$Scientific,
                                   ranks='Species',
                                   values=diet),
                 Trophic_Level=TraitDataset(reference='Wilman',
                                            trait='Trophic level',
                                            kingdom='Animalia',
                                            unit='Category',
                                            taxa=raw.data$Scientific,
                                            ranks='Species',
                                            values=trophic_level),
                 Specialisation=TraitDataset(reference='Wilman',
                                             trait='Specialisation',
                                             kingdom='Animalia',
                                             unit='Number of categories',
                                             taxa=raw.data$Scientific,
                                             ranks='Species',
                                             values=specialisation),
                 Activity_time=TraitDataset(reference='Wilman',
                                            trait='Activity_time',
                                            kingdom='Animalia',
                                            unit='Category',
                                            taxa=raw.data$Scientific,
                                            ranks='Species',
                                            values=activity),
                 Adult_wet_mass=TraitDataset(reference='Wilman',
                                             trait='Adult wet mass',
                                             kingdom='Animalia',
                                             unit='log10 (g)',
                                             taxa=raw.data$Scientific,
                                             ranks='Species',
                                             values=log10(raw.data$BodyMass.Value)),
                 Is_nocturnal=TraitDataset(reference='Wilman',
                                           trait='Is_nocturnal',
                                           kingdom='Animalia',
                                           unit='Logical',
                                           taxa=raw.data$Scientific,
                                           ranks='Species',
                                           values=raw.data$Nocturnal>0))))
}


BentleyArthropods <- function(path) {
  raw.data <- read.csv(path)

  # -999 indicates no data
  Adult_min_L<- raw.data[,c('Family', 'adult_minl')]
  Adult_min_L<- Adult_min_L[which(Adult_min_L$adult_minl > 0),]

  Adult_max_L<- raw.data[,c('Family','adult_maxl')]
  Adult_max_L<- Adult_max_L[which(Adult_max_L$adult_maxl > 0),]
  
  Adult_min_wingspan <- raw.data[,c('Family','Wing_span_min')]
  Adult_min_wingspan <- Adult_min_wingspan[which(Adult_min_wingspan$Wing_span_min>0),]
  
  Adult_max_wingspan <- raw.data[,c('Family','Wing_span_max')]
  Adult_max_wingspan <- Adult_max_wingspan[which(Adult_max_wingspan$Wing_span_max>0),]
  
  Adult_tl <- raw.data[,c('Family','adult_TL')]
  Adult_tl <- droplevels(Adult_tl[which(Adult_tl$adult_TL!=''),])
  Adult_tl$adult_TL<-paste(Adult_tl$adult_TL)
  Adult_tl$adult_TL['H' == Adult_tl$adult_TL] <- 'Herbivore'
  Adult_tl$adult_TL['C' == Adult_tl$adult_TL] <- 'Carnivore'
  Adult_tl$adult_TL['F' == Adult_tl$adult_TL] <- 'Fungivore'
  Adult_tl$adult_TL['D' == Adult_tl$adult_TL] <- 'Detritivore'
  Adult_tl$adult_TL['O.t' == Adult_tl$adult_TL] <- 'Omnivore'
  Adult_tl$adult_TL['N' == Adult_tl$adult_TL] <- 'Non-feeding'
  Adult_tl$adult_TL['O.f' == Adult_tl$adult_TL] <- NA
  Adult_tl$adult_TL<-factor(Adult_tl$adult_TL)
  Adult_tl <- droplevels(Adult_tl[which(!is.na(Adult_tl$adult_TL)),])


  levels <- .TrophicLevels()
  Adult_tl$adult_TL <- EnsureFactorLevels(Adult_tl$adult_TL, levels=levels)

  Klepto_Para<- raw.data[,c('Family','kleptoparasite')]
  Klepto_Para <- droplevels(Klepto_Para[which(Klepto_Para$kleptoparasite!=''),])
  Klepto_Para$kleptoparasite[0== Klepto_Para$kleptoparasite] <- "No kleptoparasites"
  Klepto_Para$kleptoparasite[1== Klepto_Para$kleptoparasite] <- "Some kleptoparasites"
  Klepto_Para$kleptoparasite[2== Klepto_Para$kleptoparasite] <- "All kelptoparasites"
  levels2 <- c('No kleptoparasites','Some kleptoparasites','All kelptoparasites')
  Klepto_Para$kleptoparasite <- EnsureFactorLevels(Klepto_Para$kleptoparasite, levels=levels2)
  
  Parasites<- raw.data[,c('Family','parasite')]
  Parasites <- droplevels(Parasites[which(Parasites$parasite!=''),])
  Parasites$parasite[0== Parasites$parasite] <- "No parasites"
  Parasites$parasite[1== Parasites$parasite] <- "Some parasites"
  Parasites$parasite[2== Parasites$parasite] <- "All parasites"
  levels3 <- c('No parasites','Some parasites','All parasites')
  Parasites$parasite <- EnsureFactorLevels(Parasites$parasite, levels=levels3)
  
  Parasitoids<- raw.data[,c('Family','parasitoid')]
  Parasitoids<- droplevels(Parasitoids[which(Parasitoids$parasitoid!=''),])
  Parasitoids$parasitoid[0== Parasitoids$parasitoid] <- "No parasitoids"
  Parasitoids$parasitoid[1== Parasitoids$parasitoid] <- "Some parasitoids"
  Parasitoids$parasitoid[2== Parasitoids$parasitoid] <- "All parasitoids"
  levels4 <- c('No parasitoids','Some parasitoids','All parasitoids')
  Parasitoids$parasitoid <- EnsureFactorLevels(Parasitoids$parasitoid, levels=levels4)

  Endo_Ecto_Para<- raw.data[,c('Family','endo_ecto')]
  Endo_Ecto_Para<- droplevels(Endo_Ecto_Para[which(Endo_Ecto_Para$endo_ecto!=''),])
  Endo_Ecto_Para$endo_ecto[0== Endo_Ecto_Para$endo_ecto] <- "Neither"
  Endo_Ecto_Para$endo_ecto[1== Endo_Ecto_Para$endo_ecto] <- "Endoparasites"
  Endo_Ecto_Para$endo_ecto[2== Endo_Ecto_Para$endo_ecto] <- "Ectoparasites"
  Endo_Ecto_Para$endo_ecto[3== Endo_Ecto_Para$endo_ecto] <- "Endoparasites and Ectoparasites"
  levels5 <- c('Neither','Endoparasites','Ectoparasites','Endoparasites and Ectoparasites')
  Endo_Ecto_Para$endo_ecto <- EnsureFactorLevels(Endo_Ecto_Para$endo_ecto, levels=levels5)

  Idio_Koino_biont<- raw.data[,c('Family','idio_koino')]
  Idio_Koino_biont <- droplevels(Idio_Koino_biont[which(Idio_Koino_biont$idio_koino!=''),])
  Idio_Koino_biont$idio_koino[0== Idio_Koino_biont$idio_koino] <- "Neither"
  Idio_Koino_biont$idio_koino[1== Idio_Koino_biont$idio_koino] <- "Idiobionts"
  Idio_Koino_biont$idio_koino[2== Idio_Koino_biont$idio_koino] <- "Koinobionts"
  Idio_Koino_biont$idio_koino[3== Idio_Koino_biont$idio_koino] <- "Idiobionts and Koinobionts"
  levels6 <- c('Neither','Idiobionts','Koinobionts','Idiobionts and Koinobionts')
  Idio_Koino_biont$idio_koino <- EnsureFactorLevels(Idio_Koino_biont$idio_koino, levels=levels6)


  Adult_min_L$min_mass <- .InvertebrateLengthToMass(Adult_min_L$adult_minl)
  Adult_max_L$max_mass <- .InvertebrateLengthToMass(Adult_max_L$adult_maxl)

  return (list(Adult_Min_Length=TraitDataset(reference='Bentley_Arthropods',
                                             trait='Adult Minimum Length',
                                             kingdom='Animalia',
                                             unit='mm',
                                             taxa=Adult_min_L$Family,
                                             ranks='Family',
                                             values=Adult_min_L$adult_minl),
               Adult_Max_Length=TraitDataset(reference='Bentley_Arthropods',
                                            trait='Adult Maximum Length',
                                            kingdom='Animalia',
                                            unit='mm',
                                            taxa=Adult_max_L$Family,
                                            ranks='Family',
                                            values=Adult_max_L$adult_maxl),
               Adult_Min_Mass=TraitDataset(reference='Bentley_Arthropods',
                                           trait='Adult Minimum Mass',
                                           kingdom='Animalia',
                                           unit="g",
                                           taxa=Adult_min_L$Family,
                                           ranks='Family',
                                           values=Adult_min_L$min_mass),
               Adult_Max_Mass=TraitDataset(reference='Bentley_Arthropods',
                                           trait='Adult Maximum Mass',
                                           kingdom='Animalia',
                                           unit='g',
                                           taxa=Adult_max_L$Family,
                                           ranks='Family',
                                           values=Adult_max_L$max_mass),
               Adult_Min_Wingspan=TraitDataset(reference='Bentley_Arthropods',
                                               trait='Adult Minimum Wingspan',
                                               kingdom='Animalia',
                                               unit='mm',
                                               taxa=Adult_min_wingspan$Family,
                                               ranks='Family',
                                               values=Adult_min_wingspan$Wing_span_min),
               Adult_Max_Wingspan=TraitDataset(reference='Bentley_Arthropods',
                                               trait='Adult Maximum Wingspan',
                                               kingdom='Animalia',
                                               unit='mm',
                                               taxa=Adult_max_wingspan$Family,
                                               ranks='Family',
                                               values=Adult_max_wingspan$Wing_span_max),
                Adult_Trophic_Level=TraitDataset(reference='Bentley_Arthropods',
                                        trait='Adult Trophic level',
                                        kingdom='Animalia',
                                        unit='Category',
                                        taxa=Adult_tl$Family,
                                        ranks='Family',
                                        values=Adult_tl$adult_TL),
                Kleptoparasites=TraitDataset(reference='Bentley_Arthropods',
                                        trait = 'Kleptoparasites',
                                        kingdom = 'Animalia',
                                        unit='Category',
                                        taxa=Klepto_Para$Family,
                                        ranks='Family',
                                        values=Klepto_Para$kleptoparasite),
                Parasites=TraitDataset(reference='Bentley_Arthropods',
                                            trait = 'Parasites',
                                            kingdom = 'Animalia',
                                            unit='Category',
                                            taxa=Parasites$Family,
                                            ranks='Family',
                                            values=Parasites$parasite),
               Parasitoids=TraitDataset(reference='Bentley_Arthropods',
                                      trait = 'Parasitoids',
                                      kingdom = 'Animalia',
                                      unit='Category',
                                      taxa=Parasitoids$Family,
                                      ranks='Family',
                                      values=Parasitoids$parasitoid),
              Endoparasites_or_Ectoparasites=TraitDataset(reference='Bentley_Arthropods',
                                      trait = 'Endoparasites or Ectoparasites',
                                      kingdom = 'Animalia',
                                      unit='Category',
                                      taxa=Endo_Ecto_Para$Family,
                                      ranks='Family',
                                      values=Endo_Ecto_Para$endo_ecto),
              Idiobiont_or_Koinobiont=TraitDataset(reference='Bentley_Arthropods',
                                      trait = 'Idiobiont or Koinobiont',
                                      kingdom = 'Animalia',
                                      unit='Category',
                                      taxa=Idio_Koino_biont$Family,
                                      ranks='Family',
                                      values=Idio_Koino_biont$idio_koino)))
               
}

BentleyMissingTaxa<- function(path) {
  raw.data <- read.csv(path)

  Adult_tl <- raw.data[,c('Family','adult_TL')]
  Adult_tl <- droplevels(Adult_tl[which(Adult_tl$adult_TL!=''),])
  Adult_tl$adult_TL<-paste(Adult_tl$adult_TL)
  Adult_tl$adult_TL['H' == Adult_tl$adult_TL] <- 'Herbivore'
  Adult_tl$adult_TL['C' == Adult_tl$adult_TL] <- 'Carnivore'
  Adult_tl$adult_TL['F' == Adult_tl$adult_TL] <- 'Fungivore'
  Adult_tl$adult_TL['D' == Adult_tl$adult_TL] <- 'Detritivore'
  Adult_tl$adult_TL['O.t' == Adult_tl$adult_TL] <- 'Omnivore'
  Adult_tl$adult_TL['N' == Adult_tl$adult_TL] <- 'Non-feeding'
  Adult_tl$adult_TL['O.f' == Adult_tl$adult_TL] <- NA
  Adult_tl$adult_TL<-factor(Adult_tl$adult_TL)
  Adult_tl <- droplevels(Adult_tl[which(!is.na(Adult_tl$adult_TL)),])
  
  levels <- .TrophicLevels()
  Adult_tl$adult_TL <- EnsureFactorLevels(Adult_tl$adult_TL, levels=levels)
  
  return (Adult_Trophic_Level=TraitDataset(reference='Bentley_Arthropods',
                                               trait='Adult Trophic level',
                                               kingdom='Animalia',
                                               unit='Category',
                                               taxa=Adult_tl$Family,
                                               ranks='Family',
                                               values=Adult_tl$adult_TL))
  
}

FlynnHymenoptera <- function(path, ...) {
  raw.data <- read.csv(path, ...)
  raw.data$Taxon <- paste(raw.data$Genus, raw.data$Species)

  # Mean body length
  raw.data$Body.Length <- apply(
      raw.data[, c('Body.Length.Lower.Bound', 'Body.Length.Upper.Bound')], 1, mean
  )

  raw.data$Dietary.Breadth <- Capitalize(StripWhitespace(as.character(raw.data$Dietary.Breadth)))
  raw.data$Generalism <- raw.data$Dietary.Breadth
  raw.data$Generalism['No lecty status' == raw.data$Generalism] <- ''

  raw.data$Parasitism <- raw.data$Dietary.Breadth
  raw.data$Parasitism['No lecty status' == raw.data$Parasitism] <- 'Parasitic'
  raw.data$Parasitism[raw.data$Parasitism %in% c('Polylectic', 'Oligolectic')] <- 'Not parasitic'

  # Nesting strategy
  # 1 Burrowing in soil
  # 2a  Gnawing within plants rotten wood
  # 2b  Gnawing within plants soft pithy stems
  # 2c  Gnawing within plants dense wood
  # 3a  Natural cavities
  # 3b  Natural cavities  parasites
  # 4 Exposed surfaces

  # As a logical trait
  # 3a and 3b are natural cavities
  raw.data$Nesting.Strategy.Natural <- raw.data$Nesting.Strategy %in% c('3', '3a', '3b')
  # As a scale: Soil - exposed surface - natural cavities - cavities within plants
  raw.data$Nesting.Strategy.Stratum <- ''
  raw.data$Nesting.Strategy.Stratum[raw.data$Nesting.Strategy == '1'] <- 'Burrowing in soil'
  raw.data$Nesting.Strategy.Stratum[raw.data$Nesting.Strategy == '4'] <- 'Exposed surfaces'
  raw.data$Nesting.Strategy.Stratum[raw.data$Nesting.Strategy %in% c('3', '3a', '3b')] <- 'Natural cavities'
  raw.data$Nesting.Strategy.Stratum[raw.data$Nesting.Strategy %in% c('2a', '2b', '2c')] <- 'Cavities within plants'

  # 1* Obligately Solitary
  # 2* Not Obligately Solitary
  raw.data$Obligately.solitary <- grepl('^1.*', raw.data$Sociality)
  raw.data$Obligately.solitary['' == raw.data$Sociality] <- NA

  return (list(
    Measured_intertegular_distance=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=log10(raw.data$Intertegular.Distance..mm.),
      reference='Flynn',
      trait='Measured intertegular distance',
      kingdom='Animalia',
      unit='log10 (mm)'
    ),
    Generalism=.ReadCategoricalTrait(
      raw.data,
      'Generalism',
      levels=c('Oligolectic', 'Polylectic'),
      reference='Flynn',
      trait='Generalism',
      kingdoms='Animalia'
    ),
    Parasitism=.ReadCategoricalTrait(
      raw.data,
      'Parasitism',
      levels=c('Parasitic', 'Not parasitic'),
      reference='Flynn',
      trait='Parasitism',
      kingdoms='Animalia'
    ),
    Modelled_intertegular_distance=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=log10(raw.data$Modelled.ITD..mm.),
      reference='Flynn',
      trait='Modelled intertegular distance',
      kingdom='Animalia',
      unit='log10 (mm)'
    ),
    Body_length=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=log10(raw.data$Body.Length),
      reference='Flynn',
      trait='Body length',
      kingdom='Animalia',
      unit='log10 (mm)'
    ),
    Length_derived_volume=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=3 * log10(raw.data$Body.Length),
      reference='Flynn',
      trait='Length-derived volume',
      kingdom='Animalia',
      unit='3*log10 (mm)'
    ),
    Intertegular_distance=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=log10(raw.data$Intertegular.Distance..mm.),
      reference='Flynn',
      trait='Intertegular distance',
      kingdom='Animalia',
      unit='log10 (mm)'
    ),
    Nests_in_natural_cavities=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=raw.data$Nesting.Strategy.Natural,
      reference='Flynn',
      trait='Nests in natural cavities',
      kingdom='Animalia',
      unit='Logical'
    ),
    Nesting_stratum=.ReadCategoricalTrait(
      raw.data,
      reference='Flynn',
      column='Nesting.Strategy.Stratum',
      trait='Nesting stratum',
      kingdom='Animalia',
      unit='Category',
      levels=c('Burrowing in soil', 'Exposed surfaces',
                'Natural cavities', 'Cavities within plants')
    ),
    Obligately_solitary=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=raw.data$Obligately.solitary,
      reference='Flynn',
      trait='Obligately solitary',
      kingdom='Animalia',
      unit='logical'
    )
  ))
}

WoodDiptera <- function(path, ...) {
  raw.data <- read.csv(path, ...)
  raw.data$Taxon <- paste(raw.data$Genus, raw.data$Species)
  raw.data$Adult.body.length <- MeanOfRangeOrValue(raw.data$Adult.body.length)
  raw.data$Larvae.body.length <- MeanOfRangeOrValue(raw.data$Larvae.body.length)
  raw.data$Wing.size <- MeanOfRangeOrValue(raw.data$Wing.size)

  # Pollination as a logical
  raw.data$Pollination <- tolower(raw.data$Pollination)
  stopifnot(all(raw.data$Pollination %in% c('', 'yes','no')))
  pollination <- 'yes' == raw.data$Pollination
  pollination['' == raw.data$Pollination] <- NA
  raw.data$Pollination <- pollination

  return (list(
    Adult_length=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=log10(raw.data$Adult.body.length),
      reference='Wood',
      trait='Adult length',
      kingdom='Animalia',
      unit='log10 (mm)'
    ),
    Adult_length_derived_volume=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=3 * log10(raw.data$Adult.body.length),
      reference='Wood',
      trait='Adult length-derived volume',
      kingdom='Animalia',
      unit='3*log10 (mm)'
    ),
    Larval_length=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=log10(raw.data$Larvae.body.length),
      reference='Wood',
      trait='Larval length',
      kingdom='Animalia',
      unit='log10 (mm)'
    ),
    Larval_length_derived_volume=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=3 * log10(raw.data$Larvae.body.length),
      reference='Wood',
      trait='Larval length-derived volume',
      kingdom='Animalia',
      unit='3*log10 (mm)'
    ),
    Wing_span=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=log10(raw.data$Wing.size),
      reference='Wood',
      trait='Wing span',
      kingdom='Animalia',
      unit='log10 (mm)'
    ),
    Is_pollinator=TraitDataset(
      taxa=raw.data$Taxon,
      ranks='Species',
      values=raw.data$Pollination,
      reference='Wood',
      trait='Is_pollinator',
      kingdom='Animalia',
      unit='Logical'
    ),
    Adult_niche_breath=.ReadCategoricalTrait(
      raw.data,
      reference='Wood',
      column='Adult.niche.breath',
      trait='Adult niche breadth',
      kingdom='Animalia',
      unit='Category',
      levels=c('1', '2', '3', '4', '5', '>5')
    ),
    Flight_ability=.ReadCategoricalTrait(
      raw.data,
      reference='Wood',
      column='Flight.ability',
      trait='Flight ability',
      kingdom='Animalia',
      unit='Category',
      levels=c('Poor', 'Good')
    ),
    Larval_diet_breath=.ReadCategoricalTrait(
      raw.data,
      reference='Wood',
      column='Larvae.diet.breath',
      trait='Larval diet breath',
      kingdom='Animalia',
      unit='Category'
    ),
    Larval_habitat=.ReadCategoricalTrait(
      raw.data,
      reference='Wood',
      column='Larvae.habitat',
      trait='Larval habitat',
      kingdom='Animalia',
      unit='Category'
    ),
    Larval_niche_breadth=.ReadCategoricalTrait(
      raw.data,
      reference='Wood',
      column='Larvae.niche.breath',
      trait='Larval niche breadth',
      kingdom='Animalia',
      unit='Category'
    ),
    Larval_trophic_level=.ReadCategoricalTrait(
      raw.data,
      reference='Wood',
      column='Larvae.trophic.level',
      trait='Larval trophic level',
      kingdom='Animalia',
      unit='Category'
    ),
    Lifespan=.ReadCategoricalTrait(
      raw.data,
      reference='Wood',
      column='Lifespan',
      trait='Lifespan',
      kingdom='Animalia',
      unit='Category',
      levels=c('Short', 'Medium', 'Long')
    ),
    Reproductive_capacity=.ReadCategoricalTrait(
      raw.data,
      reference='Wood',
      column='Reproductive.capacity',
      trait='Reproductive capacity',
      kingdom='Animalia',
      unit='Category',
      levels=c('Low', 'Medium', 'High')
    )
  ))
}

IUCNHabitat <- function(path, ...){
  raw.data <- read.csv(path, ...)
  
  all.habitats <- .IUCNAllHabitats()
  
  stopifnot(all(raw.data$Habitat %in% all.habitats))
  
  natural.habitats <- .IUCNNaturalHabitats()
  
  spp.natural.major <- droplevels(raw.data[which(raw.data$Habitat %in% natural.habitats),])
  spp.natural.major <- droplevels(spp.natural.major[which(spp.natural.major$Major.importance=="Yes"),])
  
  ret <- data.frame(binomial=unique(raw.data$binomial),nat_hab_spec=FALSE)
  ret$nat_hab_spec[(ret$binomial %in% spp.natural.major$binomial)] <- TRUE
  
  ret$nat_hab_spec <- factor(ret$nat_hab_spec)
  
  return(TraitDataset(reference = 'IUCN Habitat Classification',
                      trait = 'Natural habitat specialization',
                      kingdoms = 'Animalia',
                      unit = 'Category',
                      taxa = ret$binomial,
                      ranks = 'Species',
                      values = ret$nat_hab_spec))
  
  
}


IUCNForest <- function(path, ...){
  raw.data <- read.csv(path, ...)
  
  all.habitats <- .IUCNAllHabitats()
  
  stopifnot(all(raw.data$Habitat %in% all.habitats))
  
  natural.forest.habitats <- .IUCNNaturalForestHabitats()
  
  spp.natural.major <- droplevels(raw.data[which(raw.data$Habitat %in% natural.forest.habitats),])
  spp.natural.major <- droplevels(spp.natural.major[which(spp.natural.major$Major.importance=="Yes"),])
  
  ret <- data.frame(binomial=unique(raw.data$binomial),nat_hab_spec=FALSE)
  ret$nat_hab_spec[(ret$binomial %in% spp.natural.major$binomial)] <- TRUE
  
  ret$nat_hab_spec <- factor(ret$nat_hab_spec)
  
  return(TraitDataset(reference = 'IUCN Habitat Classification',
                      trait = 'Natural forest habitat specialization',
                      kingdoms = 'Animalia',
                      unit = 'Category',
                      taxa = ret$binomial,
                      ranks = 'Species',
                      values = ret$nat_hab_spec))
  
  
}


IUCNHabitat2 <- function(path, ...){
  raw.data <- read.csv(path, ...)
  
  all.habitats <- .IUCNAllHabitats()
  
  stopifnot(all(raw.data$Habitat %in% all.habitats))
  
  natural.habitats <- .IUCNNaturalHabitats()
  
  artificial.habitats <- .IUCNArtificialHabitats()
  
  spp.with.data <- droplevels(raw.data[which((raw.data$Suitability=="Suitable") | (raw.data$Suitability=="Marginal")),])
  
  spp.artificial.suitable <- droplevels(raw.data[which(raw.data$Habitat %in% artificial.habitats),])
  spp.artificial.suitable <- droplevels(spp.artificial.suitable[which(spp.artificial.suitable$Suitability=="Suitable"),])
  
  ret <- data.frame(binomial=unique(raw.data$binomial),nat_hab_spec=TRUE)
  ret$nat_hab_spec[(ret$binomial %in% spp.artificial.suitable$binomial)] <- FALSE
  ret$nat_hab_spec[!(ret$binomial %in% spp.with.data$binomial)] <- NA
  
  ret$nat_hab_spec <- factor(ret$nat_hab_spec)
  
  return(TraitDataset(reference = 'IUCN Habitat Classification',
                      trait = 'Natural habitat specialization',
                      kingdoms = 'Animalia',
                      unit = 'Category',
                      taxa = ret$binomial,
                      ranks = 'Species',
                      values = ret$nat_hab_spec))
  
  
}

IUCNHabitatBreadth <- function(path, ...){
  raw.data <- read.csv(path, ...)
  
  return(TraitDataset(reference = 'IUCN Habitat Classification',
                      trait = 'Habitat breadth',
                      kingdoms = 'Animalia',
                      unit = 'Habitat classes',
                      taxa = raw.data$Best_guess_binomial,
                      ranks = 'Species',
                      values = raw.data$Habitat_breadth_IUCN))
}

.IUCNAllHabitats <- function(){
  return(
    c(
      'Wetlands (inland) - Permanent Rivers/Streams/Creeks (includes waterfalls)',
      'Artificial/Terrestrial - Arable Land',
      'Artificial/Terrestrial - Plantations',
      'Artificial/Terrestrial - Rural Gardens',
      'Artificial/Terrestrial - Subtropical/Tropical Heavily Degraded Former Forest',
      'Forest - Subtropical/Tropical Moist Lowland',
      'Forest - Subtropical/Tropical Moist Montane',
      'Unknown',
      'Wetlands (inland) - Bogs, Marshes, Swamps, Fens, Peatlands',
      'Savanna - Moist',
      'Artificial/Terrestrial - Urban Areas',
      'Artificial/Terrestrial - Pastureland',
      'Wetlands (inland) - Permanent Freshwater Marshes/Pools (under 8ha)',
      'Wetlands (inland) - Seasonal/Intermittent Freshwater Marshes/Pools (under 8ha)',
      'Forest - Subtropical/Tropical Dry',
      'Artificial/Aquatic - Irrigated Land (includes irrigation channels)',
      'Artificial/Aquatic - Seasonally Flooded Agricultural Land',
      'Artificial/Aquatic - Canals and Drainage Channels, Ditches',
      'Savanna - Dry',
      'Wetlands (inland) - Seasonal/Intermittent/Irregular Rivers/Streams/Creeks',
      'Forest - Temperate',
      'Grassland - Subtropical/Tropical Seasonally Wet/Flooded',
      'Wetlands (inland) - Permanent Freshwater Lakes (over 8ha)',
      'Grassland - Subtropical/Tropical High Altitude',
      'Wetlands (inland) - Freshwater Springs and Oases',
      'Forest - Subtropical/Tropical Swamp',
      'Wetlands (inland) - Shrub Dominated Wetlands',
      'Shrubland - Subtropical/Tropical Moist',
      'Grassland - Subtropical/Tropical Dry',
      'Shrubland - Subtropical/Tropical Dry',
      'Artificial/Aquatic - Water Storage Areas (over 8ha)',
      'Artificial/Aquatic - Ponds (below 8ha)',
      'Artificial/Aquatic - Wastewater Treatment Areas',
      'Wetlands (inland) - Seasonal/Intermittent Freshwater Lakes (over 8ha)',
      'Grassland - Temperate',
      'Artificial/Aquatic - Excavations (open)',
      'Forest - Boreal',
      'Shrubland - Temperate',
      'Wetlands (inland) - Permanent Saline, Brackish or Alkaline Lakes',
      'Desert - Temperate',
      'Desert - Cold',
      'Rocky areas (eg. inland cliffs, mountain peaks)',
      'Shrubland - Mediterranean-type Shrubby Vegetation',
      'Wetlands (inland) - Karst and Other Subterranean Hydrological Systems (inland)',
      'Caves and Subterranean Habitats (non-aquatic) - Caves',
      'Caves and Subterranean Habitats (non-aquatic) - Other Subterranean Habitats',
      'Artificial/Aquatic - Aquaculture Ponds',
      'Introduced vegetation',
      'Forest - Subarctic',
      'Wetlands (inland) - Permanent Inland Deltas',
      'Forest - Subtropical/Tropical Mangrove Vegetation Above High Tide Level',
      'Artificial/Aquatic - Karst and Other Subterranean Hydrological Systems (human-made)',
      'Grassland - Tundra',
      'Shrubland - Subarctic',
      'Grassland - Subarctic',
      'Wetlands (inland) - Tundra Wetlands (incl. pools and temporary waters from snowmelt)',
      'Shrubland - Subtropical/Tropical High Altitude',
      'Forest - Subantarctic',
      'Wetlands (inland) - Alpine Wetlands (includes temporary waters from snowmelt)',
      'Desert - Hot',
      'Grassland - Subantarctic',
      'Wetlands (inland) - Geothermal Wetlands',
      'Marine Coastal/Supratidal - Coastal Brackish/Saline Lagoons/Marine Lakes',
      'Other',
      'Marine Coastal/Supratidal - Coastal Freshwater Lakes',
      'Wetlands (inland) - Permanent Saline, Brackish or Alkaline Marshes/Pools',
      'Marine Intertidal - Salt Marshes (Emergent Grasses)',
      'Marine Neritic - Estuaries',
      'Wetlands (inland) - Seasonal/Intermittent Saline, Brackish or Alkaline Lakes and Flats',
      'Shrubland - Subantarctic',
      'Shrubland - Boreal',
      'Artificial/Aquatic - Salt Exploitation Sites',
      'Marine Intertidal - Shingle and/or Pebble Shoreline and/or Beaches',
      'Marine Coastal/Supratidal - Coastal Sand Dunes',
      'Wetlands (inland) - Seasonal/Intermittent Saline, Brackish or Alkaline Marshes/Pools',
      'Marine Oceanic - Epipelagic (0-200m)',
      'Marine Neritic - Subtidal Rock and Rocky Reefs',
      'Marine Neritic - Subtidal Loose Rock/pebble/gravel',
      'Marine Neritic - Subtidal Sandy',
      'Marine Neritic - Subtidal Sandy-Mud',
      'Marine Neritic - Subtidal Muddy',
      'Marine Neritic - Macroalgal/Kelp',
      'Marine Intertidal - Rocky Shoreline',
      'Marine Intertidal - Sandy Shoreline and/or Beaches, Sand Bars, Spits, Etc',
      'Marine Coastal/Supratidal - Sea Cliffs and Rocky Offshore Islands',
      'Marine Coastal/supratidal - Coastal Caves/Karst',
      'Marine Oceanic - Mesopelagic (200-1000m)',
      'Marine Neritic - Pelagic',
      'Marine Oceanic - Bathypelagic (1000-4000m)',
      'Marine Intertidal - Mud Flats and Salt Flats',
      'Marine Intertidal - Mangrove Submerged Roots',
      'Outer Reef Channel',
      'Foreslope (Outer Reef Slope)',
      'Lagoon',
      'Inter-Reef Soft Substrate',
      'Inter-Reef Rubble Substrate',
      'Marine Neritic - Seagrass (Submerged)',
      '',
      'Marine Intertidal - Tidepools',
      'Marine Neritic - Coral Reef',
      'Back Slope',
      'Artificial/Marine - Marine Anthropogenic Structures',
      'Artificial/Marine - Mariculture Cages'
    )
  )
}

.IUCNNaturalHabitats <- function(){
  
  return(
    c(
      'Wetlands (inland) - Permanent Rivers/Streams/Creeks (includes waterfalls)',
      'Forest - Subtropical/Tropical Moist Lowland',
      'Forest - Subtropical/Tropical Moist Montane',
      'Wetlands (inland) - Bogs, Marshes, Swamps, Fens, Peatlands',
      'Savanna - Moist',
      'Wetlands (inland) - Permanent Freshwater Marshes/Pools (under 8ha)',
      'Wetlands (inland) - Seasonal/Intermittent Freshwater Marshes/Pools (under 8ha)',
      'Forest - Subtropical/Tropical Dry',
      'Savanna - Dry',
      'Wetlands (inland) - Seasonal/Intermittent/Irregular Rivers/Streams/Creeks',
      'Forest - Temperate',
      'Grassland - Subtropical/Tropical Seasonally Wet/Flooded',
      'Wetlands (inland) - Permanent Freshwater Lakes (over 8ha)',
      'Grassland - Subtropical/Tropical High Altitude',
      'Wetlands (inland) - Freshwater Springs and Oases',
      'Forest - Subtropical/Tropical Swamp',
      'Wetlands (inland) - Shrub Dominated Wetlands',
      'Shrubland - Subtropical/Tropical Moist',
      'Grassland - Subtropical/Tropical Dry',
      'Shrubland - Subtropical/Tropical Dry',
      'Wetlands (inland) - Seasonal/Intermittent Freshwater Lakes (over 8ha)',
      'Grassland - Temperate',
      'Forest - Boreal',
      'Shrubland - Temperate',
      'Wetlands (inland) - Permanent Saline, Brackish or Alkaline Lakes',
      'Desert - Temperate',
      'Desert - Cold',
      'Rocky areas (eg. inland cliffs, mountain peaks)',
      'Shrubland - Mediterranean-type Shrubby Vegetation',
      'Wetlands (inland) - Karst and Other Subterranean Hydrological Systems (inland)',
      'Caves and Subterranean Habitats (non-aquatic) - Caves',
      'Caves and Subterranean Habitats (non-aquatic) - Other Subterranean Habitats',
      'Forest - Subarctic',
      'Wetlands (inland) - Permanent Inland Deltas',
      'Forest - Subtropical/Tropical Mangrove Vegetation Above High Tide Level',
      'Grassland - Tundra',
      'Shrubland - Subarctic',
      'Grassland - Subarctic',
      'Wetlands (inland) - Tundra Wetlands (incl. pools and temporary waters from snowmelt)',
      'Shrubland - Subtropical/Tropical High Altitude',
      'Forest - Subantarctic',
      'Wetlands (inland) - Alpine Wetlands (includes temporary waters from snowmelt)',
      'Desert - Hot',
      'Grassland - Subantarctic',
      'Wetlands (inland) - Geothermal Wetlands',
      'Marine Coastal/Supratidal - Coastal Brackish/Saline Lagoons/Marine Lakes',
      'Other',
      'Marine Coastal/Supratidal - Coastal Freshwater Lakes',
      'Wetlands (inland) - Permanent Saline, Brackish or Alkaline Marshes/Pools',
      'Marine Intertidal - Salt Marshes (Emergent Grasses)',
      'Marine Neritic - Estuaries',
      'Wetlands (inland) - Seasonal/Intermittent Saline, Brackish or Alkaline Lakes and Flats',
      'Shrubland - Subantarctic',
      'Shrubland - Boreal',
      'Marine Intertidal - Shingle and/or Pebble Shoreline and/or Beaches',
      'Marine Coastal/Supratidal - Coastal Sand Dunes',
      'Wetlands (inland) - Seasonal/Intermittent Saline, Brackish or Alkaline Marshes/Pools',
      'Marine Oceanic - Epipelagic (0-200m)',
      'Marine Neritic - Subtidal Rock and Rocky Reefs',
      'Marine Neritic - Subtidal Loose Rock/pebble/gravel',
      'Marine Neritic - Subtidal Sandy',
      'Marine Neritic - Subtidal Sandy-Mud',
      'Marine Neritic - Subtidal Muddy',
      'Marine Neritic - Macroalgal/Kelp',
      'Marine Intertidal - Rocky Shoreline',
      'Marine Intertidal - Sandy Shoreline and/or Beaches, Sand Bars, Spits, Etc',
      'Marine Coastal/Supratidal - Sea Cliffs and Rocky Offshore Islands',
      'Marine Coastal/supratidal - Coastal Caves/Karst',
      'Marine Oceanic - Mesopelagic (200-1000m)',
      'Marine Neritic - Pelagic',
      'Marine Oceanic - Bathypelagic (1000-4000m)',
      'Marine Intertidal - Mud Flats and Salt Flats',
      'Marine Intertidal - Mangrove Submerged Roots',
      'Outer Reef Channel',
      'Foreslope (Outer Reef Slope)',
      'Lagoon',
      'Inter-Reef Soft Substrate',
      'Inter-Reef Rubble Substrate',
      'Marine Neritic - Seagrass (Submerged)',
      'Marine Intertidal - Tidepools',
      'Marine Neritic - Coral Reef'
    )
  )
}

.IUCNNaturalForestHabitats <- function(){
  
  return(
    c(
      'Forest - Subtropical/Tropical Moist Lowland',
      'Forest - Subtropical/Tropical Moist Montane',
      'Forest - Subtropical/Tropical Dry',
      'Forest - Temperate',
      'Forest - Subtropical/Tropical Swamp',
      'Forest - Boreal',
      'Forest - Subarctic',
      'Forest - Subtropical/Tropical Mangrove Vegetation Above High Tide Level',
      'Forest - Subantarctic'
    )
  )
}


.IUCNArtificialHabitats <- function(){
  
  return(
    c(
      'Artificial/Terrestrial - Arable Land',
      'Artificial/Terrestrial - Plantations',
      'Artificial/Terrestrial - Rural Gardens',
      'Artificial/Terrestrial - Subtropical/Tropical Heavily Degraded Former Forest',
      'Artificial/Terrestrial - Urban Areas',
      'Artificial/Terrestrial - Pastureland',
      'Artificial/Aquatic - Irrigated Land (includes irrigation channels)',
      'Artificial/Aquatic - Seasonally Flooded Agricultural Land',
      'Artificial/Aquatic - Canals and Drainage Channels, Ditches',
      'Artificial/Aquatic - Water Storage Areas (over 8ha)',
      'Artificial/Aquatic - Ponds (below 8ha)',
      'Artificial/Aquatic - Wastewater Treatment Areas',
      'Artificial/Aquatic - Excavations (open)',
      'Artificial/Aquatic - Aquaculture Ponds',
      'Introduced vegetation',
      'Artificial/Aquatic - Karst and Other Subterranean Hydrological Systems (human-made)',
      'Artificial/Aquatic - Salt Exploitation Sites',
      'Artificial/Marine - Marine Anthropogenic Structures',
      'Artificial/Marine - Mariculture Cages'
    )
  )
}

GMPDHostSpecies <- function(path, ...) {
  raw.data <- read.csv(path, ...)
  raw.data$Binomial <- paste(Capitalize(paste(raw.data$Host.genus)),
                             raw.data$Host.species)
  
  parasites <- tapply(X = raw.data$Parasite,INDEX = factor(raw.data$Binomial),
                               function(x) paste(x,collapse=';'))
  
  raw.data <- raw.data[!duplicated(raw.data$Binomial),]
  
  raw.data$parasites <- paste(parasites[match(raw.data$Binomial,names(parasites))])
  
  raw.data$IsHost <- 1
  
  raw.data.fungi <- raw.data[(raw.data$Host_Type=="Fungi"),]
  raw.data.animal <- raw.data[(raw.data$Host_Type %in% c(
    "Amphibian","Arthropod","Aves","Carnivore","Cnidaria",
    "Domestic","Fish","Helminth","Mammal","Mollusca","Primate",
    "Reptile","Rodent","Segmented Worm","Ungulate")),]
  
  if(nrow(raw.data.fungi)==0){
    return(list(AnimalHosts=TraitDataset(reference = "GMPD",
                                         trait = "Is a host",
                                         unit = "binary",
                                         taxa = raw.data.animal$Binomial,
                                         ranks = 'Species',
                                         kingdoms = 'Animalia',
                                         values = raw.data.animal$IsHost)))
  } else {
    return(list(AnimalHosts=TraitDataset(reference = "GMPD",
                                         trait = "Is a host",
                                         unit = "binary",
                                         taxa = raw.data.animal$Binomial,
                                         ranks = 'Species',
                                         kingdoms = 'Animalia',
                                         values = raw.data.animal$IsHost),
                FungiHosts=TraitDataset(reference = "GMPD",
                                        trait = "Is a host",
                                        unit = "binary",
                                        taxa = raw.data.fungi$Binomial,
                                        ranks = 'Species',
                                        kingdoms = 'Fungi',
                                        values = raw.data.fungi$IsHost)))
  }
  
  
}

Reptiles <- function(path, sep=',', ...) {
  raw.data <- read.csv(path,sep=sep,...)
  # Create a dataset for each column
  # A mapping from raw.data column name to a more descriptive name
  diet.cols <- c('Invertebrate'='Diet_invert',
                 'Vertebrate_gen'='Diet_vert_gen',
                 'Vert_Mam_Bir'='Diet_vert_endo',
                 'Vert_Rep_Amp'='Diet_vert_ecto',
                 'Vert_Fis'='Diet_vert_fish',
                 'Vert_Scav'='Diet_scavenging',
                 'Fruit'='Diet_fruit',
                 'Nectar'='Diet_nectar',
                 'Seed'='Diet_seeds',
                 'Plant'='Diet_plant_other',
                 'Funghi'='Diet_fungi')
  
  diet.bad <- 0==rowSums(raw.data[,names(diet.cols)])
  
  diets <- lapply(names(diet.cols), function(col) {
    print(col)
    TraitDataset(reference='Osborne-Tonner Reptiles',
                 trait=diet.cols[col],
                 kingdom='Animalia',
                 unit='Percent',
                 taxa=raw.data$Species[which(!diet.bad)],
                 ranks='Species',
                 values=raw.data[which(!diet.bad),col])
  })
  names(diets) <- unname(diet.cols)
  
  # A diet category derived from the ten diet columns
  diet <- rep('', nrow(raw.data))
  plant.cols <- c('Fruit', 'Nectar', 'Seed', 'Plant','Funghi')
  carnivore.cols <- c('Vertebrate_gen', 'Vert_Mam_Bir','Vert_Rep_Amp','Vert_Fis','Vert_Scav')
  animal.cols <- c('Vertebrate_gen', 'Invertebrate','Vert_Mam_Bir','Vert_Rep_Amp','Vert_Fis','Vert_Scav')
  all.cols<-c('Fruit', 'Nectar', 'Seed', 'Plant','Funghi', 'Vertebrate_gen', 'Invertebrate','Vert_Mam_Bir','Vert_Rep_Amp','Vert_Fis','Vert_Scav')
  
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Fruit==1] <- 'Frugivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Nectar==1] <- 'Nectarivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Seed==1] <- 'Granivore'
  diet[1<=rowSums(raw.data[,plant.cols]) & rowSums(raw.data[,animal.cols])==0] <- 'Herbivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Plant==1]<-'Herbivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Invertebrate==1] <- 'Invertivore'
  diet[rowSums(raw.data[,carnivore.cols])>=1 & rowSums(raw.data[,plant.cols])==0] <- 'Carnivore'
  
  # Anything not fitting one of the above definitins is an omnivore
  diet['' == diet] <- 'Omnivore'
  
  stopifnot(all('' != diet))
  diet <- EnsureFactorLevels(diet, levels=c(
    'Frugivore', 'Nectarivore', 'Granivore', 'Herbivore', 'Plant',
    'Omnivore', 'Invertivore', 'Carnivore', 'Animal'))
  
  
  trophic_level <- paste(diet)
  trophic_level[(trophic_level=="Frugivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Nectarivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Granivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Herbivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Omnivore")] <- "Omnivore"
  trophic_level[(trophic_level=="Invertivore")] <- "Carnivore"
  trophic_level[(trophic_level=="Carnivore")] <- "Carnivore"
  trophic_level <- factor(trophic_level)
  
  
  specialisation <- rep('', nrow(raw.data)) 
  specialisation<- rowSums(raw.data[,all.cols]>0)
  specialisation[0==specialisation]<-NA
  
  return (c(diets,
            list(Diet=TraitDataset(reference='Osborne-Tonner Reptiles',
                                   trait='Diet',
                                   kingdom='Animalia',
                                   unit='Category',
                                   taxa=raw.data$Species,
                                   ranks='Species',
                                   values=diet),
                 Trophic_Level=TraitDataset(reference='Osborne-Tonner Reptiles',
                                            trait='Trophic Level',
                                            kingdom='Animalia',
                                            unit='Category',
                                            taxa=raw.data$Species,
                                            ranks='Species',
                                            values=trophic_level),
                 Specialisation=TraitDataset(reference='Osborne-Tonner Reptiles',
                                             trait='Specialisation',
                                             kingdom='Animalia',
                                             unit='Number of categories',
                                             taxa=raw.data$Species,
                                             ranks='Species',
                                             values=specialisation))))
}



Amphibians <- function(path, sep=',', ...) {
  raw.data <- read.csv(path,sep=sep,...)
  
  # Create a dataset for each column
  # A mapping from raw.data column name to a more descriptive name
  diet.cols <- c('Invertebrate'='Diet_invert',
                 'Vertebrate'='Diet_vert_gen',
                 'Fruit'='Diet_fruit',
                 'Nectar'='Diet_nectar',
                 'Seed'='Diet_seeds',
                 'Plant'='Diet_plant_other')
  
  
  diet.bad <- 0==rowSums(raw.data[,names(diet.cols)])
  
  diets <- lapply(names(diet.cols), function(col) {
    print(col)
    TraitDataset(reference='Osborne-Tonner Amphibians',
                 trait=diet.cols[col],
                 kingdom='Animalia',
                 unit='Percent',
                 taxa=raw.data$Scientific[which(!diet.bad)],
                 ranks='Species',
                 values=raw.data[which(!diet.bad),col])
  })
  names(diets) <- unname(diet.cols)
  
  # A diet category derived from the ten diet columns
  diet <- rep('', nrow(raw.data))
  plant.cols <- c('Fruit', 'Nectar', 'Seed', 'Plant')
  animal.cols <- c('Vertebrate', 'Invertebrate')
  all.cols<-c('Fruit', 'Nectar', 'Seed', 'Plant', 'Vertebrate', 'Invertebrate')
  
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Fruit==1] <- 'Frugivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Nectar==1] <- 'Nectarivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Seed==1] <- 'Granivore'
  diet[1<rowSums(raw.data[,plant.cols]) & rowSums(raw.data[,animal.cols])==0] <- 'Herbivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Plant==1]<-'Herbivore'
  diet[rowSums(raw.data[,all.cols])==1 & raw.data$Invertebrate==1] <- 'Invertivore'
  diet[rowSums(raw.data[,animal.cols])>=1 & rowSums(raw.data[,plant.cols])==0] <- 'Carnivore'
  
  
  # Anything not fitting one of the above definitins is an omnivore
  diet['' == diet] <- 'Omnivore'
  
  stopifnot(all('' != diet))
  diet <- EnsureFactorLevels(diet, levels=c(
    'Frugivore', 'Nectarivore', 'Granivore', 'Herbivore', 'Plant',
    'Omnivore', 'Invertivore', 'Carnivore', 'Animal'))
  
  
  trophic_level <- paste(diet)
  trophic_level[(trophic_level=="Frugivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Nectarivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Granivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Herbivore")] <- "Herbivore"
  trophic_level[(trophic_level=="Omnivore")] <- "Omnivore"
  trophic_level[(trophic_level=="Invertivore")] <- "Carnivore"
  trophic_level[(trophic_level=="Carnivore")] <- "Carnivore"
  trophic_level[(trophic_level=="NA")]<-NA
  trophic_level <- factor(trophic_level)
  
  
  specialisation <- rep('', nrow(raw.data)) 
  specialisation<- rowSums(raw.data[,all.cols]>0)
  specialisation[0==specialisation]<-NA
  
  
  return (c(diets,
            list(Diet=TraitDataset(reference='Osborne-Tonner Amphibians',
                                   trait='Diet',
                                   kingdom='Animalia',
                                   unit='Category',
                                   taxa=raw.data$Scientific,
                                   ranks='Species',
                                   values=diet),
                 Trophic_Level=TraitDataset(reference='Osborne-Tonner Amphibians',
                                            trait='Trophic Level',
                                            kingdom='Animalia',
                                            unit='Category',
                                            taxa=raw.data$Scientific,
                                            ranks='Species',
                                            values=trophic_level),
                 Specialisation=TraitDataset(reference='Osborne-Tonner Amphibians',
                                             trait='Specialisation',
                                             kingdom='Animalia',
                                             unit='Number of categories',
                                             taxa=raw.data$Scientific,
                                             ranks='Species',
                                             values=specialisation))))
}


