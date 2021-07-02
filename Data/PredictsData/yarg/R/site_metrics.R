.SiteTotalAbundance <- function(diversity, site.abundance) {
  .Log("Computing total abundance\n")
  ta <- rep(NA, length(site.abundance))
  ta[site.abundance] <- tapply(diversity$Measurement[diversity$Is_abundance], 
                               droplevels(diversity$SSS[diversity$Is_abundance]), 
                               sum)
  return (ta)
}

.SiteSpeciesRichness <- function(diversity, site.abundance, site.occurrence, 
                                 site.species.richness, sites.are.unique) {

  interesting <- diversity$Is_abundance | diversity$Is_occurrence
  if(sites.are.unique) {
    .Log('Computing species richness\n')
    # The correct way - taxa are unique within sites
    # The number of non-zero abundances and/or counts
    sr <- rep(NA, length(site.abundance))
    sr[site.abundance | site.occurrence] <- 
      tapply(diversity$Measurement[interesting], 
             droplevels(diversity$SSS[interesting]), 
             function(m) sum(m>0))
  } else {
    # The broken way - taxa are not unique within sites.
    # The number of unique taxon names entered with a non-zero measurement
    # This is a temporary hack to accomodate combining sites in a way that 
    # breaks the guarantee that sites are unique within studies.
    .Log('HACK: computing species richness using taxon name entered\n')
    v <- by(diversity[interesting ,c('Measurement','Taxon_name_entered')], 
             droplevels(diversity$SSS[interesting]), 
             function(rows) length(unique(rows$Taxon_name_entered[rows$Measurement>0])))
    sr <- rep(NA, length(site.abundance))
    sr[site.abundance | site.occurrence] <- unlist(as.list(v))
  }

  sr[site.species.richness] <- 
    tapply(diversity$Measurement[diversity$Is_species_richness], 
           droplevels(diversity$SSS[diversity$Is_species_richness]), 
           sum)
  
  return (sr)
}

.SiteWeightedSpeciesRichness <- function(diversity, site.abundance, site.occurrence, 
                                 site.species.richness, richWeights, sites.are.unique) {
  
  interesting <- diversity$Is_abundance | diversity$Is_occurrence
  
  stopifnot(sites.are.unique)
  
  diversity$occur <- ifelse(diversity$Measurement>0.0,1,0)
  diversity$weightedOccur <- diversity$occur * richWeights
  
  .Log('Computing weighted species richness\n')
  # The correct way - taxa are unique within sites
  # The number of non-zero abundances and/or counts
  sr <- rep(NA, length(site.abundance))
  sr[site.abundance | site.occurrence] <- 
    tapply(diversity$weightedOccur[interesting], 
           droplevels(diversity$SSS[interesting]), 
           function(m) sum(m,na.rm=TRUE))
  
  return (sr)
}


.SiteSimpsonsDiversity <- function(diversity, site.abundance) {
  .Log("Computing Simpson's diversity\n")
  sd <- rep(NA, length(site.abundance))
  sd[site.abundance] <- tapply(diversity$Measurement[diversity$Is_abundance], 
                               droplevels(diversity$SSS[diversity$Is_abundance]), 
                               function(m) {
                                 if(any(m>0)) 1/(sum((m/sum(m))^2))
                                 else         NA
                               })
  return (sd)
}

.SiteChao <- function(diversity, site.abundance) {
  .Log("Computing Chao\n")
  # Compute Chao only for studies where Diversity_metric_is_suitable_for_Chao 
  # is TRUE and all measurements within the study are integers.
  
  # http://stackoverflow.com/questions/3476782/how-to-check-if-the-number-is-integer
  # TODO all.equal
  study.chao.suitable <- 
    tapply(diversity$Measurement[diversity$Diversity_metric_is_suitable_for_Chao],
           diversity$SS[diversity$Diversity_metric_is_suitable_for_Chao], 
           function(m) all(floor(m)==m))
  
  # Studies for which Diversity_metric_is_suitable_for_Chao is FALSE will have 
  # a value of NA in study.chao.suitable
  study.chao.suitable[is.na(study.chao.suitable)] <- FALSE
  site.chao.suitable <- study.chao.suitable[tapply(diversity$SS, diversity$SSS, unique)]
  
  diversity$Is_Chao_suitable <- site.chao.suitable[diversity$SSS]
  
  # Only abundance metrics should be marked as suitable for Chao 
  stopifnot(!any(site.chao.suitable & !site.abundance))
  chao <- rep(NA, length(site.abundance))
  chao[site.chao.suitable] <- 
    tapply(diversity$Measurement[diversity$Is_Chao_suitable], 
           droplevels(diversity$SSS[diversity$Is_Chao_suitable]), 
           function(m) sum(m>0) + (((sum(m==1) * (sum(m==1)-1)) / (2*(sum(m==2)+1)))))
  return (chao)
}

.SiteRarefiedRichness<-function(diversity, site.abundance){
  .Log("Computing Rarefied Species Richness\n")
  # Compute rarefied richness only for studies where Diversity_metric_is_suitable_for_Chao 
  # is TRUE and all measurements within the study are integers.
  
  # http://stackoverflow.com/questions/3476782/how-to-check-if-the-number-is-integer
  # TODO all.equal
  study.rsr.suitable <- 
    tapply(diversity$Measurement[diversity$Diversity_metric_is_suitable_for_Chao],
           diversity$SS[diversity$Diversity_metric_is_suitable_for_Chao], 
           function(m) all(floor(m)==m))
  
  # Studies for which Diversity_metric_is_suitable_for_Chao is FALSE will have 
  # a value of NA in study.rsr.suitable
  study.rsr.suitable[is.na(study.rsr.suitable)] <- FALSE
  site.rsr.suitable <- study.rsr.suitable[tapply(diversity$SS, diversity$SSS, unique)]
  
  diversity$Is_RSR_suitable <- site.rsr.suitable[diversity$SSS]
  
  # Only abundance metrics should be marked as suitable for rarefied richness 
  stopifnot(!any(site.rsr.suitable & !site.abundance))
  
  rsrich <- rep(NA, length(site.abundance))
  
  if (any(diversity$Is_RSR_suitable)){
    siteTotalAbund<-aggregate(Measurement~SSS,data=droplevels(
      diversity[diversity$Is_RSR_suitable,]),FUN=sum)
    siteTotalAbund$SS<-diversity$SS[match(siteTotalAbund$SSS,diversity$SSS)]
    studyMinSiteAbund<-aggregate(Measurement~SS,data=siteTotalAbund,FUN=
                                   function(x) min(x[x>0]))
    diversity$studyMinSiteAbund<-studyMinSiteAbund$Measurement[match(
      diversity$SS,studyMinSiteAbund$SS)]
    
    sample.n <- split(diversity$studyMinSiteAbund[diversity$Is_RSR_suitable], 
                      droplevels(diversity$SSS[diversity$Is_RSR_suitable]))
    values <- split(diversity$Measurement[diversity$Is_RSR_suitable],
                    droplevels(diversity$SSS[diversity$Is_RSR_suitable]))
    
    rsr<-matrix(nrow=length(which(site.rsr.suitable)),ncol=1000)
    for (i in 1:1000){
      rsr[,i]<-mapply(function(vals,samp){
        sp<-as.character(1:length(vals))
        ind<-rep(sp,vals)
        if (length(ind)==0){
          sr<-0
        } else {
          rare.sp<-sample(ind,samp[1])
          sr<-length(table(rare.sp))
        }
        return(sr)
      },values,sample.n)
    }
    
    rsrich[site.rsr.suitable]<-apply(rsr,1,mean)
    
  }
  
  
  
  
  return(rsrich)
  
}

.AppendSiteCWMTrait <- function(site.metrics, diversity, traits,
                                ignore.abund=FALSE,removeOutliers=NULL) {
  # Computes community-weighted mean trait values for each trait in traits.
  # Returns the data.frame site.metrics with a column per trait appended.
  # Also records proportion with trait data
  if(length(traits)>0) {
    # A list with a vector of values per site
    prop.found <- function(values, weights){
      tot.wt <- sum(weights)
      tot.found <-sum(weights[!is.na(values)])
      result <- tot.found/tot.wt
      return(result)
    }
    if(ignore.abund){
      diversity$occur <- ifelse(diversity$Measurement>0,1,0)
      weights <- split(diversity$occur,diversity$SSS)
    } else {
      weights <- split(diversity$Measurement, diversity$SSS)
    }
    for(trait in traits) {
      .Log("Computing community-weighted mean [", trait, "]\n", sep='')
      stopifnot(trait %in% colnames(diversity))
      stopifnot(class(diversity[,trait]) %in% c('numeric','integer','logical'))

      values <- split(diversity[,trait], diversity$SSS)
      
      if(!is.null(removeOutliers)){
        stopifnot(class(removeOutliers) %in% c('numeric','integer'))
        lowerProb <- (removeOutliers/100)/2
        upperProb <- 1 - (removeOutliers/100)/2
        
        values <- lapply(X = values,FUN = function(vals){
          
          lowerQuantile <- quantile(x = vals,probs=lowerProb)
          upperQuantile <- quantile(x = vals,probs=upperProb)
          
          vals[vals<lowerQuantile] <- NA
          vals[vals>upperQuantile] <- NA
          
          return(vals)
          
        })
        
      }
      
      cwm <- mapply(weighted.mean, values, weights, na.rm=TRUE)
      prop.trait.data <- mapply(prop.found, values, weights)
      site.metrics[,paste('CWM_', trait, sep='')] <- cwm
      site.metrics[,paste('PropWithData_', trait, sep='')] <- prop.trait.data
    }
  }
  return (site.metrics)
}

.AppendSiteTraitIPR<- function(site.metrics, diversity, traits,
                               centralPercentile) {
  # Computes trait inter-percentile range for each trait in traits. 
  # Returns the data.frame site.metrics with a column per trait appended.
  if(length(traits>0)) {
    stopifnot(0<centralPercentile && centralPercentile<100)
    if (centralPercentile < 1) .Log("WARNING: percentile should be % value, not fraction\n")
    centPercent <- centralPercentile/100
    outerPercent <- ((1-centPercent)/2)

    # A list with a vector of measurements per site
    measurements <- split(diversity$Measurement, diversity$SSS)

    for(trait in traits) {
      .Log("Computing [", trait, "] inter-percentile range\n", sep='')
      stopifnot(trait %in% colnames(diversity))
      trait.values <- split(diversity[,trait], diversity$SSS)
      ipr <- mapply(measurements, trait.values, FUN=function(m, tv) {
        # Do not consider measurements for which the trait is NA
        m <- m[!is.na(tv)]
        tv <- tv[!is.na(tv)]

        # Order measurements and trait values by trait values
        o <- order(tv)
        m <- m[o]
        tv <- tv[o]

        sum.N<-sum(m)
        cum.N<-cumsum(m)

        N.lower <- sum.N*outerPercent
        N.upper <- sum.N*(centPercent+outerPercent)

        trait.lower <- tv[which(cum.N>=N.lower)[1]]
        trait.upper <- tv[which(cum.N>=N.upper)[1]]

        return(trait.upper-trait.lower)
      })
      site.metrics[,paste('IPR_', trait, sep='')] <- ipr
    }
  }
  return (site.metrics)
}

.AppendSiteTraitPercentiles<- function(site.metrics, diversity, traits,
                               centralPercentile) {
  # Computes trait inter-percentile range for each trait in traits. 
  # Returns the data.frame site.metrics with a column per trait appended.
  if(length(traits>0)) {
    stopifnot(0<centralPercentile && centralPercentile<100)
    if (centralPercentile < 1) .Log("WARNING: percentile should be % value, not fraction\n")
    centPercent <- centralPercentile/100
    outerPercent <- ((1-centPercent)/2)
    
    # A list with a vector of measurements per site
    measurements <- split(diversity$Measurement, diversity$SSS)
    
    for(trait in traits) {
      .Log("Computing [", trait, "] percentiles\n", sep='')
      stopifnot(trait %in% colnames(diversity))
      trait.values <- split(diversity[,trait], diversity$SSS)
      lower <- mapply(measurements, trait.values, FUN=function(m, tv) {
        # Do not consider measurements for which the trait is NA
        m <- m[!is.na(tv)]
        tv <- tv[!is.na(tv)]
        
        # Order measurements and trait values by trait values
        o <- order(tv)
        m <- m[o]
        tv <- tv[o]
        
        sum.N<-sum(m)
        cum.N<-cumsum(m)
        
        N.lower <- sum.N*outerPercent
        
        trait.lower <- tv[which(cum.N>=N.lower)[1]]
        
        return(trait.lower)
      })
      upper <- mapply(measurements, trait.values, FUN=function(m, tv) {
        # Do not consider measurements for which the trait is NA
        m <- m[!is.na(tv)]
        tv <- tv[!is.na(tv)]
        
        # Order measurements and trait values by trait values
        o <- order(tv)
        m <- m[o]
        tv <- tv[o]
        
        sum.N<-sum(m)
        cum.N<-cumsum(m)
        
        N.upper <- sum.N*(centPercent+outerPercent)
        
        trait.upper <- tv[which(cum.N>=N.upper)[1]]
        
        return(trait.upper)
      })
      site.metrics[,paste('LowerPercentile_', trait, sep='')] <- lower
      site.metrics[,paste('UpperPercentile_', trait, sep='')] <- upper
    }
  }
  return (site.metrics)
}

SiteMetrics <- function(diversity, extra.cols=NULL, traits=NULL,
                        CWMtraits.ignore.abund=FALSE,
                        centralPercentile=95, sites.are.unique=TRUE,
                        srEstimators = c("Chao","Rare"),richWeights=NULL,
                        removeOutliers=NULL) {
  # Takes a data.frame of diversity measurements and returns a data.frame 
  # of site-level matrics.
  # diversity should have columns
  #   Source_ID
  #   Study_number
  #   Site_number
  #   SS (Source_ID Study_number)
  #   SSS (Source_ID Study_number Site_number)
  #   Diversity_metric_is_valid
  #   Diversity_metric_type
  #   Diversity_metric_is_effort_sensitive
  #   Diversity_metric_is_suitable_for_Chao
  #   Sampling_effort
  #   Sampling_method_is_valid
  #   Measurement
  #   All columns in 'traits'
  
  # The Diverity_* and Sampling_* columns are each guaranteed to be unique 
  # within studies. Because sites are nested within studies, these columns are 
  # also guaranteed to be unique within sites.
  .Log(paste('Computing site metrics for', nrow(diversity), 'measurements\n'))

  stopifnot(all(diversity$Diversity_metric_is_valid))
  stopifnot(all(!is.na(diversity$Sampling_effort)))

  site.cols <- c('Source_ID','Study_number','Study_name','Site_number',
                 'Site_name','Block','Predominant_habitat','Use_intensity',
                 'Longitude','Latitude','Sample_start_earliest',
                 'Sample_end_latest','Diversity_metric_type')
  
  if(!is.null(extra.cols)) site.cols <- union(site.cols, extra.cols)
  
  # Drop names of columns that are not present
  site.cols <- intersect(site.cols, colnames(diversity))

  .Log('The data contain', nlevels(diversity$Source_ID), 'sources,', 
       nlevels(diversity$SS), 'studies and',
       nlevels(diversity$SSS), 'sites\n')
  
  # Each site must have the same values for all columns in site.cols
  if(FALSE) {
    .Log('Checking for duplicated within-site values\n')
    # TODO This code is too slow. Faster to check column by column?
    if(!all(1==by(diversity[,site.cols], diversity$SSS, function(rows) nrow(unique(rows))))) {
      stop('Some values in site.cols are not unique within sites')
    }
  } else {
    .Log('TODO fix within-site uniqueness check\n')
  }
  
  if(any('Diversity index'==diversity$Diversity_metric_type)) {
    stop('Diversity indices are not yet supported')
  }

  bad <- setdiff(diversity$Diversity_metric_type, c('Abundance','Occurrence','Species richness'))
  if(length(bad)>0) {
    stop('Unrecognied diversity metrics ', paste(bad, collapse=','))
  }
  
  .Log('Computing site-level values\n')
  
  # Measurement-level values
  diversity$Is_abundance <- 'Abundance'==diversity$Diversity_metric_type
  diversity$Is_occurrence <- 'Occurrence'==diversity$Diversity_metric_type
  diversity$Is_species_richness <- 'Species richness'==diversity$Diversity_metric_type
  
  # Site-level values
  site.abundance <- tapply(diversity$Is_abundance, diversity$SSS, unique)
  site.occurrence <- tapply(diversity$Is_occurrence, diversity$SSS, unique)
  site.species.richness <- tapply(diversity$Is_species_richness, diversity$SSS, unique)
  
  total.abundance <- .SiteTotalAbundance(diversity, site.abundance)
  species.richness <- .SiteSpeciesRichness(diversity, site.abundance, 
                                           site.occurrence, 
                                           site.species.richness, 
                                           sites.are.unique)
  simpsons <- .SiteSimpsonsDiversity(diversity, site.abundance)
  
  if(!is.null(richWeights)){
    weighted.species.richness <- .SiteWeightedSpeciesRichness(diversity, site.abundance, 
                                                              site.occurrence, 
                                                              site.species.richness,
                                                              richWeights,
                                                              sites.are.unique)
  }
  
  if ("Chao" %in% srEstimators){
    chao <- .SiteChao(diversity, site.abundance)
  }
  
  if ("Rare" %in% srEstimators){
    rsrich<-.SiteRarefiedRichness(diversity,site.abundance)
  }
  
  # Extract site-level values
  .Log('Assembling site-level values\n')
  res <- cbind(diversity[!duplicated(diversity$SSS),c(site.cols,'SS','SSS')], 
               Total_abundance=total.abundance, 
               Species_richness=species.richness, 
               Simpson_diversity=simpsons, 
               ChaoR=if ("Chao" %in% srEstimators) chao else NA,
               Richness_rarefied=if ("Rare" %in% srEstimators) rsrich else NA,
               Weighted_richness=if (!is.null(richWeights)) weighted.species.richness else NA)
  rownames(res) <- NULL
  
  # Add community-weighted mean trait values
  res <- .AppendSiteCWMTrait(res, diversity, traits,ignore.abund=CWMtraits.ignore.abund,
                             removeOutliers = removeOutliers)
  
  # Add trait IPR
  res <- .AppendSiteTraitIPR(res, diversity, traits, centralPercentile)
  
  # Add trait percentiles
  res <- .AppendSiteTraitPercentiles(res, diversity, traits, centralPercentile)
  
  # The order of results will not be sensible, e.g. site 11 will come before 
  # site 2. Reorder sensibly.
  res <- res[order(res$Source_ID, res$Study_number, res$Site_number),]
  
  # Sanity checks
  stopifnot(all.equal(as.character(res$SS), 
                      paste(res$Source_ID, res$Study_number)))
  stopifnot(all.equal(as.character(res$SSS), 
                      paste(res$Source_ID, res$Study_number, res$Site_number)))
  
  return (res)
}

