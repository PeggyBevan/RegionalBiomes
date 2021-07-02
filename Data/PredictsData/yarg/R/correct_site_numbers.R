CorrectSiteNumbers <- function(diversity, sites) {
  # A faster implementation of correct_site_numbers
  # Combines sites that have the same coordinates, sampling dates, habitat and 
  # land use.

  # Ensure that SS is present
  diversity <- .AddSummaryColumns(diversity)
  sites <- .AddSummaryColumns(sites)

  unique.cols <- c('Longitude','Latitude','Predominant_habitat',
                   'Use_intensity', 'Sample_start_earliest', 
                   'Sample_end_latest','Sample_date_resolution')
  # Sites within each study
  sites.within.studies <- split(sites[,unique.cols], sites$SS)
  suspect <- sapply(sites.within.studies, function(rows) nrow(rows)!=nrow(unique(rows)))
  for(study in names(which(suspect))) {
    rows <- sites[sites$SS==study,unique.cols]
    unique.rows <- unique(rows)
    .Log('From', nrow(rows), 'sites to', nrow(unique.rows), 'sites in', study, '\n')
   
    # Paste together those columns that should be unique
    ids <- sapply(1:nrow(rows), function(n) paste(do.call('c', rows[n,unique.cols]), collapse=' '))

    # Split site numbers by id
    duplicated.sites <- split(sites[sites$SS==study,'Site_number'], ids)

    # Consider only sites that are duplicated
    duplicated.sites <- duplicated.sites[sapply(duplicated.sites, function(s) length(s)>1)]

    # Set all measurements for the duplicated sites to the first of the site numbers
    for(n in duplicated.sites) {
      diversity[diversity$SS==study & diversity$Site_number %in% n,'Site_number'] <- n[1]
    }
  }
  return (diversity)
}

correct_site_numbers<-function(dataset){
  dataset$StudyID<-paste(dataset$Source_ID,dataset$Study_number,sep="_")
  
  for (s in unique(dataset$StudyID)){
    temp.data<-subset(dataset,StudyID==s)
    
    real.sites<-unique(paste(temp.data$Longitude,temp.data$Latitude,temp.data$
                               Predominant_habitat,temp.data$Use_intensity,
                             temp.data$Sample_start_earliest,temp.data$
                               Sample_end_latest,temp.data$Sample_date_resolution))
    
    if (length(real.sites)!=length(unique(temp.data$Site_number))){
      real.sites<-data.frame(real.sites)
      real.sites$id<-1:dim(real.sites)[1]
      
      match.sites<-paste(temp.data$Longitude,temp.data$Latitude,temp.data$
                           Predominant_habitat,temp.data$Use_intensity,
                         temp.data$Sample_start_earliest,temp.data$
                           Sample_end_latest,temp.data$Sample_date_resolution)
      dataset$Site_number[which(dataset$StudyID==s)]<-real.sites$id[match(match.sites,real.sites$real.sites)]
      
      
    }
    
  }
  
  return(dataset)
  
}


