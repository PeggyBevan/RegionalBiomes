compile_best_svl<-function(amphibDataDir){
  
  #Set the working directory
  setwd(amphibDataDir)
  
  # Read in data from Cooper et al. (2008)
  CooperData<-read.csv("Cooper_svl_data.csv")
  
  # Read in data for other species in PREDICTS database, which are not covered
  # by Cooper et al.
  SeniorData<-read.csv("Senior_svl_data.csv")
  
  # Extract best estimates from supplementary dataset
  SeniorData$best.svl<-NA
  
  # Define ordering of best to worst SVL estimates (calculated midrange,calculated
  # mean,given midrange,given mean,given min, given max,female midrange, male 
  # midrange, female mean,male mean, female min, female max, male min, male max)
  
  col.order<-c(14,13,18,16,17,15,12,8,9,5,10,11,6,7)
  
  for (c in col.order){
    SeniorData$best.svl[is.na(SeniorData$best.svl)]<-
      SeniorData[,c][is.na(SeniorData$best.svl)]
  }
  
  # Add column for identifying source of data
  
  CooperData$origin<-"Cooper_et_al_2008"
  SeniorData$origin<-"Senior_svl_data_2013"
  
  # SeniorData currently has data for 96 species not covered by Cooper et al. (57
  # still lacking SVL, 12 of which have a congener in Cooper et al.)
  
  # Create new dataset by adding new data onto that from Cooper et al.
  
  # Subset by necessary columns and get rid of NAs
  best.data<-na.omit(CooperData[c("binomial","svl","origin")]) 
  SeniorData<-na.omit(SeniorData[c("Taxon","best.svl","origin")])
  
  # Match column names
  colnames(SeniorData)<-c("binomial","svl","origin")
  
  # Combine and return new dataset
  best.data<-rbind(best.data,SeniorData)
  return(best.data)
  
}