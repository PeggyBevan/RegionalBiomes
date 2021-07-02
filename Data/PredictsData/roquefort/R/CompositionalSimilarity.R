
CompositionalSimilarity<-function(data,metric){
  
  # Take subsets of the data depending on the selected metric
  # (records of abundance only for abundance-weighted Sorensen,
  # and records recorded as number of individuals for corrected Sorensen)
  if (metric=="SorAbd"){
    data<-data[data$Diversity_metric_type=="Abundance",]
  }
  if (metric=="SorCorr"){
    data<-data[((data$Diversity_metric=="abundance") & 
                  (data$Diversity_metric_unit=="individuals")),]
  }
  
  # Make a column for land use and rename
  data$LandUse<-paste(data$Predominant_habitat)
  data$LandUse[which(data$LandUse=="Primary forest")]<-"Primary Vegetation"
  data$LandUse[which(data$LandUse=="Primary non-forest")]<-"Primary Vegetation"
  data$LandUse[which(data$LandUse=="Secondary vegetation (indeterminate age)")]<-NA
  data$LandUse[which(data$LandUse=="Secondary non-forest")]<-NA
  data$LandUse[which(data$LandUse=="Cannot decide")]<-NA
  data$LandUse<-factor(data$LandUse)
  data$LandUse<-relevel(data$LandUse,ref="Primary Vegetation")
  
  # Subset the data to just the required columns and remove NAs
  data<-subset(data,select=c("SS","SSBS","Measurement","Taxon_name_entered",
                             "LandUse"))
  data<-na.omit(data)
  
  # If the selected metric is corrected Sorensen, only use studies where all
  # measurements are integers
  if (metric=="SorCorr"){
    study.all.int.meas <- 
      tapply(data$Measurement,data$SS, 
             function(m) all(floor(m)==m))
    
    int.meas<-study.all.int.meas[match(data$SS,names(study.all.int.meas))]
    
    data<-data[int.meas,]
    
  }
  
  # Get a list of unique land uses
  all.lu<-unique(paste(data$LandUse))
  
  # Make matrices to store the final results, temporary sums and counts
  all.results<-matrix(nrow=length(all.lu),ncol=length(all.lu))
  sum.matrix<-matrix(0,nrow=length(all.lu),ncol=length(all.lu))
  count.matrix<-matrix(0,nrow=length(all.lu),ncol=length(all.lu))
  
  # Convert the final results matrix to a data frame and name rows and columns
  all.results<-data.frame(all.results)
  names(all.results)<-all.lu
  row.names(all.results)<-all.lu
  
  # Crate integers to count studies and used studies
  all.studies.count<-0
  used.studies.n<-0
  
  # Create character vector to store list of used studies
  used.studies<-character(0)
  
  # Loop over studies
  for (st in unique(data$SS)){
    # Increment the count of studies
    all.studies.count<-all.studies.count+1
    cat(paste("Processing ",st,"\n",sep=""))
    
    # Take the data for the current study
    sub.data<-data[data$SS==st,]
    
    # If the metric is the corrected Sorensen, only use non-zero measurements
    if (metric != "SorCorr"){
      sub.data<-sub.data[sub.data$Measurement>0,]
    }
    
    # Make a site-by-site matrix
    sites.matrix<-matrix(nrow=length(unique(sub.data$SSBS)),ncol=length(unique(sub.data$SSBS)))
    
    # Loop over all possible pairs of sites
    i1<-1
    for (s1 in unique(sub.data$SSBS)){
      i2<-1
      for (s2 in unique(sub.data$SSBS)){
        # Calculate similarity using the desired metric
        if (metric=="Sor"){
          
          u<-length(union(sub.data$Taxon_name_entered[sub.data$SSBS==s1],
                          sub.data$Taxon_name_entered[sub.data$SSBS==s2]))
          i<-length(intersect(sub.data$Taxon_name_entered[sub.data$SSBS==s1],
                              sub.data$Taxon_name_entered[sub.data$SSBS==s2]))
          
          sor<-(2*i)/((2*i)+(u-i))
          
          
        } else if (metric=="SorAbd"){
          s1.sum<-sum(sub.data$Measurement[(sub.data$SSBS==s1)])
          s2.sum<-sum(sub.data$Measurement[(sub.data$SSBS==s2)])
          
          u<-sum(sub.data$Measurement[(sub.data$SSBS==s1) & 
                                        (sub.data$Taxon_name_entered %in% 
                                           intersect(sub.data$Taxon_name_entered[
                                             sub.data$SSBS==s1],sub.data$
                                               Taxon_name_entered[sub.data$
                                                                    SSBS==s2]))]/s1.sum)
          v<-sum(sub.data$Measurement[(sub.data$SSBS==s2) & 
                                        (sub.data$Taxon_name_entered %in% 
                                           intersect(sub.data$Taxon_name_entered[
                                             sub.data$SSBS==s1],sub.data$
                                               Taxon_name_entered[sub.data$
                                                                    SSBS==s2]))]/s2.sum)
          sor<-(2*u*v)/(u+v)
        } else if (metric=="SorCorr"){         
          
          n<-sum(sub.data$Measurement[sub.data$SSBS==s1])
          m<-sum(sub.data$Measurement[sub.data$SSBS==s2])
          
          if ((n > 0) & (m > 0)){
            xi<-sub.data$Measurement[sub.data$SSBS==s1][(match(union(
              sub.data$Taxon_name_entered[sub.data$SSBS==s1],
              sub.data$Taxon_name_entered[sub.data$SSBS==s2]),
              sub.data$Taxon_name_entered[sub.data$SSBS==s1]))]
            yi<-sub.data$Measurement[sub.data$SSBS==s2][(match(union(
              sub.data$Taxon_name_entered[sub.data$SSBS==s1],
              sub.data$Taxon_name_entered[sub.data$SSBS==s2]),
              sub.data$Taxon_name_entered[sub.data$SSBS==s2]))]
            
            xi[is.na(xi)]<-0
            yi[is.na(yi)]<-0
            
            f1.<-length(which((xi==1) & (yi>0)))
            f2.<-max(1,length(which((xi==2) & (yi>0))))
            f.1<-length(which((xi>0) & (yi==1)))
            f.2<-max(1,length(which((xi>0) & (yi==2))))
            
            p1<-sum(xi[yi>0]/n)
            p2<-((m-1)/m)*(f.1/(2*f.2))
            p3<-sum(xi[yi==1]/n)
            
            u<-min(1,p1+p2*p3)
            
            q1<-sum(yi[xi>0]/m)
            q2<-((n-1)/n)*(f1./(2*f2.))
            q3<-sum(yi[xi==1]/m)
            
            v<-min(1,q1+q2*q3)
            
            
            if ((u>0) & (v>0)){
              sor<-(2*u*v)/(u+v)
            } else {
              sor<-0
            }
            
          } else {
            sor <- 0
          }
          
          
        } else {
          stop("Error: specfied dissimilarity metric is not supported")
        }
        
        if (s1!=s2) sites.matrix[i1,i2]<-sor
        i2<-i2+1
      }
      i1<-i1+1
    }
    
    # Get the land uses for all of the sites in this study
    site.lu<-paste(sub.data$LandUse)[match(unique(sub.data$SSBS),sub.data$SSBS)]
    
    # Make a matrix to hold average compositional similarity for each pair of land uses
    lu.matrix<-matrix(nrow=length(unique(site.lu)),ncol=length(unique(site.lu)))
    
    # Calculate average compositional similarity for each pair of land uses
    for (lu1 in 1:length(unique(site.lu))){
      for (lu2 in 1:length(unique(site.lu))){
        lu.matrix[lu1,lu2]<-mean(sites.matrix[which(site.lu==unique(site.lu)[lu1]),
                                              which(site.lu==unique(site.lu)[lu2])],na.rm=T)
        
      }
    }
    
    # Only consider a study if it sampled a primary vegetation site
    # and if there was more than one site in at least some land uses
    pind<-which(unique(site.lu)=="Primary Vegetation")
    if (("Primary Vegetation" %in% site.lu) & (FALSE %in% is.na(lu.matrix))){
      if (!is.na(lu.matrix[pind,pind])){
        # Increment the counter of used studies and add the study name to the list
        used.studies.n<-used.studies.n+1
        used.studies<-c(used.studies,st)
        # Rescale the matrix of similarities by land use to be 1 for primary
        # vegetation compared with itself
        lu.matrix<-lu.matrix/lu.matrix[which(unique(site.lu)=="Primary Vegetation"),
                                       which(unique(site.lu)=="Primary Vegetation")]
        # Increment the matrix of counts for land use pairs represented in this study
        count.matrix[match(unique(site.lu),all.lu),match(unique(site.lu),all.lu)]<-
          count.matrix[match(unique(site.lu),all.lu),match(unique(site.lu),all.lu)]+
          ifelse(is.na(lu.matrix),0,1)
        # Replace NAs in the matrix of similarities by land use pairs with 0
        # (to avoid converting values in the sum matrix to NA)
        lu.matrix[is.na(lu.matrix)]<-0
        # Add similarities between land uses from this study to the sum matrix
        # for all studies
        sum.matrix[match(unique(site.lu),all.lu),match(unique(site.lu),all.lu)]<-
          sum.matrix[match(unique(site.lu),all.lu),match(unique(site.lu),all.lu)]+
          lu.matrix
        
      }
      
    }
    
  }
  
  # Calculate average similarities for pairs of land uses across all studies
  all.results<-sum.matrix/count.matrix
  
  # Convert the matrix of average values to a data frame
  all.results<-data.frame(all.results)
  
  # Name the rows and columns in the table of average values by land use
  names(all.results)<-all.lu
  row.names(all.results)<-all.lu
  
  all.results[is.na(all.results)]<-NA
  
  cat(paste("Considered ",used.studies.n," of ",all.studies.count,"\n",sep=""))
  
  return(list(cd=all.results,n=count.matrix,studies=used.studies))
}
