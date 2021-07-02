
CompositionalSimilarityData<-function(data,metric,nIters,weights=NULL,
                                      luClassification="lu.only",adjust=FALSE,
                                      envVars=NULL,
                                      other.diffs=NULL){
  
  if(!is.null(envVars)){
    cat("Calculating environmental distance matrix\n")
    sites<-diversity[,c('SSBS',envVars)]
    sites<-unique(sites)
    site_env<-sites[,envVars]
    row.names(site_env)<-sites$SSBS
    env.dists<-as.matrix(gowdis(site_env))
  }
  
  
  if ((metric=="SorWeight") & (is.null(weights))){
    stop("If metric is weighted Sorensen, weights column must be specified")
  }
  
  for (other.diff in other.diffs){
    data.cols<-names(data)
    stopifnot(other.diff %in% data.cols)
  }
  
  # Take subsets of the data depending on the selected metric
  # (records of abundance only for abundance-weighted Sorensen,
  # and records recorded as number of individuals for corrected Sorensen)
  if ((metric=="SorAbd") | (metric=="JaccAbdAsymm") | (metric=="BC")){
    data<-data[data$Diversity_metric_type=="Abundance",]
  }
  if (metric=="SorCorr"){
    data<-data[((data$Diversity_metric=="abundance") & 
                  (data$Diversity_metric_unit=="individuals")),]
  }
  if (metric=="SorWeight"){
    data<-data[(!is.na(data[,weights])),]
  }
  
  if (luClassification=="lu.only"){
    # Make a column for land use and rename
    data$LandUse<-paste(data$Predominant_habitat)
    data$LandUse[which(data$LandUse=="Primary forest")]<-"Primary Vegetation"
    data$LandUse[which(data$LandUse=="Primary non-forest")]<-"Primary Vegetation"
    data$LandUse[which(data$LandUse=="Secondary vegetation (indeterminate age)")]<-NA
    data$LandUse[which(data$LandUse=="Secondary non-forest")]<-NA
    data$LandUse[which(data$LandUse=="Cannot decide")]<-NA
    data$LandUse<-factor(data$LandUse)
    data$LandUse<-relevel(data$LandUse,ref="Primary Vegetation")
  } else if (luClassification=="lu.pa"){
    data$LandUse<-paste(data$Predominant_habitat)
    data$LandUse[which(data$LandUse=="Primary forest")]<-"Primary Vegetation"
    data$LandUse[which(data$LandUse=="Primary non-forest")]<-"Primary Vegetation"
    data$LandUse[which(data$LandUse=="Secondary vegetation (indeterminate age)")]<-"Secondary Vegetation"
    data$LandUse[which(data$LandUse=="Secondary non-forest")]<-"Secondary Vegetation"
    data$LandUse[which(data$LandUse=="Young secondary vegetation")]<-"Secondary Vegetation"
    data$LandUse[which(data$LandUse=="Intermediate secondary vegetation")]<-"Secondary Vegetation"
    data$LandUse[which(data$LandUse=="Mature secondary vegetation")]<-"Secondary Vegetation"
    data$LandUse[which(data$LandUse=="Cannot decide")]<-NA
    
    data$LandUse<-paste(data$LandUse,data$Within_PA)
    
    data$LandUse[grep("Primary",data$LandUse)]<-"Primary Vegetation"
    data$LandUse[grep("NA",data$LandUse)]<-NA
    
    data$LandUse<-factor(data$LandUse)
    data$LandUse<-relevel(data$LandUse,ref="Primary Vegetation")
    
  } else if (luClassification=="lu.allsec"){
    
    data$LandUse<-paste(data$Predominant_habitat)
    data$LandUse[which(data$LandUse=="Primary forest")]<-"Primary Vegetation"
    data$LandUse[which(data$LandUse=="Primary non-forest")]<-"Primary Vegetation"
    data$LandUse[which(data$LandUse=="Secondary vegetation (indeterminate age)")]<-"Secondary Vegetation"
    data$LandUse[which(data$LandUse=="Secondary non-forest")]<-"Secondary Vegetation"
    data$LandUse[which(data$LandUse=="Young secondary vegetation")]<-"Secondary Vegetation"
    data$LandUse[which(data$LandUse=="Intermediate secondary vegetation")]<-"Secondary Vegetation"
    data$LandUse[which(data$LandUse=="Mature secondary vegetation")]<-"Secondary Vegetation"
    data$LandUse[which(data$LandUse=="Cannot decide")]<-NA
    
    data$LandUse<-factor(data$LandUse)
    data$LandUse<-relevel(data$LandUse,ref="Primary Vegetation")
    
  } else {
    stop("Land-use classification not recognized")
  }
  
  
  
  # Create list to store counts of sites in different land uses
  count.lus<-as.list(rep(0,length(na.omit(unique(data$LandUse)))))
  names(count.lus)<-na.omit(unique(data$LandUse))
  
  # Subset the data to just the required columns and remove NAs
  if(metric=="SorWeight"){
    data<-subset(data,select=c("SS","SSBS","Measurement","Taxon_name_entered",
                               "LandUse","Longitude","Latitude",weights,other.diffs,
                               envVars))
  } else {
    data<-subset(data,select=c("SS","SSBS","Measurement","Taxon_name_entered",
                               "LandUse","Longitude","Latitude",other.diffs,
                               envVars))
  }
  
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
  
  # Make a list to store all final compositional similarity data in
  all.data<-list()
  
  # Create a data frame for each random subset of the data to be created
  for (i in 1:nIters){
    all.data[[i]]<-data.frame(sim=numeric(0),
                              lu.comparison=character(0),
                              SS=character(0))
    
    for (j in other.diffs){
      all.data[[i]][,paste(j,"diff",sep=".")]=numeric(0)
    }
    
  }
  
  
  # Loop over studies
  for (st in unique(data$SS)){
    # Increment the count of studies
    all.studies.count<-all.studies.count+1
    cat(paste("Processing ",st,"\n",sep=""))
    
    # Take the data for the current study
    sub.data<-data[data$SS==st,]
    
    # If the metric is the corrected Sorensen, only use non-zero measurements
    if ((metric != "SorCorr") & (metric != "BC")){
      sub.data<-sub.data[sub.data$Measurement>0,]
    }
    
    if(dim(sub.data)[1]>0){
      # Make a site-by-site matrix
      sites.matrix<-matrix(nrow=length(unique(sub.data$SSBS)),
                           ncol=length(unique(sub.data$SSBS)))
      
      
      if(!is.null(other.diffs)){
        sm.other<-list()
        for (j in 1:length(other.diffs)){
          sm.other[[j]]<-matrix(nrow=length(unique(sub.data$SSBS)),
                                ncol=length(unique(sub.data$SSBS)))
        }
      }
      
      
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
            
            if(adjust){
              
              div.sub.sub<-droplevels(sub.data[((sub.data$SSBS==s1) | (sub.data$SSBS==s2)),])
              
              sp.rich<-tapply(div.sub.sub$Taxon_name_entered,div.sub.sub$SSBS,
                              function(x) return(length(unique(x))))
              
              max.spp<-max(sp.rich)
              min.spp<-min(sp.rich)
              
              spp.ratio<-min.spp/max.spp
              
              sorMax<-(2*spp.ratio)/((2*spp.ratio)+(1-spp.ratio))
              
              sor<-sor/sorMax
              
            }
            
          } else if (metric=="Sim"){
            a<-length(union(sub.data$Taxon_name_entered[sub.data$SSBS==s1],
                            sub.data$Taxon_name_entered[sub.data$SSBS==s2]))
            b<-length(which(!(sub.data$Taxon_name_entered[sub.data$SSBS==s1] %in% 
                                sub.data$Taxon_name_entered[sub.data$SSBS==s2])))
            c<-length(which(!(sub.data$Taxon_name_entered[sub.data$SSBS==s2] %in% 
                                sub.data$Taxon_name_entered[sub.data$SSBS==s1])))
            
            sor<-(min(b,c))/(min(b,c)+a)
            
            
          } else if (metric=="SorWeight"){
            a.spp<-intersect(sub.data$Taxon_name_entered[sub.data$SSBS==s1],
                             sub.data$Taxon_name_entered[sub.data$SSBS==s2])
            b.spp<-setdiff(sub.data$Taxon_name_entered[sub.data$SSBS==s1],
                           sub.data$Taxon_name_entered[sub.data$SSBS==s2])
            c.spp<-setdiff(sub.data$Taxon_name_entered[sub.data$SSBS==s2],
                           sub.data$Taxon_name_entered[sub.data$SSBS==s1])
            
            a<-sum(sub.data[,weights][match(a.spp,sub.data$Taxon_name_entered)])
            b<-sum(sub.data[,weights][match(b.spp,sub.data$Taxon_name_entered)])
            c<-sum(sub.data[,weights][match(c.spp,sub.data$Taxon_name_entered)])
            
            sor<-(2*a)/(2*a+b+c)
            
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
          } else if (metric=="BC"){
            
            if (!all(sub.data$Taxon_name_entered[(sub.data$SSBS==s1)]==
                     sub.data$Taxon_name_entered[(sub.data$SSBS==s2)])){
              stop("Taxon names don't match")
            } 
            
            sor<-1-((sum(abs(sub.data$Measurement[(sub.data$SSBS==s1)]-
                               sub.data$Measurement[(sub.data$SSBS==s2)])))/(
                                 sum(sub.data$Measurement[(sub.data$SSBS==s1)])+
                                   sum(sub.data$Measurement[(sub.data$SSBS==s2)])))
            
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
            
          } else if (metric=="JaccAsymm"){
            c <- length(setdiff(sub.data$Taxon_name_entered[sub.data$SSBS==s2],
                       sub.data$Taxon_name_entered[sub.data$SSBS==s1]))
            a <- length(intersect(sub.data$Taxon_name_entered[sub.data$SSBS==s2],
                           sub.data$Taxon_name_entered[sub.data$SSBS==s1]))
            
            sor <- a/(a+c)
            
          } else if (metric=="JaccAbdAsymm"){
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
            
            sor <- (u*v)/u
          } else {
            stop("Error: specfied dissimilarity metric is not supported")
          }
          
          if (s1!=s2) sites.matrix[i1,i2]<-sor
          
          if(!is.null(other.diffs)){
            for (j in 1:length(other.diffs)){
              sm.other[[j]][i1,i2]<-abs(sub.data[which(sub.data$SSBS==s1)[1],other.diffs[j]]-
                                          sub.data[which(sub.data$SSBS==s2)[1],other.diffs[j]])
            }
          }
          
          
          i2<-i2+1
        }
        i1<-i1+1
      }
      
      # Loop over required number of random subsets of the compositional similarity data
      for (i in 1:nIters){
        # Generate a random order for the sites
        randorder<-sample(x = dim(sites.matrix)[1],size = dim(sites.matrix)[1],replace = FALSE)
        
        # Rearrange the sites matrix into this random order
        sites.matrix.rand<-sites.matrix[randorder,randorder]
        
        # Keep only the first and second off-diagonals as independent comparisons 
        sites.matrix.rand[row(sites.matrix)-col(sites.matrix)<1]<-NA
        sites.matrix.rand[row(sites.matrix)-col(sites.matrix)>1]<-NA
        
        # Rearrange the matrices of other specified differences
        if(!is.null(other.diffs)){
          sm.other.rand<-list()
          for (j in 1:length(other.diffs)){
            sm.other.rand[[j]]<-sm.other[[j]][randorder,randorder]
          }
        }
        
        # Get the land uses for all of the sites in this study, reordered into the random order
        site.lu<-paste(sub.data$LandUse)[match(unique(sub.data$SSBS),sub.data$SSBS)][
          randorder]
        
        # Make a matrix to hold average compositional similarity for each pair of land uses
        lu.matrix<-matrix(paste(rep(site.lu,length(site.lu)),rep(site.lu,each=length(site.lu)),
                                sep="-"),nrow=length(site.lu),ncol=length(site.lu))
        
        # Get the environmental distances between all pairs of sites
        if(!is.null(envVars)){
          env.dist.matrix<-env.dists[unique(sub.data$SSBS)[randorder],
                                     unique(sub.data$SSBS)[randorder]]
        }
        
        # Create a data frame of unique sites, reordered into the random order
        site.data<-sub.data[(!duplicated(sub.data$SSBS)),][randorder,]
        # Calculate geographic distances between all pairs of sites
        dist.mat<-distm(x=site.data[,c('Longitude','Latitude')])
        
        # Increment the counts of sites in different land uses (only on 1st iteration)
        if (i==1){
          lu.counts<-summary(factor(site.lu))
          count.lus[names(lu.counts)]<-unlist(count.lus[names(lu.counts)])+lu.counts
        }
        
        # Make a data frame with compositional similarity and distance data for this study
        dframe<-data.frame(sim=sites.matrix.rand[lower.tri(sites.matrix.rand)],
                           lu.comparison=lu.matrix[lower.tri(lu.matrix)],
                           SS=rep(sub.data$SS[1],length(sites.matrix.rand[
                             lower.tri(sites.matrix.rand)])),
                           dist=dist.mat[lower.tri(dist.mat)])
        if(!is.null(envVars)){
          dframe$env_dist=env.dist.matrix[lower.tri(env.dist.matrix)]
        }
        # Add any other requested differences among sites
        if(!is.null(other.diffs)){
          for (j in 1:length(other.diffs)){
            dframe[,paste(other.diffs[j],"diff",sep=".")]<-sm.other.rand[[j]][
              lower.tri(sm.other.rand[[j]])]
          }
          
        }
        
        # Remove NAs
        dframe<-na.omit(dframe)
        
        # Add the data frame to the data for all other studies under the appropriate random iteration
        all.data[[i]]<-merge(all.data[[i]],dframe,all=TRUE)
      }
      
      
    }
    
    
  }
  
  return(list(data=all.data,counts=count.lus))
}
