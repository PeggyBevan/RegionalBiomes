
LanduseSiteDistances<-function(sites,luColumn="Predominant_habitat",method=distHaversine){
  
  all.studies<-unique(sites$SS)
  
  all.lu<-unique(paste(sites[,luColumn]))
  all.results<-array(dim=c(length(all.lu),length(all.lu),length(all.studies)))
  
  i<-1
  for (st in all.studies){
    cat(paste("\rProcessing study (",i," of ",length(all.studies),"): ",st,paste(rep(" ",30),collapse=" "),sep=""))
    
    study.sites<-sites[(sites$SS==st),]
  
    dist.mat<-distm(x=study.sites[,c('Longitude','Latitude')])
    
    diag(dist.mat)<-NA
    
    site.lu<-paste(study.sites[,luColumn])
    
    lu.matrix<-matrix(nrow=length(unique(site.lu)),ncol=length(unique(site.lu)))
    
    for (lu1 in 1:length(unique(site.lu))){
      for (lu2 in 1:length(unique(site.lu))){
        lu.matrix[lu1,lu2]<-mean(dist.mat[which(site.lu==unique(site.lu)[lu1]),
                                              which(site.lu==unique(site.lu)[lu2])],na.rm=T)
        
      }
    }
    
    if (-Inf %in% lu.matrix) stop()
    
    all.results[match(unique(site.lu),all.lu),match(unique(site.lu),all.lu),i]<-lu.matrix
    
    i<-i+1
    
  }
  
  lu.means<-apply(X=all.results,MARGIN=c(1,2),FUN=function(x) mean(log(x+1),na.rm=TRUE))
  lu.ses<-apply(X=all.results,MARGIN=c(1,2),FUN=function(x) sd(log(x+1),na.rm=TRUE)/sqrt(length(x)))
  
  lu.means<-data.frame(lu.means)
  lu.ses<-data.frame(lu.ses)
  names(lu.means)<-all.lu
  row.names(lu.means)<-all.lu
  names(lu.ses)<-all.lu
  row.names(lu.ses)<-all.lu
  
  
  return(list(Landuses=all.lu,LanduseMeans=lu.means,LanduseStdErrs=lu.ses,FullArray=all.results))
  
}
