
map_study_randoms<-function(model,site.data,overlayMap=NULL,slopeLevel="(Intercept)"){
  
  
  eval(substitute(randoms<-ranef(model)$StudyID$X,list(X=slopeLevel)))
  
  source.names<-row.names(ranef(model)$StudyID)
  
  longs<-aggregate(Longitude~StudyID,data=site.data,FUN=mean)
  lats<-aggregate(Latitude~StudyID,data=site.data,FUN=mean)
  
  sources<-merge(longs,lats,all.x=T,all.y=T)
  sources$raneff<-randoms[match(sources$StudyID,source.names)]
  
  sources<-na.omit(sources)
  
  cols<-brewer.pal(11,"RdYlBu")
  
  par(mar=c(0,0,0,0))
  
  if(!is.null(overlayMap)){
    plot(overlayMap)
  }
  
  points(sources$Long,sources$Lat,pch=20,col=as.character(cut(sources$raneff,
                                                              breaks=11,
                                                              labels=cols)))
  
  text(0,90,slopeLevel)
  
}
