
SpatialAutocorrelationTest<-function(model,all.data,siteModel=TRUE){
  
  if ("data" %in% names(model)){
    model.data<-model$data
    model.data$res<-residuals(model$model)
  } else {
    model.data<-model@frame
    model.data$res<-residuals(model)
  }
  
  
  model.data$Longitude<-all.data$Longitude[match(model.data$SSBS,all.data$SSBS)]
  model.data$Latitude<-all.data$Latitude[match(model.data$SSBS,all.data$SSBS)]
  
    
  studies<-character()
  failed<-character()
  moran.i<-numeric()
  moran.p<-numeric()
  
  i=1
  for (ss in unique(model.data$SS)){
    cat(paste("\rProcessing study ",i," of ",length(unique(model.data$SS)),sep=""))
    data.sub<-droplevels(model.data[model.data$SS==ss,])

    if(!siteModel){
      resids<-tapply(data.sub$res,data.sub$SSBS,mean)
      long<-tapply(data.sub$Longitude,data.sub$SSBS,mean)
      lat<-tapply(data.sub$Latitude,data.sub$SSBS,mean)
      data.sub<-data.frame(res=resids,Longitude=long,Latitude=lat)
    }
    
    ds.nb<-try(dnearneigh(cbind(data.sub$Longitude,data.sub$Latitude),
                          d1=0.00000001,d2=10),silent=TRUE)
    ds.listw<-try(nb2listw(ds.nb),silent=TRUE)
    mt<-tryCatch(moran.test(data.sub$res,ds.listw),silent=TRUE,error=function(e) e, 
                 warning=function(w) w)
    
    if(class(mt)[1]=="htest"){
      if ((!is.na(mt$statistic))){
        studies<-c(studies,ss)
        moran.i<-c(moran.i,mt$statistic)
        moran.p<-c(moran.p,mt$p.value)
      } else {
        failed<-c(failed,ss)
      }
      
    } else {
      failed<-c(failed,ss)
    }
    
    i<-i+1
  }
  
  return(list(studies=studies,I=moran.i,P=moran.p,failed=failed))
  
  
  
}

SpatialAutocorrelationTestCombined<-function(model,all.data){
  
  if ("data" %in% names(model)){
    model.data<-model$data
  } else {
    model.data<-model@frame
  }
  model.data$res<-residuals(model$model)
  model.data$Longitude<-all.data$Longitude[match(model.data$SSBS,all.data$SSBS)]
  model.data$Latitude<-all.data$Latitude[match(model.data$SSBS,all.data$SSBS)]
  
  ds.nb<-try(dnearneigh(cbind(model.data$Longitude,model.data$Latitude),
                        d1=0.00000001,d2=10),silent=TRUE)
  ds.listw<-try(nb2listw(ds.nb),silent=TRUE)
  mt<-try(moran.test(model.data$res,ds.listw),silent=TRUE)
  
  if(class(mt)=="htest"){
    if ((!is.na(mt$statistic))){
      moran.i<-mt$statistic
      moran.p<-mt$p.value
    } else {
      moran.i<-NA
      moran.p<-NA
    }
    
  } else {
    moran.i<-NA
    moran.p<-NA
  }
    
  return(list(I=moran.i,P=moran.p))
  
}
