# TODO: make model_select return the family of the model to avoid having to pass it here
DisaggregateMixedModel<-function(model,divFactor,fitFamily,data=NULL,fixedFactors=
                          character(0),fixedTerms=list(),
                        fixedInteractions=character(0),
                        fitInteractions=FALSE,randomStruct=NULL){
  cat("Assembling model call\n")
  
  if ((!("data" %in% names(model))) & (is.null(data))){
    stop("If model is not a roquefort model, then data must be supplied")
  }
  
  if ("model" %in% names(model)){
    responseVar<-strsplit(as.character(model$model@call[[2]]),'~')[[2]][1]
  } else {
    responseVar<-strsplit(as.character(model@call[[2]]),'~')[[2]][1]
  }
  
  
  allTerms<-character(0)
  fixedStruct<-""
  for (i in 1:length(fixedFactors)){
    fixedStruct<-paste(fixedStruct,fixedFactors[i],sep="")
    allTerms<-c(allTerms,fixedFactors[i])
    if ((i != length(fixedFactors)) | (length(fixedTerms)>0) | 
          ((length(fixedTerms)==0) & (
            length(fixedInteractions)>0))){
      fixedStruct<-paste(fixedStruct,"+",sep="")
    }
  }
  if (length(fixedTerms)>0){
    for (i in 1:length(fixedTerms)){
      term<-paste("poly(",names(fixedTerms)[i],
                  ",",fixedTerms[i],")",sep="")
      fixedStruct<-paste(fixedStruct,term,sep="")
      allTerms<-c(allTerms,term)
      if ((i != length(fixedTerms)) | (length(fixedInteractions)>0)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }
  
  if (fitInteractions){
    fixedStruct<-paste(fixedStruct,"+")
    
    mainTerms<-allTerms
    
    for (i in 1:(length(mainTerms)-1)){
      for (j in (i+1):length(mainTerms)){
        term<-paste(mainTerms[i],mainTerms[j],sep=":")
        fixedStruct<-paste(fixedStruct,term)
        allTerms<-c(allTerms,term)
      }
    }
    
  }
  
  if (length(fixedInteractions)>0){
    for (i in 1:length(fixedInteractions)){
      fixedStruct<-paste(fixedStruct,fixedInteractions[i],sep="")
      allTerms<-c(allTerms,fixedInteractions[i])
      if (i != length(fixedInteractions)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }
  
  if(!is.null(randomStruct)) randomStruct<-gsub(" ","",randomStruct)
  
  model.call<-paste(responseVar,"~",allTerms[1],sep="")
  if (length(allTerms)>1){
    for (t in 2:length(allTerms)){
      model.call<-paste(model.call,"+",allTerms[t],sep="")
    }
  }
  if(!is.null(randomStruct)){
    model.call<-paste(model.call,"+",randomStruct,sep="")
  } 
  
  cat("Preparing data\n")
  if ("data" %in% names(model)){
    modelData<-model$data
  } else {
    modelData<-data
  }
  
  levels<-data.frame(level=unique(modelData[,divFactor]))
  
  nFolds<-dim(levels)[1]
  levels$fold<-1:nFolds
    
  cat("Running overall model\n")
  if(is.null(randomStruct)){
    model<-try(glm(formula=model.call,data=modelData,family=fitFamily),
               silent=TRUE)
    all.results<-data.frame(matrix(data=NA,nrow=nFolds,ncol=length(coefficients(model))))
    names(all.results)<-names(coefficients(model))
  } else {
    model<-try(glmer(formula=model.call,data=modelData,family=fitFamily),
               silent=TRUE)
    all.results<-data.frame(matrix(data=NA,nrow=nFolds,ncol=length(fixef(model))))
    names(all.results)<-names(fixef(model))
  }
  
  row.names(all.results)<-levels$level
  
  cat("Performing disaggregation\n")
  for (f in 1:nFolds){
    cat(paste("\rPerforming iteration ",f," of ",nFolds,sep=""))
    
    sub.levels<-levels[levels$fold==f,]
    
    sub.data<-modelData[(modelData[,divFactor] %in% sub.levels$level),]
    
    if(is.null(randomStruct)){
      sub.model<-try(glm(formula=model.call,data=sub.data,family=fitFamily),
                     silent=TRUE)
      if(class(sub.model)!="try-error"){
        all.results[f,match(names(coefficients(sub.model)),names(all.results))]<-coefficients(sub.model)
      }
    } else {
      sub.model<-try(glmer(formula=model.call,data=sub.data,family=fitFamily),
                     silent=TRUE)
      if (class(sub.model)!="try-error"){
        all.results[f,match(names(fixef(sub.model)),names(all.results))]<-fixef(sub.model)
      }
      
    }
    
    
    
  }
  
  fixef.quantiles<-apply(all.results,2,quantile,probs=c(0.025,0.975),na.rm=TRUE)
  
  return(list(all.results=all.results,confidence=fixef.quantiles))
}