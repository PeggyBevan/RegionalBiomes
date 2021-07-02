predictGLMER<-function(model,newdata=model@frame,se.fit=FALSE,seMultiplier=1.96,logLink="n"){
  
  
  mm<-model.matrix(terms(model),newdata)
  pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)),mm))
  y<-mm %*% fixef(model)
  if(se.fit){
    yplus<-y+sqrt(pvar1)*seMultiplier
    yminus<-y-sqrt(pvar1)*seMultiplier
  }
  
  if (logLink=="e"){
    y<-(exp(y))
    if(se.fit){
      yplus<-(exp(yplus))
      yminus<-(exp(yminus))
    }
  }  else if (logLink=="10"){
    y<-(10^(y))
    if(se.fit){
      yplus<-(10^(yplus))
      yminus<-(10^(yminus))
    }
  } else if (logLink=="n"){
    
  } else {
    stop("Error: the specified log link is not supported")
  }
      
  if(se.fit){
    r<-data.frame(y=y,yplus=yplus,yminus=yminus)
  } else {
    r<-y[,1]
  }
  
  return(r)

}
