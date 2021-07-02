OverdispersionTest<-function(model){
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  resid.df <- nrow(model.frame(model))-model.df
  
  resids<-residuals(model,type="pearson")
  resid.dev<-sum(residuals(model,type="pearson")^2)
  
  ratio<-resid.dev/resid.df
  
  p<-pchisq(resid.dev,df=resid.df,lower.tail=FALSE)
  
  cat(paste("\nResidual deviance =",resid.dev))
  cat(paste("\nResidual DF =",resid.df))
  cat(paste("\nRatio =",ratio))
  cat(paste("\nChi-square: P =",p))
  
}