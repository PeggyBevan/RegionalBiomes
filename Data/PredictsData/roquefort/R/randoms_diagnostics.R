RandomsDiagnostics<-function(model,randomName,fname,wDir=getwd()){
  
  setwd(wDir)
  
  pdf(fname,width=8.5/2.54,height=6.5/2.54)
  
  par(cex.lab=0.7)
  par(cex.axis=0.4)
  par(mgp=c(2,1,0))
  par(mar=c(3,3,1,1))
  
  eval(substitute(ranefs<-ranef(model)$x,list(x=randomName)))
  eval(substitute(ran.names<-names(ranef(model)$x),list(x=randomName)))
       
  for (r in ran.names){
    eval(substitute(qqnorm(ranefs$x,main=x),list(x=r)))
    eval(substitute(qqline(ranefs$x),list(x=r)))
  }
  
  dev.off()
  
}