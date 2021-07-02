
relate_randoms<-function(model,randomName,data,colName,outDir=NULL){
  
  eval(substitute(data$ranef<-ranef(model)$x$'(Intercept)'[match(data$x,row.names(ranef(model)$x))],list(x=randomName)))
  
  eval(substitute(data2<-aggregate(data$ranef~data$x+data$y,FUN=max),list(x=randomName,y=colName)))
  
  names(data2)<-gsub("data[$]","",names(data2))
  
  if(!is.null(outDir)){
    png(paste(outDir,"/",names(model@frame)[1]," ",randomName," randoms by ",colName,".png",
              sep=""),width=20/2.54,height=14/2.54,units="in",res=1200)
  }
  
  par(las=2)
  par(mar=c(14,2,1,1))
  
  eval(substitute(boxplot(data2$ranef~data2$x),list(x=colName)))
  
  if(!is.null(outDir)){
    dev.off()
  }
  
  eval(substitute(m1<-lm(data2$ranef~factor(data2$x)),list(x=colName)))
  
  return(list(adj.r.squared=summary(m1)$adj.r.squared,
              r.squared=summary(m1)$r.squared,F=anova(m1)$'F value'[1],
         P=anova(m1)$'Pr(>F)'[1]))

  
}

