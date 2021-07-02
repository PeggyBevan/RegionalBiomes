TraitPhylogeneticSignal<-function(data,commonLevel,trait){
  
  
  taxa<-unique(data[,c('Kingdom','Phylum','Class','Order','Family','Genus','Best_guess_binomial')])
  taxa<-taxa[(taxa$Best_guess_binomial!=""),]
  
  if(commonLevel=="Kingdom"){
    phylog<-as.phylo(~Phylum/Class/Order/Family/Genus/Best_guess_binomial,data=taxa)
  } else if (commonLevel=="Phylum"){
    phylog<-as.phylo(~Class/Order/Family/Genus/Best_guess_binomial,data=taxa)
  } else if (commonLevel=="Class"){
    phylog<-as.phylo(~Order/Family/Genus/Best_guess_binomial,data=taxa)
  } else if (commonLevel=="Order"){
    phylog<-as.phylo(~Family/Genus/Best_guess_binomial,data=taxa)
  } else if (commonLevel=="Family"){
    phylog<-as.phylo(~Genus/Best_guess_binomial,data=taxa)
  }
  
  phylogFinal<-compute.brlen(phylog,method="Grafen")
  phylogFinal<-multi2di(phylogFinal,random=TRUE)
  
  data<-data[(data$Best_guess_binomial!=""),]
  
  traitVals<-aggregate(data[,trait]~data$Best_guess_binomial,FUN=mean)

  nam<-traitVals[,1]
  traitVals<-traitVals[,2]
  names(traitVals)<-nam
  
  
  lambdaFit<-fitContinuous(phy=phylogFinal,traitVals,model="lambda",control=list(hessian=TRUE))
  
  return(list(lambda=lambdaFit$opt$lambda,
              lower=lambdaFit$opt$CI['lb','lambda'],
              upper=lambdaFit$opt$CI['ub','lambda']))
  
}