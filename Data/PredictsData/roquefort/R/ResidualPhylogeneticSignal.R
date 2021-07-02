ResidualPhylogeneticSignal<-function(model,data,commonLevel){
  
  model.data<-cbind(model@frame,data[match(model@frame$Taxon_name_entered,data$Taxon_name_entered),
                                     c('Kingdom','Phylum','Class','Order','Family','Genus','Best_guess_binomial')])
  
  model.taxa<-unique(model.data[,c('Kingdom','Phylum','Class','Order','Family','Genus','Best_guess_binomial')])
  model.taxa<-model.taxa[(model.taxa$Best_guess_binomial!=""),]
  
  if(commonLevel=="Kingdom"){
    phylog<-as.phylo(~Phylum/Class/Order/Family/Genus/Best_guess_binomial,data=model.taxa)
  } else if (commonLevel=="Phylum"){
    phylog<-as.phylo(~Class/Order/Family/Genus/Best_guess_binomial,data=model.taxa)
  } else if (commonLevel=="Class"){
    phylog<-as.phylo(~Order/Family/Genus/Best_guess_binomial,data=model.taxa)
  } else if (commonLevel=="Order"){
    phylog<-as.phylo(~Family/Genus/Best_guess_binomial,data=model.taxa)
  } else if (commonLevel=="Family"){
    phylog<-as.phylo(~Genus/Best_guess_binomial,data=model.taxa)
  }
  
  phylogFinal<-compute.brlen(phylog,method="Grafen")
  phylogFinal<-multi2di(phylogFinal,random=TRUE)
  
  model.data$residual<-residuals(model)

  model.data<-model.data[(model.data$Best_guess_binomial!=""),]

  res<-aggregate(model.data$residual~model.data$Best_guess_binomial,FUN=mean)
  nam<-res[,1]
  res<-res[,2]
  names(res)<-nam
  
  
  lambdaFit<-fitContinuous(phy=phylogFinal,res,model="lambda",control=list(hessian=TRUE))
  
  return(list(lambda=lambdaFit$opt$lambda,
              lower=lambdaFit$opt$CI['lb','lambda'],
              upper=lambdaFit$opt$CI['ub','lambda']))
  
}