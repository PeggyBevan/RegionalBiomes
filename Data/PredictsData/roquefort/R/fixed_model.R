FixedModel<-function(all.data,responseVar,fitFamily,fixedStruct,
                      randomStruct,saveVars=character(0),REML=TRUE,
                     optimizer="Nelder_Mead",maxIters=10000){
  
  call.best<-construct_call(responseVar,fixedStruct,randomStruct)
  
  allTerms<-unlist(strsplit(call.best,"[+]"))
  allTerms<-unlist(strsplit(allTerms,"[-]"))
  allTerms<-unlist(strsplit(allTerms,"[~]"))
  allTerms<-unlist(strsplit(allTerms,"[(]"))
  allTerms<-unlist(strsplit(allTerms,"[)]"))
  allTerms<-unlist(strsplit(allTerms,"[|]"))
  allTerms<-unlist(strsplit(allTerms,"[:]"))
  allTerms<-unlist(strsplit(allTerms,"[*]"))
  allTerms<-unique(allTerms)
  allTerms<-gsub(" ","",allTerms)
  allTerms<-allTerms[allTerms!=""]
  allTerms<-allTerms[allTerms!="1"]
  allTerms<-allTerms[allTerms!="0"]
  allTerms<-allTerms[!grepl("poly",allTerms)]
  
  allTerms<-gsub("[,][0-9]","",allTerms)
  
  allTerms<-c(allTerms,saveVars)
  
  allTerms<-unique(allTerms)
  
  model.data<-subset(all.data,select=c(allTerms))
  model.data<-na.omit(model.data)
  
  if (fitFamily=="gaussian"){
    eval(substitute(m<-lmer(cb,data=model.data,lmerControl(optimizer = optimizer,
                                                           optCtrl = list(maxfun=maxIters)),
                            REML=REML),
                    list(cb=call.best)))
  } else {
    eval(substitute(m<-glmer(cb,family=fitFamily,data=model.data,
                             control=glmerControl(optimizer = optimizer,
                                                  optCtrl = list(maxfun=maxIters))),
                    list(cb=call.best)))
  }
  
  
  
  return(list(model=m,data=model.data))
  
  
}