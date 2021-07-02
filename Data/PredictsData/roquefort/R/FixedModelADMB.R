FixedModelADMB<-function(all.data,responseVar,fitFamily,fixedStruct,
                     randomStruct,saveVars=character(0)){
  
  call.best<-construct_call_admb(responseVar,fixedStruct)
  
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
  
  model.data<-subset(all.data,select=c(allTerms,"SS","SSB","SSBS"))
  model.data<-na.omit(model.data)
  
  eval(substitute(m<-glmmadmb(formula = cb,data=model.data,family = fitFamily,random = randomStruct),
                  list(cb=formula(call.best))))
  
  return(list(model=m,data=model.data))
  
  
}