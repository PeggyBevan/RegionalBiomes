construct_call<-function(responseVar,fixedStruct,randomStruct){
  return(paste(responseVar,"~",fixedStruct,"+",randomStruct,sep=""))
}

construct_call_admb<-function(responseVar,fixedStruct){
  return(paste(responseVar,"~",fixedStruct,sep=""))
}