CompareRandoms<-function(dataset,responseVar,fitFamily,fixedFactors=character(0),
                          fixedTerms=list(),fixedInteractions=character(0),
                          siteRandom=FALSE,blockRandom=TRUE,
                         otherRandoms=character(0),
                          fitInteractions=FALSE,verbose=FALSE,
                         optimizer="Nelder_Mead",randomSlopes=TRUE,
                         specificRandomSlopes = NULL,
                         maxIters=10000){
  
  if ((length(fixedInteractions)>0) & (fitInteractions)){
    stop("Error: specifying particular interactions and all two-way interactions will not work!")
  }
  
  inter.split<-unlist(strsplit(fixedInteractions,":"))
  inter.split<-gsub("poly[(]","",inter.split)
  inter.mains<-gsub("[,][:0-9:][)]","",inter.split)
  
  # Subset the data frame, retaining just the relevant columns
  if (("UI" %in% fixedFactors) & (!("LUTax" %in% fixedFactors))){
    model.data<-subset(dataset,select=c("SS","SSB","SSBS",fixedFactors,
                                         "LandUse","UseIntensity",
                                         names(fixedTerms),responseVar,
                                        otherRandoms))
  } else if (("LUTax" %in% fixedFactors) & (!("UI" %in% fixedFactors))){
    model.data<-subset(dataset,select=c("SS","SSB","SSBS",fixedFactors,
                                         "LandUse","Taxon",
                                         names(fixedTerms),responseVar,
                                        otherRandoms))
  } else if (("UI" %in% fixedFactors) & ("LUTax" %in% fixedFactors)){
    model.data<-subset(dataset,select=c("SS","SSB","SSBS",fixedFactors,
                                         "LandUse","UseIntensity","Taxon",
                                         names(fixedTerms),responseVar,
                                        otherRandoms))
  } else {
    model.data<-subset(dataset,select=c("SS","SSB","SSBS",fixedFactors,
                                         names(fixedTerms),responseVar,
                                        otherRandoms))
  }
  model.data<-na.omit(model.data)
  cat<-sapply(model.data,is.factor)
  model.data[cat]<-lapply(model.data[cat],factor)
  
  # Create a list to store the results
  results<-list(ranef=character(),AIC=numeric())
  
  # Construct the fixed-effects part of the original model call
  fixedStruct<-""
  for (i in 1:length(fixedFactors)){
    fixedStruct<-paste(fixedStruct,fixedFactors[i],sep="")
    if ((i != length(fixedFactors)) | (length(fixedTerms)>0) | 
          ((length(fixedTerms)==0) & (
            length(fixedInteractions)>0))){
      fixedStruct<-paste(fixedStruct,"+",sep="")
    }
  }
  if (length(fixedTerms)>0){
    for (i in 1:length(fixedTerms)){
      fixedStruct<-paste(fixedStruct,"poly(",names(fixedTerms)[i],
                         ",",fixedTerms[i],")",sep="")
      if ((i != length(fixedTerms)) | (length(fixedInteractions)>0)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }
  
  if (fitInteractions) fixedStruct<-paste("(",fixedStruct,")^2",sep="")
  
  if (length(fixedInteractions)>0){
    for (i in 1:length(fixedInteractions)){
      fixedStruct<-paste(fixedStruct,fixedInteractions[i],sep="")
      if (i != length(fixedInteractions)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }
  
  
  if (("UI" %in% fixedFactors) & (!("LUTax" %in% fixedFactors))){
    fixedFactors<-fixedFactors[-which(fixedFactors=="UI")]
    fixedFactors<-c("LandUse","UseIntensity",fixedFactors)
  } else if (("LUTax" %in% fixedFactors) & (!("UI" %in% fixedFactors))){
    fixedFactors<-fixedFactors[-which(fixedFactors=="LUTax")]
    fixedFactors<-c("LandUse","Taxon",fixedFactors)
  } else if (("UI" %in% fixedFactors) & ("LUTax" %in% fixedFactors)){
    fixedFactors<-fixedFactors[-which(fixedFactors=="UI")]
    fixedFactors<-fixedFactors[-which(fixedFactors=="LUTax")]
    fixedFactors<-c("LandUse","UseIntensity","Taxon",fixedFactors)
  } else {
  }
  
  print(paste("Fixed structure:",fixedStruct))
  
  # Set the original value for the best model AIC
  best.aic<-Inf
  
  # First, try fitting just study as a random effect
  print(paste("Comparing random structure 1 of ",length(fixedFactors)+
                length(fixedTerms)+length(otherRandoms)+
                2+siteRandom,sep=""))
  
  new.random<-"(1|SS)"
  best.random<-new.random
  
  # Construct the complete model call
  new.call<-construct_call(responseVar,fixedStruct,new.random)
  if (verbose) print(new.call)
  
  # Run the model
  if (fitFamily=="gaussian"){
    mod<-try(lmer(new.call,data=model.data,REML=FALSE,
                  lmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
  } else {
    mod<-try(glmer(new.call,family=fitFamily,data=model.data,
                   control=glmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
  }
  
  
  # Check whether this model has the best AIC value and update the results
  if ((class(mod)!="try-error")){
    if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
      results$ranef<-c(results$ranef,new.random)
      results$AIC<-c(results$AIC,AIC(mod))
      
      if(AIC(mod)<best.aic){
        best.aic<-AIC(mod)
        best.random<-new.random
      }
    }
  }
  
  if (blockRandom){
    # Try adding block to the random structure
    print(paste("Comparing random structure 2 of ",length(fixedFactors)+
                  length(fixedTerms)+length(otherRandoms)+
                  2+siteRandom,sep=""))
    
    new.random<-paste(best.random,"+(1|SSB)",sep="")
    
    # Construct the complete model call
    new.call<-construct_call(responseVar,fixedStruct,new.random)
    if (verbose) print(new.call)
    
    # Run the model
    if (fitFamily=="gaussian"){
      mod<-try(lmer(new.call,data=model.data,REML=FALSE,
                    lmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
    } else {
      mod<-try(glmer(new.call,family=fitFamily,data=model.data,
                     control=glmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
    }
    
    
    # Check whether this model has the best AIC value and update the results
    if ((class(mod)!="try-error")){
      if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
        results$ranef<-c(results$ranef,new.random)
        results$AIC<-c(results$AIC,AIC(mod))
        
        if(AIC(mod)<best.aic){
          best.aic<-AIC(mod)
          best.random<-new.random
        }
      }
    }
    
  }
  
  # If the consideration of a site-level random has been specified, then try fitting it
  if(siteRandom){
    print(paste("Comparing random structure 3 of ",length(fixedFactors)+
                  length(fixedTerms)+length(otherRandoms)+
                  2+siteRandom,sep=""))
    
    new.random<-paste(best.random,"+(1|SSBS)",sep="")
    if (verbose) print(new.call)
    
    # Construct the complete model call
    new.call<-construct_call(responseVar,fixedStruct,new.random)
    
    # Run the model
    if (fitFamily=="gaussian"){
      mod<-try(lmer(new.call,data=model.data,REML=FALSE,
                    lmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
    } else {
      mod<-try(glmer(new.call,family=fitFamily,data=model.data,
                     control=glmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
    }
    
    
    # Check whether this model has the best AIC value and update the results
    if ((class(mod)!="try-error")){
      if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
        results$ranef<-c(results$ranef,new.random)
        results$AIC<-c(results$AIC,AIC(mod))
        
        if(AIC(mod)<best.aic){
          best.aic<-AIC(mod)
          best.random<-new.random
        }
      }
    }
    
    
  }
  
  # Try adding any other specified randoms
  i<-1
  for (re in otherRandoms){
    print(paste("Comparing random structure ",2+siteRandom+i," of ",
                length(fixedFactors)+length(fixedTerms)+
                  length(otherRandoms)+2+siteRandom,sep=""))
    
    new.random<-paste(best.random,"+(1|",re,")",sep="")
    if (verbose) print(new.call)
    
    # Construct the complete model call
    new.call<-construct_call(responseVar,fixedStruct,new.random)
    
    # Run the model
    if (fitFamily=="gaussian"){
      mod<-try(lmer(new.call,data=model.data,REML=FALSE,
                    lmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
    } else {
      mod<-try(glmer(new.call,family=fitFamily,data=model.data,
                     control=glmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
    }
    
    
    # Check whether this model has the best AIC value and update the results
    if ((class(mod)!="try-error")){
      if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
        results$ranef<-c(results$ranef,new.random)
        results$AIC<-c(results$AIC,AIC(mod))
        
        if(AIC(mod)<best.aic){
          best.aic<-AIC(mod)
          best.random<-new.random
        }
      }
    }
    
    i<-i+1
  }
  
  if (randomSlopes){
    # For each specified factor, try adding the factor in a random 
    # slopes-and-intercepts model
    
    SlopesToTest <- fixedFactors
    
    if(!is.null(SlopesToTest)){
      stopifnot(all(specificRandomSlopes %in% SlopesToTest))
      
      SlopesToTest <- SlopesToTest[SlopesToTest %in% specificRandomSlopes]
      
    }
    
    i<-1
    for (f in SlopesToTest){
      print(paste("Comparing random structure ",2+length(otherRandoms)+siteRandom+i," of ",
                  length(SlopesToTest)+length(fixedTerms)+
                    length(otherRandoms)+2+siteRandom,sep=""))
      
      new.random<-gsub("[|]SS)",paste("+",f,"|SS)",sep=""),best.random)
      
      # Construct the complete model call
      new.call<-construct_call(responseVar,fixedStruct,new.random)
      
      if (verbose) print(new.call)
      
      # Run the model
      if (fitFamily=="gaussian"){
        mod<-try(lmer(new.call,data=model.data,REML=FALSE,
                      lmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
      } else {
        mod<-try(glmer(new.call,family=fitFamily,data=model.data,
                       control=glmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
      }
      
      
      # Check whether this model has the best AIC value and update the results
      if ((class(mod)!="try-error")){
        if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
          results$ranef<-c(results$ranef,new.random)
          results$AIC<-c(results$AIC,AIC(mod))
          
          if(AIC(mod)<best.aic){
            best.aic<-AIC(mod)
            best.random<-new.random
          }
        }
      }
      
      
    }
    
    # For each specified continuous term, try adding the factor in a random 
    # slopes-and-intercepts model
    if (length(fixedTerms)>0){
      for (i in 1:length(fixedTerms)){
        print(paste("Comparing random structure ",2+siteRandom+length(otherRandoms)+length(fixedFactors)+i," of ",
                    length(otherRandoms)+length(fixedFactors)+length(fixedTerms)+2+siteRandom,sep=""))
        
        new.random<-gsub("[|]SS)",paste("+poly(",names(fixedTerms)[i],",",
                                        fixedTerms[i],
                                        ")|SS)",sep=""),best.random)
        
        # Construct the complete model call
        new.call<-construct_call(responseVar,fixedStruct,new.random)
        
        if (verbose) print(new.call)
        
        # Run the model
        if (fitFamily=="gaussian"){
          mod<-try(lmer(new.call,data=model.data,REML=FALSE,
                        lmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
        } else {
          mod<-try(glmer(new.call,family=fitFamily,data=model.data,
                         control=glmerControl(optimizer=optimizer,optCtrl=list(maxfun=maxIters))))
        }
        
        
        # Check whether this model has the best AIC value and update the results
        if ((class(mod)!="try-error")){
          if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
            results$ranef<-c(results$ranef,new.random)
            results$AIC<-c(results$AIC,AIC(mod))
            
            if(AIC(mod)<best.aic){
              best.aic<-AIC(mod)
              best.random<-new.random
            }
          }
        }
        
        
      }
      
    }
    
  }
  
  
  return(list(full.results=results,best.random=best.random))
  
}