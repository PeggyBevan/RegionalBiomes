compile_metric_registry<-function(registryPath){

  warning('This function is deprecated - units are now included in database ',
          'extracts')
  # Read in the metrics and units registry from file, and return
  metricsReg<-read.csv(registryPath)
  
  return(metricsReg)
  
}


