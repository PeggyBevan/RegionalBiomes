DropInvalidMetricsAndMethods <- function(extract) {
  # Drop rows with diversity metrics or sampling methods that are not valid
  keep <- extract$Diversity_metric_is_valid & extract$Sampling_method_is_valid
  .Log('Dropping', sum(!keep), 'rows with metrics that are not valid\n')
  extract <- droplevels(extract[keep,])
  .Log(paste(nrow(extract), 'rows remain\n'))
  return (extract)
}

CorrectSamplingEffort <- function(diversity) {
  # Assume that any missing values of sample effort mean equal sampling 
  # effort was applied.
  missing <- is.na(diversity$Sampling_effort)
  .Log('Correcting', sum(missing), 'missing sampling effort values\n')
  diversity$Sampling_effort[missing] <- 1
  
  # TODO Check this logic with Tim
  # Rescale sampling effort to have a maximum of 1 in each source
  .Log('Rescaling sampling effort\n')
  diversity$Sampling_effort <- do.call('c', tapply(diversity$Sampling_effort, 
                                                   diversity$SS, function(se) return (se/max(se))))
  
  # Correct for sensitivity to sampling effort
  sensitive <- diversity$Diversity_metric_is_effort_sensitive
  .Log('Correcting', sum(sensitive), 'values for sensitivity to sampling', 
       'effort\n')
  diversity$Measurement[sensitive] <- diversity[sensitive,'Measurement'] / 
    diversity[sensitive,'Sampling_effort']
  return (diversity)
}
