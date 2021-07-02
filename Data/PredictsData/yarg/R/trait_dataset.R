# TraitDataset. An S3 class that represents values:
#   for a single trait
#   from a single source
#   with a single unit
#   for one or more taxa in one of more kingdoms

is.TraitDataset <- function(x) {
  return (inherits(x, 'TraitDataset'))
}

"[<-.TraitDataset" <- function(x, i, value) {
  if(!is.TraitDataset(x)) stop('Not a TraitDataset')
  stop("Can't assign to a TraitDataset")
}

"[[<-.TraitDataset" <- function(x, i, value) {
  if(!is.TraitDataset(x)) stop('Not a TraitDataset')
  stop("Can't assign to a TraitDataset")
}

'$<-.TraitDataset' <- function(object, x, value) {
  if(!is.TraitDataset(object)) stop('Not a TraitDataset')
  stop("Can't assign to a TraitDataset")
}

"names<-.TraitDataset" <- function(x, value) {
  if(!is.TraitDataset(x)) stop('Not a TraitDataset')
  stop("Can't assign to a TraitDataset")
}

"length<-.TraitDataset" <- function(x, value) {
  if(!is.TraitDataset(x)) stop('Not a TraitDataset')
  stop("Can't change length of a TraitDataset")
}

"levels<-.TraitDataset" <- function(x, value) {
  if(!is.TraitDataset(x)) stop('Not a TraitDataset')
  stop("Can't change levels of a TraitDataset")
}

"dim<-.TraitDataset" <- function(x, value) {
  if(!is.TraitDataset(x)) stop('Not a TraitDataset')
  stop("Can't change dim of a TraitDataset")
}

print.TraitDataset <- function(x, ...) {
  if(!is.TraitDataset(x)) stop('Not a TraitDataset')
  cat(FormatThousands(nrow(x$Values)), ' measurements of ', x$Trait, ' (', 
      x$Unit, ') from ', x$Reference, '\n', sep='')
  invisible(x)
}

summary.TraitDataset <- function(object, ...) {
  if(!is.TraitDataset(object)) stop('Not a TraitDataset')
  if(is.factor(object$Values$Value)) {
    return (table(object$Values$Value))
  } else {
    return (summary(object$Values$Value))
  }
}

plot.TraitDataset <- function(x, ...) {
  if(!is.TraitDataset(x)) stop('Not a TraitDataset')
  title <- paste(FormatThousands(nrow(x$Values)), 'measurements of', x$Trait, 
                 'from', x$Reference)
  if(is.factor(x$Values$Value)) {
    plot(x$Values$Value, xlab=x$Unit, main=title, ...)
  } else if(is.logical(x$Values$Value)) {
    barplot(table(x$Values$Value), xlab='', main=title, ...)
  } else {
    hist(x$Values$Value, xlab=x$Unit, main=title, ...)
    m <- median(x$Values$Value)
    abline(v=m, lty=2, ...)
    mtext(side=3, at=m, sprintf('median=%.2f', m), ...)
  }
}

TraitDataset <- function(reference, trait, unit, taxa, ranks, kingdoms, values) {
  stopifnot(1==length(reference) && is.character(reference) && 0<nchar(reference))
  stopifnot(1==length(trait) && is.character(trait) && 0<nchar(trait))
  stopifnot(1==length(unit) && is.character(unit) && 0<nchar(unit))
  stopifnot(0<length(taxa))
  stopifnot(is.character(taxa) || is.factor(taxa))
  stopifnot(all(nchar(as.character(taxa))>0))
  stopifnot(!any(duplicated(paste(kingdoms, taxa))))
  stopifnot(length(values)==length(taxa))
  stopifnot(is.numeric(values) || is.factor(values) || is.character(values) || is.logical(values))
  stopifnot(1==length(kingdoms) || length(kingdoms)==length(taxa))
  stopifnot(1==length(ranks) || length(ranks)==length(taxa))
  ranks <- EnsureFactorLevels(ranks, TaxonomicRanks())

  if(is.character(values)) values <- factor(values)

  values <- data.frame(Taxon=taxa, Rank=ranks, Kingdom=kingdoms, Value=values)

  to.remove <- is.na(values$Value) | is.infinite(values$Value)
  if(any(to.remove)) {
    .Log('Removing ', sum(to.remove), ' NA or Inf measurements of ', trait , 
         ' (', unit, ') from ', reference, '\n', sep='')
    values <- values[!to.remove,]
    # Drop levels for Kingdoms and Taxon, not for ranks or values
    values[,c('Kingdom','Taxon')] <- droplevels(values[,c('Kingdom','Taxon')])
  }

  stopifnot(nrow(values)>0)

  .Log('New TraitDataset of ', nrow(values), ' measurements of ', trait , ' (', 
       unit, ') from ', reference, '\n', sep='')
  self <- list(Reference=reference, 
               Trait=trait,
               Unit=unit,
               Values=values)
  class(self) <- c('TraitDataset', 'list')
  return(self)
}

TraitColname <- function(self) {
  # A character suitable for use as a column name
  if(!is.TraitDataset(self)) stop('Not a TraitDataset')
  n <- paste(self$Trait, self$Unit)
  n <- gsub('[- ]', '_', n, perl=TRUE)  # yarg standard is to use _ rather than .
  n <- gsub('[-\\(\\)\\[\\]\\*]*', '', n, perl=TRUE)
  return (make.names(n))  # Replace remaining illegal characters with .
}

.MatchedValues <- function(a, b, a.taxon, b.taxon, a.kingdom, b.kingdom) {
  # Add values in b to a, where not already present in a and where taxa and 
  # kingdoms match.
  use <- match(paste(a.kingdom, a.taxon), paste(b.kingdom, b.taxon))

  # Don't match values for taxa that have an existing value
  use[!is.na(a)] <- NA

  .Log(paste('Using', sum(!is.na(use)), 'values\n'))

  a[!is.na(use)] <- b[use[!is.na(use)]]
  return (a)
}

AddTraits <- function(values, ...) {
  # Add traits in ... to values. Arguments should be objects of class 
  # TraitDataset or lists of TraitDatasets.
  stopifnot('Kingdom' %in% colnames(values))

  # Taxon first, to catch infraspecies matches, then best-guess binomial
  candidate.cols <- c('Taxon','Best_guess_binomial')
  stopifnot(any(candidate.cols %in% colnames(values)))
  F <- function(trait, values) {
    if(!is.TraitDataset(trait)) stop('Not a TraitDataset')

    # Adds values in traits to values by matching taxon names

    .Log('Adding', trait$Trait, 'in units of', trait$Unit, 'from', 
         trait$Reference, '\n')

    col <- TraitColname(trait)
    if(!col %in% colnames(values)) {
      # Create empty column
      # TraitDataset guarantees that trait$Values$Value is either a numeric, 
      # logical or factor
      if(is.factor(trait$Values$Value)) {
        values[,col] <- factor(levels=levels(trait$Values$Value))
      } else {
        # Logical or numeric
        values[,col] <- NA
        class(values[,col]) <- class(trait$Values$Value)
      }
    }

    # Match on best buess binomial taxon before best-guess binomial so that 
    # we match infraspecies.
    for(match.col in intersect(candidate.cols, colnames(values))) {
      .Log('Matching using', match.col, '\n')
      values[,col] <- .MatchedValues(values[,col], trait$Values$Value, 
                                     values[,match.col], trait$Values$Taxon, 
                                     values$Kingdom, trait$Values$Kingdom)
    }
    return (values)
  }

  dots <- list(...)
  stopifnot(all(unlist(sapply(dots, class)) %in% c('list', 'TraitDataset')))

  for(x in dots) {
    if(is.TraitDataset(x)) {
      # Add traits from this dataset
      values <- F(x, values)
    } else {
      # Add traits from this list of TraitDatasets
      for(trait in x) {
        values <- F(trait, values)
      }
    }
  }

  return (values)
}

AddGenusAveragedTraits <- function(values, ...) {
  # Adds genus-level averaged values in ... to values. Arguments should be 
  # objects of class TraitDataset or lists of TraitDatasets.

  # It is tempting to match genera by using the first word of Parsed_name or of 
  # Best_guess_binomial but this is dangerous because we might use a name 
  # that was not intended as a genus but that matches a genus in a traits 
  # dataset.

  stopifnot('Kingdom' %in% colnames(values))
  stopifnot('Genus' %in% colnames(values))

  F <- function(trait, values) {
    if(!is.TraitDataset(trait)) stop('Not a TraitDataset')

    stopifnot('numeric'==class(trait$Values$Value))

    # Adds values in traits to values by matching taxon names
    .Log('Adding genus-averaged', trait$Trait, 'from', trait$Reference, 
         'in units of', trait$Unit, '\n')

    col <- TraitColname(trait)
    if(!col %in% colnames(values)) {
      # Create empty column
      values[,col] <- NA
    }

    if(any(is.na(values[,col]))) {
      # Genus-level means
      trait.genus <- sapply(strsplit(as.character(trait$Values$Taxon), ' '), '[', 1)
      genus.mean <- tapply(trait$Values$Value, trait.genus, mean, na.rm=TRUE)[trait.genus]

      values[,col] <- .MatchedValues(values[,col], genus.mean, 
                                     values$Genus, trait.genus, 
                                     values$Kingdom, trait$Values$Kingdom)
    }
    return (values)
  }

  dots <- list(...)
  stopifnot(all(unlist(sapply(dots, class)) %in% c('list', 'TraitDataset')))

  for(x in dots) {
    if(is.TraitDataset(x)) {
      # Add traits from this dataset
      values <- F(x, values)
    } else {
      # Add traits from this list of TraitDatasets
      for(trait in x) {
        values <- F(trait, values)
      }
    }
  }

  return (values)
}

AddGenusInterpolatedCategoricalTraits <- function(values,threshold=0.95, ...) {
  # Adds genus-level averaged values in ... to values. Arguments should be 
  # objects of class TraitDataset or lists of TraitDatasets.
  
  # It is tempting to match genera by using the first word of Parsed_name or of 
  # Best_guess_binomial but this is dangerous because we might use a name 
  # that was not intended as a genus but that matches a genus in a traits 
  # dataset.
  
  stopifnot('Kingdom' %in% colnames(values))
  stopifnot('Genus' %in% colnames(values))
  
  F <- function(trait, values) {
    if(!is.TraitDataset(trait)) stop('Not a TraitDataset')
    
    stopifnot('factor'==class(trait$Values$Value))
    stopifnot(1>=threshold)
    stopifnot(0.5<threshold)
    
    # Adds values in traits to values by matching taxon names
    .Log('Adding genus-averaged', trait$Trait, 'from', trait$Reference, 
         'in units of', trait$Unit, '\n')
    
    col <- TraitColname(trait)
    if(!col %in% colnames(values)) {
      # Create empty column
      values[,col] <- NA
    }
    
    if(any(is.na(values[,col]))) {
      # Genus-level majority values
      trait.genus <- sapply(strsplit(as.character(trait$Values$Taxon), ' '), '[', 1)
      genus.majority<-tapply(trait$Values$Value, trait.genus,function(x){
        
        trait.propns<-summary(x)/sum(summary(x))
        if(max(trait.propns>threshold)){
          return(names(trait.propns[trait.propns>threshold]))
        } else {
          return(NA)
        }
        
        
      })[trait.genus]
      
      values[,col] <- .MatchedValues(values[,col], genus.majority, 
                                     values$Genus, trait.genus, 
                                     values$Kingdom, trait$Values$Kingdom)
    }
    return (values)
  }
  
  dots <- list(...)
  stopifnot(all(unlist(sapply(dots, class)) %in% c('list', 'TraitDataset')))
  
  for(x in dots) {
    if(is.TraitDataset(x)) {
      # Add traits from this dataset
      values <- F(x, values)
    } else {
      # Add traits from this list of TraitDatasets
      for(trait in x) {
        values <- F(trait, values)
      }
    }
  }
  
  return (values)
}


InferTraitValues <- function(a, b, r.2.threshold=0.9, show=FALSE) {
  # Infers values of a from values of b by fitting an ordinary linear 
  # regression. a and b should be data.frames with columns Taxon and Value.

  if(!is.TraitDataset(a)) stop('Not a TraitDataset')
  if(!is.TraitDataset(b)) stop('Not a TraitDataset')

  stopifnot(a$Trait!=b$Trait)

  stopifnot(r.2.threshold<1 && r.2.threshold>0.5)

  stopifnot('numeric'==class(a$Values$Value))
  stopifnot('numeric'==class(b$Values$Value))

  .Log('Inferring', a$Trait, 'in units of', a$Unit, 'from', a$Reference, 
       'from', b$Trait, 'in units of', b$Unit, 'from', b$Reference, '\n')

  # Assemble a data.frame with columns Taxon, A and B
  values <- a$Values[,c('Taxon','Value')]
  colnames(values)[2] <- 'A'
  values$B <- b$Values[match(a$Values$Taxon, b$Values$Taxon), 'Value']

  # Drop values for which we do not have both A and B
  values <- na.omit(values)

  stopifnot(nrow(values)>1)

  # Infer values for taxa that are not in a
  infer.for.taxa <- !as.character(b$Values$Taxon) %in% as.character(a$Values$Taxon)

  .Log(paste('Inferring values for', sum(infer.for.taxa), 'values from', 
             'data for', nrow(values), 'taxa\n'))

  # Fit an ordinary linear regression
  m <- lm(A~B, data=values)

  # R-squared must be greater than or equal to threshold
  if(summary(m)$r.squared<r.2.threshold) {
    stop('R-squared is less than ', sprintf('%.2f', r.2.threshold))
  } else {
    .Log(paste('Inferring values based on a regression with an R-squared of',summary(m)$r.squared))
    # Create a data.frame of predictors from the taxa that we are interested in
    predictors <- b$Values[infer.for.taxa,'Value',drop=FALSE]
    colnames(predictors) <- 'B'

    # Predicted values
    inferred.values <- unname(predict(m, predictors))

    if(show) {
      plot(A~B, data=values, pch=19, col='grey')
      abline(m)
      points(predictors$B, inferred.values, pch=1)
    }

    reference <- paste('Inferred from', paste(b$Reference, b$Trait))
    return (TraitDataset(reference=reference,
                         trait=a$Trait,
                         unit=a$Unit,
                         taxa=droplevels(b$Values[infer.for.taxa,'Taxon']),
                         ranks=droplevels(b$Values[infer.for.taxa,'Rank']), 
                         kingdoms=droplevels(b$Values[infer.for.taxa,'Kingdom']), 
                         value=inferred.values))
  }
}
