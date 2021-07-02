# Useful bits and bobs
.Log <- function(...) {
  if(getOption('PREDICTSverbose', TRUE)) cat(...)
}

SubHeading <- function(index, headings=paste0(c(letters,LETTERS), ')'), side=2,
                       at=par('usr')[4], line=2.5, las=1, ...) {
  # Plots a subheading
  mtext(headings[index], side=side, at=at, line=line, las=las, ...)
  invisible (1+index)
}

DegreeCellAreaKM <- function(lat, height, width) {
  # Returns the area in km squared of a grid cell in degrees of arc
  # lat - the latitudinal centre of the cell
  # height, width - the size of the grid cell in degrees

  # TODO Unit test
  # TODO Reference for this method

  radians <- function(theta) theta*pi/180.0

  # Convert the latitude into radians
  lat.rad <- radians(lat)

  # The equatorial and polar radii of the Earth in km
  eq.radius <-  6378137
  pol.radius <- 6356752.3142

  # Calculate cell area
  angular.eccentricity <- acos(radians(pol.radius/eq.radius))
  ecc.sq <- sin(radians(angular.eccentricity))^2
  flattening <- 1-cos(radians(angular.eccentricity))
  temp.val <- (eq.radius*cos(lat.rad))^2+(pol.radius*sin(lat.rad))^2
  m.phi <- ((eq.radius*pol.radius)^2)/(temp.val^1.5)
  n.phi <- (eq.radius^2)/sqrt(temp.val)
  lat.length <- pi/180*m.phi/1000
  long.length <- pi/180*cos(lat.rad)*n.phi/1000
  return (lat.length*height*long.length*width)
}

StripWhitespace <- function (v) {
  # Returns v with whitespace stripped from the start and the end
  return(gsub("^\\s+|\\s+$", "", v, perl = TRUE))
}

FormatThousands <- function(v, big.mark=",", ...) {
  # Returns v formatted as a string with a mark evert thousand
  return (StripWhitespace(format(v, big.mark=big.mark, ...)))
}

FormatPercent <- function(v) {
  # Returns v formatted as a string to 2 decimal places with a percentage sign
  f <- sprintf('%.2f%%', v)
  f[v>0 & v<0.005] <- '<0.01%'
  return (f)
}

Capitalize <- function(v) {
  # Returns the string v with the first letter in uppercase and the remainder 
  # of the letters in lowercase
  return (paste(toupper(substr(v, 1, 1)), 
                tolower(substr(v, 2, nchar(v))), 
                sep=''))
}

ParseRangeOrValue <- function(v, print.errors=TRUE) {
    # Takes a character vector v - values should be empty, a number with an
    # optional decimal point or two numbers, separated by a -, both with an 
    # optional decimal point. Returns a data.frame with three columns: lower, 
    # upper and value.

    ok <- grepl('(^\\s*$)|(^\\s*\\d+\\.?\\d*\\s*$)|(^\\s*\\d+\\.?\\d*\\s*-\\s*\\d+\\.?\\d*\\s*$)', v)
    if(!all(ok)) {
        if(print.errors) {
          cat('Bad values:\n')
          print(cbind(paste('row', which(!ok)), as.character(v[!ok])))
        }
        stop('Bad values')
    }

    contains.hyphen <- grepl('-', v)
    split <- strsplit(as.character(v), '-')

    # single.value will be a vector of logicals, containing TRUE where v
    # contains a single value, FALSE where v contains a hyphen.
    single.value <- 1==sapply(split, length)

    res <- matrix(NA, ncol=3, nrow=length(v))
    res[!single.value,1] <- as.numeric(sapply(split[!single.value], '[', 1))

    res[!single.value,2] <- as.numeric(sapply(split[!single.value], '[', 2))

    # If a range, both lower and upper should be given
    # lower should be <= upper or both should be NA
    ok <- res[,1]<res[,2] | (is.na(res[,1]) & is.na(res[,1]))
    if(!all(ok)) {
        if(print.errors) {
          cat('Bad values:\n')
          print(cbind(paste('row', which(!ok)), as.character(v[!ok])))
        }
        stop('Bad values')
    }

    res[single.value,3] <- as.numeric(sapply(split[single.value], '[', 1))

    return (res)
}

MeanOfRangeOrValue <- function(v, ...) {
    res <- ParseRangeOrValue(v, ...)
    v <- res[,3]
    # Compute mid-point of range
    v[is.na(v)] <- apply(res[is.na(v),1:2,drop=FALSE], 1, mean)
    return(v)
}

DiffOfRangeOrValue <- function(v, ...) {
    res <- ParseRangeOrValue(v, ...)
    v <- res[,3]
    # Compute difference of range
    v[is.na(v)] <- res[is.na(v),2]-res[is.na(v),1]
    return(v)
}

AddAlpha <- function(colour, alpha=0xcc, maxColorValue=0xff) {
    rgb <- col2rgb(colour)
    return (apply(rgb, 2, function(v) rgb(v[1],v[2],v[3], alpha, maxColorValue=maxColorValue)))
}
