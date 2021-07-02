\name{CheckCollinearity}
\alias{CheckCollinearity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Check collinearity among variables}
\description{
Check collinearity among pairs of variables in a data frame
}
\usage{
  CheckCollinearity(data,allTerms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A data frame
  }
  \item{allTerms}{
    The terms among which to check collinearity
  }
  }
\details{
    }
    \references{
    }
    \author{
    Tim Newbold <Tim.Newbold@unep-wcmc.org>
    }
    \note{
    %%  ~~further notes~~
    }
    
    %% ~Make other sections like Warning with \section{Warning }{....} ~
    
    \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
    }
    \examples{
# Load example PREDICTS site-level data for the UK
data(sites)

# Load human population density map-header and data for the UK
data(hpd)

# Convert to a SpatialGridDataFrame object
hpd<-SpatialGridDataFrame(grid=GridTopology(
  cellcentre.offset=c(hpdHeader[3,2]+(hpdHeader[5,2]/2),
                      hpdHeader[4,2]++(hpdHeader[5,2]/2)),
  cellsize=c(hpdHeader[5,2],hpdHeader[5,2]),
  cells.dim=c(hpdHeader[1,2],hpdHeader[2,2])),
  data=data.frame(band1=c(t(hpdData))))

# Add human population density estimates to the sites data frame
sites<-AddGridData(gridData=hpd,dataFrame=sites,
                   columnName="HPD")
    
# Log estimates of human population density
sites$LogHPD<-log(sites$HPD+1)
  
# Test collinearity among land use, use-intensity and human population density
CheckCollinearity(data = sites,allTerms = c("LandUse","Use_intensity","LogHPD"))

    }
    % Add one or more standard keywords, see file 'KEYWORDS' in the
    % R documentation directory.
    \keyword{ ~kwd1 }
    \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
    