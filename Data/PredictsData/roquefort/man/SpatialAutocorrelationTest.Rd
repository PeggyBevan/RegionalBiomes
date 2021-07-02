\name{SpatialAutocorrelationTest}
\alias{SpatialAutocorrelationTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial autocorrelation test (individual studies)}
\description{
Tests for spatial autocorrelation in the residuals of a PREDICTS model. Tests residuals for each PREDICTS study separately
}
\usage{
SpatialAutocorrelationTest(model,all.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A PREDICTS model with study identity (SS) as one of the random intercepts
  }
  \item{all.data}{
  The original dataset used to build the model, with 'Longitude' and 'Latitude' columns
  }
}
\details{
}
\value{
A list containing 4 items:

studies: all of the studies for which a Moran's test could be performed

I: the Moran's I statistic of all successful tests

P: the P-value of all successful tests

failed: the studies for which the Moran's test failed for some reason
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
sites$UseIntensity<-sites$Use_intensity
sites$UI<-factor(paste(sites$LandUse,sites$Use_intensity))
sites$UI<-relevel(sites$UI,ref="Natural Minimal use")

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

# Build a model
m1<-FixedModel(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",
fixedStruct="LandUse+poly(LogHPD,2)+LandUse:poly(LogHPD,2)",
randomStruct="(1|SS)+(1|SSB)",saveVars="SSBS")

# Perform spatial autocorrelation test on the model residuals
sat<-SpatialAutocorrelationTest(model=m1,all.data=sites)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
    