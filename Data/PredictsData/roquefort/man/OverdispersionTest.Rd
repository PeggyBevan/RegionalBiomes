\name{OverdispersionTest}
\alias{OverdispersionTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test for over-dispersion}
\description{
Tests for over-dispersion in the residuals of a mixed-effects model
}
\usage{
OverdispersionTest(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A mixed-effects model
  }
}
\details{
Compares residual deviance and residual degrees of freedom of a mixed-effects model to a chi-sqaure distribution. Adapted from code taken from the GLMM Wiki (see references) 
}
\value{
Reports the residual deviance, residual degrees of freedom and the ratio of these, and also a the P-value of a chi-square test comparing the residual deviance and degrees of freedom to a chi-square distribution
}
\references{
http://glmm.wikidot.com/faq
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
randomStruct="(1|SS)+(1|SSB)")

# Perform test for over-dispersion
OverdispersionTest(m1$model)

# Add an observation-level (in this case site-level) random effect to the model
m2<-FixedModel(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",
fixedStruct="LandUse+poly(LogHPD,2)+LandUse:poly(LogHPD,2)",
randomStruct="(1|SS)+(1|SSB)+(1|SSBS)")

# Test for overdispersion again
OverdispersionTest(m2$model)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
    