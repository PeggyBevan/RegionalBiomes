
\name{RandomsDiagnostics}
\alias{RandomsDiagnostics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random-effect diagnostics}
\description{
Performs diagnostic tests on modelled random effects from the models: currently only QQ plots 
}
\usage{
RandomsDiagnostics(model,randomName,fname,wDir=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The model whose random effects will be checked
}
  \item{randomName}{
The name of random effect for which to plot diagnostics
}
  \item{fname}{
The name of the file in which to save the QQ plots
}
  \item{wDir}{
The directory in which to save the plots, defaults to the current working directory
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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

# Perform random-effects diagnostics
RandomsDiagnostics(model=m1$model,randomName="SS",
fname="RandomsDiagnosticsExample.pdf")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
