\name{FixedModel}
\alias{FixedModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a fixed-structure model}
\description{
Fit a model with fixed random- and fixed-effects structures
}
\usage{
FixedModel(all.data,responseVar,fitFamily,fixedStruct,
                      randomStruct,saveVars=character(0),
                      REML=TRUE,optimizer="Nelder_Mead")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all.data}{
A data frame containing the response variable, and all fixed and random effects in the specified model structure
}
  \item{responseVar}{
The response variable to fit in the model
}
  \item{fitFamily}{
The family to use for the generalized linear mixed effects model
}
  \item{fixedStruct}{
The fixed effects to include, in the format of a glmer model-call   
}
  \item{randomStruct}{
The random effects to include, in the format of a glmer call  
}
  \item{saveVars}{
Any variables in the original data frame to retain in the model data frame for later analysis  
}
  \item{REML}{
Whether to use Restricted Maximum Likelihood for fitting the model. Alternative is simple maximum likelihoof. Default is to use REML  
}
  \item{optimizer}{
The GLMER optimizer to use. Options are 'Nelder_Mead' (the default) and 'bobyqa'
}
}
\details{
}
\value{
model: the model

data: the dataset used in to fit the model, i.e. the original data frame with any rows containing NAs removed
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

# An example model
m1<-FixedModel(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",
fixedStruct="LandUse+poly(LogHPD,2)+LandUse:poly(LogHPD,2)",
randomStruct="(1|SS)+(1|SSBS)")

}