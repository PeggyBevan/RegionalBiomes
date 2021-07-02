\name{ModelSelect}
\alias{ModelSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Backward stepwise selection of fixed effects}
\description{
Performs backward stepwise selection of fixed effects in a model. Tests interaction terms first, and then drops them to test main effects. Main effects that are part of interaction terms will be retained, regardless of their significance as main effects
}
\usage{
ModelSelect(all.data,responseVar,fitFamily,fixedFactors=
                         character(0),fixedTerms=list(),
                       fixedInteractions=character(0),
                       randomStruct,siteRandom=FALSE,
                       fitInteractions=FALSE,verbose=FALSE,
                       otherRandoms=character(0),saveVars=character(0),
                       optimizer="Nelder_Mead")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all.data}{
A data frame containing the response variable, all fixed effects to be considered, and all terms in the specified random-effects structure. Only works with PREDICTS study, block and site designations and works best with the data frames returned by SiteMetrics
}
  \item{responseVar}{
The response variable to fit in the model
}
  \item{fitFamily}{
The family to use for the generalized linear mixed effects model
}
  \item{fixedFactors}{
The fixed-effect factors to consider in the model
}
  \item{fixedTerms}{
The fixed-effect continuous variables to consider in the model
}
  \item{fixedInteractions}{
Specific interaction terms to consider in the model
}
  \item{randomStruct}{
The random-effects structure to use  
}
  \item{siteRandom}{
Whether to consider a site-level random effect in the model. Default is FALSE
}
  \item{fitInteractions}{
Whether to fit all two-way interactions between the fixed effects in the model. Default is FALSE
}
  \item{verbose}{
Whether to report progress in detail
}
  \item{otherRandoms}{
Any random intercepts besides standard PREDICTS study, block and site designations that are contained in the specified random-effects structure  
}
  \item{saveVars}{
Any variables in the original data frame to retain in the model data frame for later analysis  
}
  \item{optimizer}{
The GLMER optimizer to use. Options are 'Nelder_Mead' (the default) and 'bobyqa'
}
}
\details{
With the specified random-effects structures, starts with the most complex fixed-effects structure and performs backward stepwise selection to obtain the minimum adequate model. Comparison of the fit of different models is based on likelihood-ratio tests. Interaction terms are tested first, and then removed to test main effects. All main effects that are part of significant interaction terms are retained in the final model regardless of their significance as main effects. ML estimation is used during selection of terms, and then REML is used for fitting the final model.
}
\value{
model: the final minimum adequate model

data: the dataset used in fitting the models, i.e. the original data frame with any rows containing NAs removed

stats: a table of statistics relating to each term considered

final.call: the call used to generate the final model
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

# Simple model of species richness as a function of land use
# First, compare random-effects structures
cr<-CompareRandoms(dataset=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors="LandUse",siteRandom=TRUE,verbose=TRUE)
# Then perform model selection
m1<-ModelSelect(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors="LandUse",
randomStruct=cr$best.random,siteRandom=TRUE,verbose=TRUE,
saveVars="Latitude")

# Model with an interactive effect of land use and human population density
cr<-CompareRandoms(dataset=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors=c("LandUse"),
fixedTerms=list(LogHPD=1),
fixedInteractions="LandUse:poly(LogHPD,1)",siteRandom=TRUE,verbose=TRUE)
m1<-ModelSelect(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors="LandUse",
fixedTerms=list(LogHPD=1),fixedInteractions="LandUse:poly(LogHPD,1)",
randomStruct=cr$best.random,siteRandom=TRUE,verbose=TRUE)

# Model with an interactive effect of land use and use-intensity
# Note, that the interaction must be fitted by pasting land use
# and use-intensity together to avoid rank deficiencies
sites$UI<-factor(paste(sites$LandUse,sites$Use_intensity))
sites$UI<-relevel(sites$UI,ref="Natural Minimal use")
cr<-CompareRandoms(dataset=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors=c("UI"),siteRandom=TRUE,verbose=TRUE)
m1<-ModelSelect(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors="UI",
randomStruct=cr$best.random,siteRandom=TRUE,verbose=TRUE)

# Model of Simpson diversity (Gaussian family)
sites$LogSimpson<-log(sites$Simpson_diversity)
cr<-CompareRandoms(dataset=sites,responseVar="LogSimpson",
fitFamily="gaussian",fixedFactors=c("UI"),fixedTerms=list(LogHPD=1),
fixedInteractions="LandUse:poly(LogHPD,1)",siteRandom=FALSE,verbose=TRUE)
m1<-ModelSelect(all.data=sites,responseVar="LogSimpson",
fitFamily="gaussian",fixedFactors="UI",fixedTerms=list(LogHPD=1),
fixedInteractions="LandUse:poly(LogHPD,1)",
randomStruct=cr$best.random,siteRandom=FALSE,verbose=TRUE)

}