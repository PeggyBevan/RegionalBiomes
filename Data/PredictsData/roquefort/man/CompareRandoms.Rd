\name{CompareRandoms}
\alias{CompareRandoms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compare different random-effects structures}
\description{
Compares different possible random effects structures and selects the one that gives the lowest AIC value
}
\usage{
CompareRandoms(dataset, responseVar, fitFamily, fixedFactors = character(0), 
    fixedTerms = list(), fixedInteractions = character(0), siteRandom = FALSE, 
    otherRandoms = character(0), fitInteractions = FALSE, verbose = FALSE, 
    optimizer = "Nelder_Mead", randomSlopes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
A data frame. Only works with PREDICTS study, block and site designations and works best with the data frames returned by SiteMetrics
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
  \item{siteRandom}{
Whether to consider a site-level random effect in the model. Default is FALSE
}
  \item{otherRandoms}{
Any other random intercepts to consider  
}
  \item{fitInteractions}{
Whether to fit all two-way interactions between the fixed effects in the model. Default is FALSE
}
  \item{verbose}{
Whether to report progress in detail
}
  \item{optimizer}{
The optimizer to use in fitting the generalized linear mixed effects models. One of  "Nelder_Mead" (the default) or "bobyqa".  
}
  \item{randomSlopes}{
Whether to consider random slopes in the random-effects structures considered. Default is to consider random slopes.  
}
}
\details{
The most complex model in terms of fixed-effect structure is fitted - including two-way interactions if specified - 
and then different random effects (random intercepts, then random slopes and intercepts) are added, being retained if they lead to a decrease in model AIC.
}
\value{
A list containing two items:

full.results: a list with all random-effects structures considered ($ranef) and corresponding AIC values ($AIC), 

best.random: the random-effects structure with the lowest AIC value
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

# Compare random-effects structures for a simple model 
# of species richness as a function of land use
cr<-CompareRandoms(dataset=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors="LandUse",siteRandom=TRUE,verbose=TRUE)

# Try comparison for a model with land use and use-intensity
cr<-CompareRandoms(dataset=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors=c("LandUse","UseIntensity"),
siteRandom=TRUE,verbose=TRUE)

# Try comparison for a model with an interactive effect
# of land use and human population density
cr<-CompareRandoms(dataset=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors=c("LandUse"),
fixedTerms=list(LogHPD=1),
fixedInteractions="LandUse:poly(LogHPD,1)",siteRandom=TRUE,verbose=TRUE)

# Try comparison for a model with an interactive effect 
# of land use and use-intensity
# Note, that the interaction must be fitted by pasting land use
# and use-intensity together to avoid rank deficiencies
sites$UI<-factor(paste(sites$LandUse,sites$Use_intensity))
sites$UI<-relevel(sites$UI,ref="Natural Minimal use")
cr<-CompareRandoms(dataset=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors=c("UI"),siteRandom=TRUE,verbose=TRUE)

# Try comparison including biome as a random intercept
cr<-CompareRandoms(dataset=sites,responseVar="Species_richness",
fitFamily="poisson",fixedFactors=c("UI"),siteRandom=TRUE,
otherRandoms="Biome",verbose=TRUE)

}