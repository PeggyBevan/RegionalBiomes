\name{DisaggregateMixedModel}
\alias{DisaggregateMixedModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Disaggregate a mixed-effects model by random factor level}
\description{
Runs a model for each level in a specified random factor and calculates the variation in estimated coefficient values
}
\usage{
DisaggregateMixedModel(model,divFactor,fitFamily,data=NULL,fixedFactors=
                          character(0),fixedTerms=list(),
                        fixedInteractions=character(0),
                        fitInteractions=FALSE,randomStruct)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The model to disaggregate. Can be the full structure returned by the FixedModel and ModelSelect function in roquefort, or a standard lmer or glmer model from lme4
}
  \item{divFactor}{
The random factor to divide the data by
}
  \item{fitFamily}{
The family of GLMER that was used for the model  
}
  \item{data}{
The data used to fit the model. Not necessary if the FixedModel or ModelSelect functions in roquefort were used to generate the model
}
  \item{fixedFactors}{
Factors to include in the model during disaggregation
}
  \item{fixedTerms}{
Continuous effects to include in the model during disaggregation
}
  \item{fixedInteractions}{
Interaction terms to include in the model during disaggregation
}
  \item{fitInteractions}{
Whether to include all pairwise interaction terms in the model during disaggregation  
}
  \item{randomStruct}{
The random-effects structure to be used in the model for disaggregation. Remeber to remove the grouping factor used for disaggregation. 
}
}
\details{
}
\value{
all.results: a data frame of coefficient estimates with rows for each random factor level and columns for each possible coefficient

confidence: 95-percent confidence intervals around the estimated coefficients across factor levels

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

# An example roquefort model
m1<-FixedModel(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",
fixedStruct="LandUse",
randomStruct="(1|SS)+(1|SSBS)")

# Perform disaggregation across levels of SS
da.m1<-DisaggregateMixedModel(model=m1,divFactor="SS",fitFamily="poisson",fixedFactors="LandUse",randomStruct="(1|SSBS)")

# An example glmer model
m2<-glmer(Species_richness~LandUse+(1|SS)+(1|SSBS),family="poisson",data=sites)

da.m2<-DisaggregateMixedModel(model=m2,divFactor="SS",fitFamily="poisson",data=sites,fixedFactors="LandUse",randomStruct="(1|SSBS)")

}