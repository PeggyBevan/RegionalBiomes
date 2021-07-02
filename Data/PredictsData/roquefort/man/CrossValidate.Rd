\name{CrossValidate}
\alias{CrossValidate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross validate a mixed-effects model}
\description{
Performs leave-one-out cross validation of a mixed effects model, where entire levels in a random intercept are dropped together
}
\usage{
CrossValidate(model,nFolds,divFactor,fitFamily,data=NULL,fixedFactors=
                          character(0),fixedTerms=list(),
                        fixedInteractions=character(0),
                        fitInteractions=FALSE,randomStruct)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The model to cross-validate. Can be the full structure returned by the FixedModel and ModelSelect function in roquefort, or a standard lmer or glmer model from lme4
}
  \item{nFolds}{
The number of subsets to divide the data into, -1 gives all levels of the specified random factor
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
Factors to include in the model during cross-validation
}
  \item{fixedTerms}{
Continuous effects to include in the model during cross-validation
}
  \item{fixedInteractions}{
Interaction terms to include in the model during cross-validation  
}
  \item{fitInteractions}{
Whether to include all pairwise interaction terms in the model during cross-validation  
}
  \item{randomStruct}{
The random-effects structure to be used in the model for cross-validation. Should be the same as the random-effects structure in the original model  
}
}
\details{
}
\value{
all.results: a data frame of coefficient estimates with rows for each fold in the cross-validation and columns for each possible coefficient

confidence: 95-percent confidence intervals around the estimated coefficients across folds

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

# Cross-validate the model for all levels of SS
cv.m1<-CrossValidate(model=m1,nFolds=-1,divFactor="SS",fitFamily="poisson",fixedFactors="LandUse",randomStruct="(1|SS)+(1|SSBS)")

# 10-fold cross-validation
cv.m1<-CrossValidate(model=m1,nFolds=10,divFactor="SS",fitFamily="poisson",fixedFactors="LandUse",randomStruct="(1|SS)+(1|SSBS)")

# An example glmer model
m2<-glmer(Species_richness~LandUse+(1|SS)+(1|SSBS),family="poisson",data=sites)

cv.m2<-CrossValidate(model=m2,nFolds=-1,divFactor="SS",fitFamily="poisson",data=sites,fixedFactors="LandUse",randomStruct="(1|SS)+(1|SSBS)")

}