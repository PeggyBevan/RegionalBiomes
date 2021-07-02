\name{CompositionalSimilarityData}
\alias{CompositionalSimilarityData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate independent comparisons of compositional similarity among sites}
\description{
  Calculates compositional similarity within pairs of sites in a study, and generates independent comparisons
}
\usage{
  CompositionalSimilarity(data,metric,nIters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A PREDICTS measurement-level data frame
  }
  \item{metric}{
    The metric of compositional similarity to use. One of "Sor" (Sorensen), "SorAbd" (abundance-weighted Sorensen) and "SorCorr" (Sorensen corrected for incompleteness of sampling)
  }
  \item{nIters}{
    The number of random independent comparisons to create
  }
}
\details{
}
\value{
data: a list of length equal to the number of iterations specified, where each is a data fram of independent comparisons

counts: the number of sites considered for each land use
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
  # Load example PREDICTS data for the UK
  data(diversity)
  
  # Calculate compositional similarity using the Sorensen Index
  cd<-CompositionalSimilarityData(data=diversity,metric="Sor",nIters=5)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
