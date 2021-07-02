\name{CompositionalSimilarity}
\alias{CompositionalSimilarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate compositional similarity among land uses}
\description{
  Calculates compositional similarity within pairs of sites in a study, and averages across pairs of land uses
}
\usage{
  CompositionalSimilarity(data,metric)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A PREDICTS measurement-level data frame
  }
  \item{metric}{
    The metric of compositional similarity to use. One of "Sor" (Sorensen), "SorAbd" (abundance-weighted Sorensen) and "SorCorr" (Sorensen corrected for incompleteness of sampling)
  }
}
\details{
}
\value{
cd: the matrix of average similarities in pairs of land uses

n: the number of site comparisons considered for each pair of land uses

studies: the studies considered in generating the estimates
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
  cd<-CompositionalSimilarity(data=diversity,metric="Sor")
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
