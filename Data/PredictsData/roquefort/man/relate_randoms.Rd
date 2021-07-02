\name{relate_randoms}
\alias{relate_randoms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Relate random effects to a factor}
\description{
Relates modelled random effects to an associated factor
}
\usage{
relate_randoms(model,randomName,data,colName,figurePath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The model from which to get the random effects
}
  \item{randomName}{
The name of random effect to use
}
  \item{data}{
The data frame that was used to fit the model
}
  \item{colName}{
The name of the column containing the factor to relate the random effects to
}
  \item{outDir}{
Optional directory to write the graph to. If not specified, then the graph will be written to the native graphics device
}
}
\details{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
