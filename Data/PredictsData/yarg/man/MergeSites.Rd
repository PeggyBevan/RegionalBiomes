\name{MergeSites}
\alias{MergeSites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Merge nearby and similar sites}
\description{
Merge sites that have the same coordinates, predominant habitat, land-use intensity and sampling dates
}
\usage{
MergeSites(diversity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{diversity}{
A PREDICTS database extract
}
}
\details{
}
\value{
The PREDICTS database extract with updated site identifiers
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
# Load example PREDICTS database extract
data(diversity)

# Run the merge sites function
diversity<-MergeSites(diversity)


}