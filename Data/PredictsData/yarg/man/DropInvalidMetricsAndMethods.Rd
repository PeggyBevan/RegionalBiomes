\name{DropInvalidMetricsAndMethods}
\alias{DropInvalidMetricsAndMethods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Drops records with invalid metrics and units}
\description{
Drops any records that have invalid combinations of diversity metric and unit
}
\usage{
DropInvalidMetricsAndMethods(extract)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{extract}{
The PREDICTS database extract from which to drop invalid records
}
}
\details{
}
\value{
A PREDICTS database extract without records with invalid diversity metrics and units
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

# Remove any invalid metrics and units
diversity<-DropInvalidMetricsAndMethods(diversity)

}