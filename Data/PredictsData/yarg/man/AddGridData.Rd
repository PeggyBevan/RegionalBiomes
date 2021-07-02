\name{AddGridData}
\alias{AddGridData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adds data from maps to a data frame}
\description{
Adds data from a map loaded as a SpatialGridDataFrame in R to a data frame with coordinate data
}
\usage{
AddGridData(gridData,dataFrame,columnName,silent=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gridData}{
The map to get data from 
}
  \item{dataFrame}{
The data frame to add the data to. Must have columns called 'Longitude' and 'Latitude'
}
  \item{columnName}{
The name to give the new column of data in the data frame
}
  \item{silent}{
Whether to suppress progress reporting. Default is FALSE.  
}
}
\details{
}
\value{
The data frame with a new column containing corresponding values from the map
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

}