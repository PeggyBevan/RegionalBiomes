\name{SiteMetrics}
\alias{SiteMetrics}
\title{Site-Level Diversity Metrics}
\description{Compute site-level diversity metrics (species richness, 
Chao1-estimated species richness, total abundance and Simpson's Diversity 
Index) given a dataset that is in standard PREDICTS structure.}

\usage{
SiteMetrics(diversity, extra.cols=NULL, traits=NULL)
}

\arguments{
  \item{diverstiy}{a \code{data.frame} containing abundance, occurrence or 
  diversity records. Must be in the standard PREDICTS format.}
  \item{extra.cols}{the names of any columns in `diversity' that will be 
  included in the returned \code{data.frame}. Values within these columns 
  should be unique within sites.}
  \item{traits}{the names of columns in `diversity' for which 
  community-weighted mean values will be included in the returned 
  \code{data.frame}.}
}
\details{
Species richness is calculated simply as the number of observed taxa at a site.

Chao1-estimated species richness is calculated using the bias-corrected Chao1 
estimator (Gotelli & Colwell, 2010): Sest = Sobs + ((f1*(f1-1))/(2*(f2+1))),
where Sest is estimated species richness, Sobs is observed species 
richness, f1 is the number of taxa represented by one individual, and f2 
is the number of taxa represented by two individuals. Chao1-estimated species 
richness is computed only for studies for which the 
Diversity_metric_is_suitable_for_Chao column is \code{TRUE} and all 
values in the Measurement column are integers.

Total abundance is simply the sum of abundances of all taxa at a site, 
regardless of how abundance was measured.

Simpson's Diversity Index is calculated as 1/D, where D = sum((a/A)^2), where 
a is the abundance of a single taxon, and A is the total abundance of all taxa 
at a site.
}

\value{A \code{data.frame} with a row for each site nested within studies  
nested within sources.}

\references{
Gotelli & Colwell (2010). Estimating Species Richness. In Biological Diversity: 
Frontiers in Measurement and Assessment (Magurran, A. & McGill, B. J., eds.). 
Oxford University Press, Oxford, UK.
}

\author{Tim Newbold and Lawrence Hudson}

\seealso{\code{\link{data.frame}}}

\examples{
# BASIC EXAMPLE
# Read in PREDICTS database extract
data(diversity)

# Correct for sampling effort
diversity <- CorrectSamplingEffort(diversity)

# Merge sites that are very close to one another and in the same land-use type
diversity <- MergeSites(diversity)

# Calculate site-level diversity metrics
sites<-SiteMetrics(diversity=diversity,
                       extra.cols=c("SSB","SSBS","Biome"),
                       sites.are.unique=TRUE)

# EXAMPLE WITH SPECIES TRAITS
# Read in PREDICTS database extract
data(diversity)

# Correct for sampling effort
diversity <- CorrectSamplingEffort(diversity)

# Merge sites that are very close to one another and in the same land-use type
diversity <- MergeSites(diversity)

# Read beetle body length data
beetle.mass<-EdgarBeetleLength(path=paste(path.package("yarg"),"/data/beetle.mass.csv.bz2",sep=""))

# Add trait data to PREDICTS database extract
diversity<-AddTraits(diversity,beetle.mass)

# Rename trait column
diversity$log.length<-diversity$Length_derived_volume_3log10_mm

# Calculate site-level diversity metrics
sites<-SiteMetrics(diversity=diversity,
                   extra.cols=c("SSB","SSBS","Biome"),
                   traits="log.length",centralPercentile=0.95,
                   sites.are.unique=TRUE)
                   
}
