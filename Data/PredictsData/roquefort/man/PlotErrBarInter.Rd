
\name{PlotErrBarInter}
\alias{PlotErrBarInter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Make an error-bar plot showing the modelled interaction between two factors}
\description{
Plots modelled effects on a specified response variable of the interaction between two factors.
}
\usage{
PlotErrBarInter<-function(model,data,responseVar,seMultiplier=1.96,outDir=NULL,
logLink="n",catInteraction=character(0),forPaper=FALSE,secdAge=FALSE,
xtext.srt=0,ylim=NA,rescale=NULL,errbar.cols=NULL,pt.pch=NULL,
params=list(),plotLabels=TRUE,plotLandUses=TRUE,cex.txt=NULL,pt.cex=1,
main=NULL,type="percent",rescaleRefLevel=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The model to derive parameter estimates from
}
  \item{data}{
A data frame containing all variables in the model
}
  \item{responseVar}{
The name of the response variable fitted, as text, to display on the y-axis of the graph
}
  \item{seMultiplier}{
The number of standard errors to use for the error bars. Default is 1.96, which shows 95 percent confidence intervals
}
  \item{outDir}{
Optional: directory to write the graph to. If not specified, then the graph will be written to the native graphics device
}
  \item{logLink}{
The link function used by the model: one of "n" (identity), "e" (log, base e), "10" (log, base 10) or "b" (logit - used in binomial models)
}
  \item{catInteraction}{
A vector of length two containing the two factors whose interaction should be plotted  
}
  \item{forPaper}{
Whether the graph is for a published paper. Default is FALSE. If TRUE, the coloured shading is removed and the legend reduced in size
}
  \item{secdAge}{
If the model contains PREDICTS land-use classes, whether secondary vegetation was also considered
}
  \item{xtext.srt}{
The number of degrees by which to rotate the x-axis labels
}
  \item{ylim}{
Fixed y-axis limits. If not specified, then resolved automatically
}
  \item{rescale}{
The amount by which to rescale y-axis values for the categorical effects by. If not specified, then no rescaling is performed
}
  \item{errbar.cols}{
Colours to overwrite the default colours with
}
  \item{pt.pch}{
The point types (using R's 'pch' values), with which to overwrite the default
}
  \item{params}{
Any R graphical parameters, which will overwrite the defaults  
}
  \item{add}{
Whether to add the plot to an existing plot; default is FALSE
}
  \item{offset}{
The amount by which to offset points on the x axis; default is zero  
}
  \item{plotLabels}{
Whether to plot individual labels for each factor level  
}
  \item{plotLandUses}{
Whether to add the main land-use labels to the x-axis  
}
  \item{cex.txt}{
Text size for additional labels (i.e. text not on the axes or in the axis labels)  
}
  \item{pt.cex}{
Point size to use. Default is 1.  
}
  \item{main}{
The title for the plot  
}
  \item{type}{
Whether to plot as a percentage change ('percent' - the default) or directly as the response variable ('response')  
}
  \item{rescaleRefLevel}{
Whether to rescale the y values so that all points for the reference level of the first factor are plotted at the same y-value 
}
}
\details{
Code for calculating predicted values and confidence intervals for continuous effects was taken from the GLMM wiki (see references).
}
\references{
http://glmm.wikidot.com/faq
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
