
\name{PlotErrBar}
\alias{PlotErrBar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Make an error-bar plot showing model coefficients}
\description{
Plots modelled effects on a specified response variable of specified factors and continuous variables. Continuous effects are plotted as fitted values at high, medium and low values of the continuous variable.
}
\usage{
PlotErrBar(model,data,responseVar,seMultiplier=1.96,outDir=NULL,
logLink="n",catEffects=NULL,contEffects=list(),contEffectsLabels=NULL,
otherCatEffects=list(),otherContEffects=character(0),forPaper=FALSE,
align=FALSE,secdAge=FALSE,xtext.srt=0,ylim=NA,order=NULL,rescale=NULL,
errbar.cols=NULL,pt.pch=NULL,params=list(),add=FALSE,offset=0,
plotLabels=TRUE,plotLandUses=TRUE,cex.txt=NULL,pt.cex=1,main=NULL,
modelNoInteractions=NULL,type="percent",matchLandUses=TRUE,
returnColours=FALSE)
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
  \item{catEffects}{
The categorical effects (factors) to plot, as a vector  
}
  \item{contEffects}{
The continuous effects to plot, as a list where the names in the list are the effects and the values are 1 (plot values of the coninutous effects from low to high), 2 (plot values of the continuous effects from high to low) or 0 (don't plot the values, but create a space in the plot - useful if you want to multiple-panel plots to align)
}
  \item{contEffectsLabels}{
The text labels to use for the continuous effects
}
  \item{otherCatEffects}{
Other categorical effects that appear in the model, but which you don't want to plot, given as a list where the names indicate the name of the factor and the values indicate the levels in the factor that you want to use as the reference level for plotting
}
  \item{otherContEffects}{
Other continuous effects that appear in the model, but which you don't want to plot, given as a character vector of the names of the terms. Only necessary if some continuous effects are plotted.  
}
  \item{forPaper}{
Whether the graph is for a published paper. Default is FALSE. If TRUE, the coloured shading is removed and the legend reduced in size
}
  \item{align}{
Whether to align recognized PREDICTS land-use classes for consistency across plots
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
  \item{order}{
The order to plot the categorical factor levels in. If not specified, then resolved automatically
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
  \item{modelNoInteractions}{
A model with no interaction terms, to be used for plotting continuous effects  
}
  \item{type}{
Whether to plot as a percentage change ('percent' - the default) or directly as the response variable ('response')  
}
  \item{matchLandUses}{
Whether to try to match and order factor levels according to the PREDICTS land uses  
}
  \item{returnColours}{
Whether to return the colours used for plotting land uses - useful for adding additional points  
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
# Load example PREDICTS site-level data for the UK
data(sites)
sites$UseIntensity<-sites$Use_intensity
sites$UI<-factor(paste(sites$LandUse,sites$Use_intensity))
sites$UI<-relevel(sites$UI,ref="Natural Minimal use")

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
    
# Log estimates of human population density
sites$LogHPD<-log(sites$HPD+1)

# Example plot of a model of species richness against land use
m1<-FixedModel(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",fixedStruct="LandUse",
randomStruct="(1|SS)+(1|SSB)")

PlotErrBar(model=m1$model,data=m1$data,responseVar="Species richness",
logLink="e",catEffects="LandUse",align=FALSE)

# To order the plot so that natural vegetation comes first
PlotErrBar(model=m1$model,data=m1$data,responseVar="Species richness",
logLink="e",catEffects="LandUse",align=FALSE,order=c(4,1,2,3))

# Example plot with a continuous effect included
m2<-FixedModel(all.data=sites,responseVar="Species_richness",
fitFamily="poisson",
fixedStruct="LandUse+poly(LogHPD,2)+LandUse:poly(LogHPD,2)",
randomStruct="(1|SS)+(1|SSB)")

PlotErrBar(model=m2$model,data=m2$data,responseVar="Species richness",
logLink="e",catEffects="LandUse",contEffects=list(LogHPD=1),
contEffectsLabels="Human popn.",align=FALSE,order=c(4,1,2,3))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
