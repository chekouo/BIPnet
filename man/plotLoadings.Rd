\name{plotLoadings}
\alias{plotLoadings}
\title{
	Plot the estimated loadings
}
\description{
	Plots, for one data type, the estimated loading of each feature on each component, either from the model of the selected features (\code{EstLoad} from \code{\link{BIP}}) or averaged over the top models (Bayesian model averaging of \code{EstLoadModel} with weights \code{PostGam}).
}
\usage{
plotLoadings(Result, which = 1, components = NULL, type = c("index", "heatmap"),
             model = c("selected", "BMA"), featureNames = NULL, labelTop = 5,
             maxFeatures = 50, orient = FALSE, dataName = NULL, main = NULL)
}
\arguments{
\item{Result}{ Output of \code{\link{BIP}}.}
\item{which}{ Number of one data type (not the response), in the order of \code{dataList}.}
\item{components}{ Components to plot; by default all of them.}
\item{type}{ "index" (default) draws one panel per component, with a thin line per feature from 0 to its loading, blue if positive and red if negative, on a y-axis shared by all panels; the panel heading gives the probability of selecting the component, its loading on the outcome and the number of non-zero loadings. "heatmap" draws features by components on a diverging red (negative) - gray (0) - blue (positive) scale with a colour key, for the features with a non-zero loading in at least one of the plotted components (at most \code{maxFeatures}, those with the largest absolute loadings).}
\item{model}{ "selected" (default): the loadings of the model made of the components and features with posterior probabilities of at least 0.5 (\code{Result$EstLoad}); the loadings of the other features are 0. "BMA": the average of the loadings of the top models (\code{Result$EstLoadModel}), weighted by their posterior probabilities (\code{Result$PostGam}).}
\item{featureNames}{ Names of the features. By default the column names of the data set given to \code{BIP} (\code{Result$FeatureNames}), otherwise "Feature 1", "Feature 2", ...}
\item{labelTop}{ Number of features with the largest absolute loadings to label by name in each panel of the index plot (default 5); features whose label would overlap one already placed are skipped. 0 for no labels.}
\item{maxFeatures}{ Maximum number of features shown by the heatmap (default 50).}
\item{orient}{ The sign of a component is arbitrary: changing the signs of its shared component and of all its loadings gives the same model. If TRUE, the sign of each component is changed for the plot so that its loading on the outcome is positive; a positive loading then means that the feature increases with the outcome through that component. Default FALSE (the signs of the fit).}
\item{dataName}{ Name of the data type, used in the titles.}
\item{main}{ Overall title; \code{""} for none.}
}
\details{
The data are standardised by \code{BIP}, so the loadings of the features of a data type are on a common scale. The loadings are point estimates; no posterior intervals are available.

The components are identified only up to their order and signs. If the labels switched during the MCMC (\code{Result$RelabelInfo}) and \code{BIP} was run with \code{relabel = FALSE}, the function warns, since the loadings of different components, or of the same component with opposite signs, are then mixed; run \code{BIP} with \code{relabel = TRUE} to align them.
}
\value{
Invisibly, the plotted loadings as a matrix of features by components (for the heatmap, only the features shown).
}
\seealso{\code{\link{BIP}}, \code{\link{plotVarSel}}, \code{\link{plotCompoSel}}}
\examples{
dat=Simulate(scenario=1,setting=1,overlap="yes",seed=1);
BA=BIP(dataList=list(dat$X1,dat$X2,dat$Y),IndicVar=c(0,0,1),groupList=list(dat$Group1,dat$Group2),
       Method="BIPnet",nbrcomp=4,sample=1000,burnin=500,relabel=TRUE);
plotLoadings(BA, which=1, dataName="Omics 1")
plotLoadings(BA, which=1, dataName="Omics 1", model="BMA", orient=TRUE)
plotLoadings(BA, which=2, dataName="Omics 2", type="heatmap", maxFeatures=30)
}
