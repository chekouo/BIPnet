\name{plotVarSel}
\alias{plotVarSel}
\alias{plotVarSelGlobal}
\title{
	Plot the posterior probabilities of selecting the features
}
\description{
	\code{plotVarSel} plots, for one data type, the marginal posterior probability that each feature is selected in each component (\code{VarSelMean} from \code{\link{BIP}}). \code{plotVarSelGlobal} plots, for each data type, the probability that each feature is selected in at least one component (\code{VarSelMeanGlobal}).
}
\usage{
plotVarSel(Result, which = 1, components = NULL, type = c("index", "heatmap"),
           featureNames = NULL, threshold = 0.5, labelTop = 5, maxFeatures = 50,
           dataName = NULL, main = NULL)
plotVarSelGlobal(Result, which = NULL, featureNames = NULL, threshold = 0.5,
                 labelTop = 5, dataNames = NULL, main = NULL)
}
\arguments{
\item{Result}{ Output of \code{\link{BIP}}.}
\item{which}{ For \code{plotVarSel}, the number of one data type (not the response), in the order of \code{dataList}. For \code{plotVarSelGlobal}, the numbers of the data types to plot; by default all omics data types (IndicVar = 0).}
\item{components}{ Components to plot; by default all of them.}
\item{type}{ "index" (default) draws one panel per component, with a thin line per feature from 0 to its probability, blue at or above \code{threshold} and gray below; the panel heading gives the probability of selecting the component and the number of features at or above the threshold. "heatmap" draws features by components on a light-to-dark blue scale, for the features at or above \code{threshold} in at least one of the plotted components (at most \code{maxFeatures}, those with the highest probabilities); the probability of selecting each component is shown under its label.}
\item{featureNames}{ Names of the features. By default the column names of the data sets given to \code{BIP} (\code{Result$FeatureNames}), otherwise "Feature 1", "Feature 2", ... For \code{plotVarSelGlobal} with several data types, a list with one vector of names per plotted data type.}
\item{threshold}{ Selection threshold, drawn as a dashed line (default 0.5); \code{NULL} for no threshold.}
\item{labelTop}{ Number of features with the highest probabilities to label by name in each panel (default 5); features that would sit too close to a label already placed are skipped. 0 for no labels.}
\item{maxFeatures}{ Maximum number of features shown by the heatmap (default 50).}
\item{dataName, dataNames}{ Name of the data type (\code{plotVarSel}), or names of all data types in the order of \code{dataList} (\code{plotVarSelGlobal}).}
\item{main}{ Overall title; \code{""} for none.}
}
\details{
The components are identified only up to their order and signs. If the labels switched during the MCMC (\code{Result$RelabelInfo}) and \code{BIP} was run with \code{relabel = FALSE}, \code{plotVarSel} warns, since the probabilities of different components are then mixed; run \code{BIP} with \code{relabel = TRUE} to align them. The global probabilities of \code{plotVarSelGlobal} do not depend on the labels.
}
\value{
Invisibly, the plotted probabilities: a matrix of features by components for \code{plotVarSel}, and a list with one named vector per data type for \code{plotVarSelGlobal}.
}
\seealso{\code{\link{BIP}}, \code{\link{plotCompoSel}}}
\examples{
dat=Simulate(scenario=1,setting=1,overlap="yes",seed=1);
BA=BIP(dataList=list(dat$X1,dat$X2,dat$Y),IndicVar=c(0,0,1),groupList=list(dat$Group1,dat$Group2),
       Method="BIPnet",nbrcomp=4,sample=1000,burnin=500,relabel=TRUE);
plotVarSel(BA, which=1, dataName="Omics 1")
plotVarSel(BA, which=1, dataName="Omics 1", type="heatmap", maxFeatures=30)
plotVarSelGlobal(BA, dataNames=c("Omics 1","Omics 2","Outcome"))
}
