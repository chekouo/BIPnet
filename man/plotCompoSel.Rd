\name{plotCompoSel}
\alias{plotCompoSel}
\title{
	Plot the posterior probabilities of selecting the components
}
\description{
	Plots, for each data type, the marginal posterior probability that each component is selected (\code{CompoSelMean} from \code{\link{BIP}}), as one bar chart per data type or as a heatmap of data types by components.
}
\usage{
plotCompoSel(Result, dataNames = NULL, which = NULL, type = c("bar", "heatmap"),
             threshold = 0.5, showValues = TRUE, main = NULL)
}
\arguments{
\item{Result}{ Output of \code{\link{BIP}}.}
\item{dataNames}{ Names of the data types, in the order of \code{dataList} in \code{BIP}. By default "Data type 1", "Data type 2", ..., with "Response" for the response and "Covariates" for clinical covariates.}
\item{which}{ Numbers of the data types to plot. By default all of them except the clinical covariates, whose components are always selected.}
\item{type}{ "bar" (default) draws one bar chart per data type, with a dashed line at \code{threshold}; "heatmap" draws a single grid of data types by components on a light-to-dark blue scale.}
\item{threshold}{ Selection threshold drawn as a dashed line on the bar charts (default 0.5, as used by \code{BIP} for the final model); \code{NULL} for no line.}
\item{showValues}{ If TRUE (default), the probabilities are printed on the bars or in the cells.}
\item{main}{ Overall title; \code{""} for none.}
}
\details{
The components are identified only up to their order and signs. If the labels switched during the MCMC (\code{Result$RelabelInfo}) and \code{BIP} was run with \code{relabel = FALSE}, the function warns, since the probabilities of different components are then mixed; run \code{BIP} with \code{relabel = TRUE} to align them.
}
\value{
Invisibly, the matrix of plotted probabilities (data types by components).
}
\seealso{\code{\link{BIP}}}
\examples{
dat=Simulate(scenario=1,setting=1,overlap="yes",seed=1);
BA=BIP(dataList=list(dat$X1,dat$X2,dat$Y),IndicVar=c(0,0,1),groupList=list(dat$Group1,dat$Group2),
       Method="BIPnet",nbrcomp=4,sample=1000,burnin=500,relabel=TRUE);
plotCompoSel(BA, dataNames=c("Omics 1","Omics 2","Outcome"))
plotCompoSel(BA, dataNames=c("Omics 1","Omics 2","Outcome"), type="heatmap")
}
