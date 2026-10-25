\name{plotComponents}
\alias{plotComponents}
\title{
	Scatter plot of two shared components, coloured by the outcome
}
\description{
	Plots the samples on two estimated shared components (by default the posterior means \code{EstU} from \code{\link{BIP}}), with each sample coloured by its value of the outcome. An arrow shows the direction in which the outcome increases according to the model.
}
\usage{
plotComponents(Result, components = c(1, 2), y = NULL, U = NULL,
               scale = c("sequential", "diverging"), orient = FALSE, arrow = TRUE,
               outcomeName = "Outcome", main = NULL)
}
\arguments{
\item{Result}{ Output of \code{\link{BIP}}.}
\item{components}{ The two components to plot: the first on the horizontal axis, the second on the vertical axis.}
\item{y}{ Outcome values used to colour the samples, one per row of \code{U}. By default the response used to fit the model (\code{Result$Y}).}
\item{U}{ Component scores to plot, a matrix with one row per sample and one column per component. By default the posterior means of the training samples (\code{Result$EstU}); for instance \code{BIPpredict(...)$Upredtest} with the outcome of the test samples in \code{y}.}
\item{scale}{ "sequential" (default): a single blue hue, from light (lowest value) to dark (highest value). "diverging": red below the mean of \code{y}, gray at the mean and blue above it.}
\item{orient}{ The sign of a component is arbitrary: changing the signs of a shared component and of all its loadings gives the same model. If TRUE, the sign of each component is changed for the plot so that its loading on the outcome is positive, as in \code{\link{plotLoadings}}. Default FALSE (the signs of the fit).}
\item{arrow}{ If TRUE (default), draws an arrow from the origin in the direction of the loadings of the outcome on the two components (from \code{Result$EstLoad}): the outcome predicted by the model increases in that direction.}
\item{outcomeName}{ Name of the outcome, used for the colour key and the titles.}
\item{main}{ Overall title; \code{""} for none.}
}
\details{
Each axis label gives the loading of the component on the outcome and the posterior probability that the component is linked to the outcome (the row of \code{CompoSelMean} for the response). The two axes have the same scale.

\code{EstU} contains posterior means, which are shrunk towards 0 compared with single MCMC draws. The components are identified only up to their order and signs; if the labels switched during the MCMC (\code{Result$RelabelInfo}) and \code{BIP} was run with \code{relabel = FALSE}, the function warns. When several true components load on the same features, individual components are hard to interpret on their own, but the arrow and the colour gradient still show how the outcome varies in the plane of the two components.
}
\value{
Invisibly, a data frame with the plotted coordinates of the two components and the outcome values.
}
\seealso{\code{\link{BIP}}, \code{\link{BIPpredict}}, \code{\link{plotLoadings}}}
\examples{
dat=Simulate(scenario=1,setting=1,overlap="yes",seed=1);
BA=BIP(dataList=list(dat$X1,dat$X2,dat$Y),IndicVar=c(0,0,1),groupList=list(dat$Group1,dat$Group2),
       Method="BIPnet",nbrcomp=4,sample=1000,burnin=500,relabel=TRUE);
plotComponents(BA, components=c(2,3))
plotComponents(BA, components=c(2,3), scale="diverging", orient=TRUE)
#### Test samples, with their predicted components
datTest=Simulate(scenario=1,setting=1,overlap="yes",seed=2);
Pred=BIPpredict(dataListNew=list(datTest$X1,datTest$X2),Result=BA)
plotComponents(BA, components=c(2,3), U=Pred$Upredtest, y=datTest$Y)
}
