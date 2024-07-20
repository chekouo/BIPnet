\name{Simulate}
\alias{Simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation of simulated data as explained in the reference manuscript.	      
}
\description{ 
	This function generates data described in the manuscript.
}

\usage{
Simulate(scenario=1,setting=NULL,overlap=NULL,seed=1)
}

\arguments{
	        \item{scenario}{ We have three scenarios as decrisbed in the paper. Scenario 1 has 5 settings, and scenarios 2 and 3 have both two choices: overlap and no overlap.}
		\item{setting}{ Settings for secnario 1 only. It's a total of 4 settings}
		\item{overlap}{ If yes, then there are some overlap between components, otherwise, there is no overlap. This input only works when scenario is either 2 or 3.}
		\item{seed}{ Seed to generate random numbers.}

}

\details{
The function will return two data matrix X1 and X2, an estimated response Y, shared component matrix U, loadings A1 and A2, and the grouping information (Group1 and Group2)
}
\references{
 Thierry Chekouo and Sandra Safo (2023), Bayesian Integrative Analysis and Prediction with Application to Atherosclerosis Cardiovascular Disease, \emph{Biostatistics}, , Volume 24, Issue 1.
}

\examples{
dat=Simulate(scenario=3,setting=NULL,overlap="no",seed=1)
str(dat)
}
