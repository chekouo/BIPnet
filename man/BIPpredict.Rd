\name{BIPpredict}
\alias{BIPpredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Prediction of the response variable and shared components
}
\description{
	This function computes  predicted response and shared components on a test data set. 
}

\usage{
BIPpredict(dataListNew=dataListNew,Result=Result,meth="BMA")
}

\arguments{
\item{dataListNew}{ List of new data sets: Do not include the response as you may not have be available for a new set of samples. 
			        }
\item{Result}{ This is a BIP object obtained using function BIP on the training set.}
 \item{meth }{ meth is either "BMA" (Bayesian model averaging) or "noBMA" for no BMA. In this case, the prediction is performed using only one model obtained using  threshold of 0.5 on the MPPs}

}

\details{
The function will return both predicted responses and predicted shared latent components on a test data.
}
\references{
	 Author names (2020), \emph{Bayesian Integrative Analysis and Prediction with Application to Atherosclerosis Cardiovascular Disease}, \emph{submitted}.
}
\examples{
### We run Scenario 1-Setting 1 in the manuscript
library(BIPnet);
dat=Simulate(scenario=1,setting=1,overlap="yes",seed=1);
dataList=list(dat$X1,dat$X2,dat$Y)
PathList=list(dat$Group1,dat$Group2)
BA=BIP(dataList=dataList,IndicVar=c(0,0,1),groupList=PathList,Method="BIPnet",nbrcomp=4, 
       sample=5000, burnin=1000);
#### Prediction performance
datTest=Simulate(scenario=1,setting=1,overlap="yes",seed=2);
X_new=list(datTest$X1,datTest$X2);
Pred=BIPpredict(dataListNew=X_new,Result=BA,meth="BMA")
### Compute mse

y_new=datTest$Y;
mse=mean((Pred$ypredict-as.vector(y_new))^2);
mse
}


