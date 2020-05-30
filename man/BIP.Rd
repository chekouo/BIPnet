\name{BIP}
\alias{BIP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	An MCMC algorithm to perform a Bayesian integrative approach for multi-omics data.  
}
\description{
The method can incorporate clinical covariates, clinical response variables and external grouping network information. It performs two methods: i) BIPnet: A Bayesian approach that measures the association between multiple data types and incorporates external prior grouping information; ii) BIP does not account for external grouping information. The algorithm computes (i) the marginal posterior probabilities (MPP) of inclusions for each component and each feature within active components; (ii) posterior mean of the share latent components; iii) Estimated Loading obtained using a threshold on the MPP of inclusions for  each feature; iv) Posterior means of group effects. 
}

\usage{
BIP(dataList=dataList,IndicVar=IndicVar, groupList=NULL,Method=Method,nbrcomp=4, 
    sample=5000, burnin=1000,nbrmaxmodels=50,priorcompselv=c(1,1),priorcompselo=c(1,1),
    priorb0=c(2,2),priorb=c(1,1),priorgrpsel=c(1,1),probvarsel=0.05)
}

\arguments{
	\item{dataList}{ List of data sets: If M is the number of omics data types, then this list isof length M+1. Otherwise, if we include the covariate set, then its of length M+2. Each omics dataset object is a data matrix  measured in a set of n samples. For each omics matrix, the rows represent samples, and the columns represent genomic features. 
	}
	\item{IndicVar}{ This outputs indicates which type of variables. For inctance, if IndicVar=c(0,0,1,2), the first two variable are omicds data that would require variable selection, the third variable is the response variable and the fourth is the set of clinical covariates.}

	\item{groupList}{ List of groups of lenght M for each omics data. For each omics data type l, group is a matrix \eqn{ P_l \times K_l} where \eqn{P_l} is the number of features, \eqn{K_l} is the number of groups and each entry is either 1 or 0, with 1 if the feature belongs to the corresponding group. If Method="BIP", then groupList=NULL, the default value. 
	}
	\item{Method}{ Method is either BIPnet (with network information) or BIP (withouth network information)}
	\item{nbrcomp}{ Number of components of the shared components}
	\item{sample }{ Total number of MCMC draws. It must be larger than burnin. }
	\item{burnin }{ 
		Number of draws to discard for burn-in
	}
	\item{nbrmaxmodels}{ Number maximum of models for prediction}
	\item{priorcompselv}{ Hyperparameters of a beta(a,b) distribution, prior distribution of the probability of selecting components}
\item{priorcompselo}{ Hyperparameters of a beta(a,b) distribution, prior distribution of the probability of selecting the response on a component.}
\item{priorb0}{ Hyperparameters of a gamma distribution; prior distribution of the paramater \eqn{b_{l0}} that controls the shrinkage of loadings when there is no grouping information }
\item{priorb}{ Hyperparameters of a gamma distribution; prior distribution of group effect coefficients \eqn{b_{lj}}, \eqn{j=1,...,K_l}.} 
\item{priorgrpsel}{ Hyperparameters of a beta(a,b) distribution, prior distribution of the probability of selecting groups.}
\item{probvarsel}{ Prior probability for variable selection}
}
\details{
	The function will return several R objects, which can be assigned to a variable. To see the results, use the ``$" operator. 
}
\value{
\item{EstU}{ Posterior mean of the shared components \eqn{U}. Important components is determined using  MPPs of \eqn{\gamma}}
\item{VarSelMean}{ MPP of feature selection for each omics data type for every component; this computes the proportion of iterations that a  variable or feature is selected in each component \eqn{l=1,...,r} }
\item{VarSelMeanGlobal}{ MPP of feature selection for each omics data type; this computes the proportion of iterations that a  variable or feature is selected in at least one component.}
\item{CompoSelMean}{ MPPs of the selection of components \eqn{l}.}
\item{GrpSelMean}{ MPPs of the selection of groups if BIPnet is chosen.}
\item{GrpEffectMean}{ Posterior mean of the group effects}
\item{....}{ The other parameters are used for prediction}
}

\references{
 Author names (2020), \emph{Bayesian Integrative Analysis and Prediction with Application to Atherosclerosis Cardiovascular Disease}, \emph{submitted}.
}

\seealso{
\code{\link{BIPpredict}}
}

\examples{
### We run Scenario 1-Setting 1 in the manuscript
library(BIPnet);
dat=Simulate(scenario=1,setting=1,overlap="yes",seed=1);

dataList=list(dat$X1,dat$X2,dat$Y)
PathList=list(dat$Group1,dat$Group2)
BA=BIP(dataList=dataList,IndicVar=c(0,0,1),groupList=PathList,Method="BIPnet",nbrcomp=4, 
       sample=5000, burnin=1000);
TrueVar1=dat$TrueVar1;
IndicVar=BA$IndicVar;
Np=length(IndicVar)
for (m in 1:Np){
if (IndicVar[m]==0){
print(ComputeVarCriteria(BA$VarSelMeanGlobal[[m]],TrueVar1))
}
}

}

