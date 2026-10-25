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
    priorb0=c(1,100),priorb=c(1,1),priorgrpsel=c(1,1),probvarsel=rep(0.05,length(dataList)),chainNbr=1,relabel=FALSE)
}

\arguments{
	\item{dataList}{ List of data sets: If M is the number of omics data types, then this list is of length M+1. Otherwise, if we include the covariate set, then its of length M+2. Each omics dataset object is a data matrix  measured in a set of n samples. For each omics matrix, the rows represent samples, and the columns represent genomic features. 
	}
	\item{IndicVar}{ This outputs indicates which type of variables. For inctance, if IndicVar=c(0,0,1,2), the first two variable are omics data that would require variable selection, the third variable is the response variable and the fourth is the set of clinical covariates.}

	\item{groupList}{ List of groups of length M for each omics data. For each omics data type l, group is a matrix \eqn{ P_l \times K_l} where \eqn{P_l} is the number of features, \eqn{K_l} is the number of groups and each entry is either 1 or 0, with 1 if the feature belongs to the corresponding group. If Method="BIP", then groupList=NULL, the default value. 
	}
	\item{Method}{ Method is either BIPnet (with network information) or BIP (withouth network information)}
	\item{nbrcomp}{ Number of components of the shared components}
	\item{sample }{ Number of MCMC draws kept after burn-in; the posterior summaries are computed from these draws. The sampler runs sample + burnin iterations in total. It must be larger than both burnin and 20. }
	\item{burnin }{ 
		Number of initial draws discarded as burn-in, run before the sample draws that are kept.
	}
	\item{nbrmaxmodels}{ Number maximum of models for prediction}
	\item{priorcompselv}{ Hyperparameters of a beta(a,b) distribution, prior distribution of the probability of selecting components}
\item{priorcompselo}{ Hyperparameters of a beta(a,b) distribution, prior distribution of the probability of selecting a component associated with the response.}
\item{priorb0}{ Hyperparameters of a gamma distribution; prior distribution of the paramater \eqn{b_{l0}} that controls the shrinkage of loadings when there is no grouping information. The two values are the shape and the rate. The default c(1,100), with prior mean 0.01, keeps \eqn{b_{l0}} small relative to the group effects, so that the grouping information drives the shrinkage. }
\item{priorb}{ Hyperparameters of a gamma distribution; prior distribution of group effect coefficients \eqn{b_{lj}}, \eqn{j=1,...,K_l}.} 
\item{priorgrpsel}{ Hyperparameters of a beta(a,b) distribution, prior distribution of the probability of selecting groups.}
\item{probvarsel}{ Prior probability for variable selection for each data type. It's a vector of size the number of data types. Please note that when a data type m is not subject of variable selection (e.g., covariate variables or outcome), then probvarsel[m]=1.}
\item{chainNbr}{ MCMC chain number. The default is 1 for one MCMC number 1. If you want to run N multiple MCMC chains, make a loop in R as for ( i in 1:N)  BIP(...,chainNbr=i)}
\item{relabel}{ If TRUE, the components are relabelled after burn-in to correct for label switching, so that the component-level summaries (CompoSelMean, VarSelMean, GrpSelMean, GrpEffectMean, IntGrpMean, AcceptB0Rate and EstU) refer to the same component in every draw. Each draw is matched to the running mean of the aligned draws of \eqn{U} by the signed permutation of its components that maximises the absolute correlations (Hungarian algorithm); the signs are applied to \eqn{U}. The sampler itself is not changed, and the models used for prediction are stored in the aligned labels. The default FALSE keeps the labels of the sampler. Components of different chains are not aligned with each other.}
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
\item{IntGrpMean}{ Posterior mean of \eqn{b_{l0}}, an \eqn{(M+1) \times r} matrix (one row per data type, one column per component).}
\item{RelabelInfo}{ Label switching after burn-in, computed whether or not relabel is TRUE: relabel (whether the summaries were relabelled), labelSwaps (number of draws where the matching of the components changed from the previous draw), signFlips (draws where only a sign changed) and pctRelabelled (percentage of draws whose matching is not the identity).}
\item{Y}{ The response used to fit the model, as given in dataList (it is not standardised); used by \code{\link{plotComponents}} to colour the samples.}
\item{FeatureNames}{ List with the column names of each data set in dataList (NULL for a data set without column names), used to label the plots of \code{\link{plotVarSel}} and \code{\link{plotVarSelGlobal}}.}
\item{AcceptB0Rate}{ Metropolis-Hastings acceptance rate of \eqn{b_{l0}} after burn-in, same layout as IntGrpMean. It is NA for a component that was never selected after burn-in, since \eqn{b_{l0}} is then drawn from its prior.}
\item{....}{ The other parameters are used for prediction}
}

\references{
 Thierry Chekouo and Sandra Safo (2023), Bayesian Integrative Analysis and Prediction with Application to Atherosclerosis Cardiovascular Disease, \emph{Biostatistics}, , Volume 24, Issue 1.
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
       sample=5000, burnin=1000, probvarsel=rep(0.05,length(dataList)));
TrueVar1=dat$TrueVar1;
IndicVar=BA$IndicVar;
Np=length(IndicVar)
for (m in 1:Np){
if (IndicVar[m]==0){
print(ComputeVarCriteria(BA$VarSelMeanGlobal[[m]],TrueVar1))
}
}

}

