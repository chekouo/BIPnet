\name{ComputeVarCriteria}
\alias{ComputeVarCriteria}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute criteria to evaluate variable selection performamce 
}
\description{
	This function computes false negative rates (FNR), false positive rate (FPR), F1-measure (F1) and AUC. 
}

\usage{
ComputeVarCriteria(pred.prob,truegroups,thres=0.5)
}

\arguments{
  \item{pred.prob}{ (Predicted) probabilities (e.g. marginal posterior probabilities) 
			        }
\item{truegroups}{ Binary vector with 1 as positive and 0 as negative.}
\item{thres}{Threshold on the probablities used to compute FPR, FNR and F1. Default is 0.5.}

}

\details{
The function returns FPR, FNR, F1 and AUC.
}
\references{
 Thierry Chekouo and Sandra Safo (2023), Bayesian Integrative Analysis and Prediction with Application to Atherosclerosis Cardiovascular Disease, \emph{Biostatistics}, , Volume 24, Issue 1.
}
