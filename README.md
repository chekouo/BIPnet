# BIPnet

**Bayesian integrative analysis and prediction with network information**

BIPnet is an R package for the joint analysis of several omics data types (for example gene
expression and metabolomics measured on the same samples), optional clinical covariates, and a
clinical outcome. It finds a small number of shared latent components that link the data types
to each other and to the outcome, and it selects, at the same time:

- the **components** that matter for each data type and for the outcome,
- the **features** (genes, metabolites, ...) that load on each component, and
- the **groups of features** (pathways or networks) that drive this selection, when grouping
  information is available (method `"BIPnet"`; method `"BIP"` ignores it).

The model is fitted by Markov chain Monte Carlo (MCMC), implemented in C. The fitted model can
then predict the outcome for new samples, using either one model or Bayesian model averaging.

Reference: Thierry Chekouo and Sandra Safo (2023). Bayesian integrative analysis and prediction
with application to atherosclerosis cardiovascular disease. *Biostatistics*, 24(1).

## Installation

BIPnet needs the **GSL** (GNU Scientific Library) C library.

| System | Install GSL with |
|---|---|
| macOS (Homebrew) | `brew install gsl` |
| Debian / Ubuntu | `sudo apt-get install libgsl-dev` |
| Fedora / RHEL | `sudo dnf install gsl-devel` |
| Windows | recent versions of Rtools include GSL |

Then install the package from GitHub:

```r
install.packages("remotes")
remotes::install_github("chekouo/BIPnet")
```

## Quick start

```r
library(BIPnet)

## Simulated data: two omics data types (500 features each), their grouping information
## (11 groups each) and a continuous outcome, for 200 samples
dat <- Simulate(scenario = 1, setting = 1, overlap = "yes", seed = 1)

## Fit the model. IndicVar gives the type of each data set: 0 = omics, 1 = outcome,
## 2 = clinical covariates. 'sample' is the number of MCMC draws kept after 'burnin'.
fit <- BIP(dataList  = list(dat$X1, dat$X2, dat$Y),
           IndicVar  = c(0, 0, 1),
           groupList = list(dat$Group1, dat$Group2),
           Method    = "BIPnet",
           nbrcomp   = 4,
           sample    = 5000, burnin = 1000,
           relabel   = TRUE)

## Which components and features were selected?
fit$CompoSelMean                 # probability of selecting each component, per data type
head(fit$VarSelMeanGlobal[[1]])  # probability of selecting each feature of data type 1

## Variable selection performance against the truth (simulated data only)
ComputeVarCriteria(fit$VarSelMeanGlobal[[1]], dat$TrueVar1)

## Predict the outcome of new samples (same true model, new samples)
test <- Simulate(scenario = 1, setting = 1, overlap = "yes", seed = 2)
pred <- BIPpredict(dataListNew = list(test$X1, test$X2), Result = fit, meth = "BMA")
mean((pred$ypredict - test$Y)^2) # test mean squared error
```

## Main functions

| Function | Purpose |
|---|---|
| `BIP()` | Fits the model by MCMC. |
| `BIPpredict()` | Predicts the outcome and the shared components for new samples, with Bayesian model averaging (`meth = "BMA"`) or the model of the selected features (`meth = "noBMA"`). |
| `ComputeVarCriteria()` | False positive and negative rates, F1 and AUC of a variable selection, given the true features. |
| `Simulate()` | Simulates data from the scenarios of the paper. |
| `plotCompoSel()` | Probability of selecting each component, per data type (bar charts or heatmap). |
| `plotVarSel()` | Probability of selecting each feature within each component, for one data type. |
| `plotVarSelGlobal()` | Probability of selecting each feature in at least one component, per data type. |
| `plotLoadings()` | Estimated loadings, from the model of the selected features or averaged over models. |
| `plotComponents()` | Samples on two shared components, coloured by their outcome. |

See the help pages (for example `?BIP`) for all arguments and outputs.

## Inputs

- `dataList`: a list of matrices (or data frames) with the same samples in the rows, one per
  data type, including the outcome as a single column. Omics data and covariates are
  standardised by `BIP()`; the outcome is not. Missing values and constant columns are not
  allowed.
- `IndicVar`: the type of each data set: `0` for omics data (with variable selection), `1` for
  the outcome (exactly one), `2` for clinical covariates (always kept).
- `groupList`: for `Method = "BIPnet"`, one 0/1 matrix per omics data type, with one row per
  feature and one column per group (a feature can belong to several groups).
- `nbrcomp`: number of shared components. With `nbrcomp = NULL` it is chosen from the singular
  values of the omics data.
- Priors: `probvarsel` (prior probability of selecting a feature, per data type) and the
  hyperparameters `priorcompselv`, `priorcompselo`, `priorb0`, `priorb` and `priorgrpsel`.

## Main outputs of `BIP()`

| Output | Content |
|---|---|
| `CompoSelMean` | Posterior probability of selecting each component, per data type. |
| `VarSelMean` | Posterior probability of selecting each feature within each component. |
| `VarSelMeanGlobal` | Posterior probability of selecting each feature in at least one component. |
| `GrpSelMean`, `GrpEffectMean` | Posterior probability of selecting each group, and mean group effect. |
| `EstU` | Posterior mean of the shared components (samples by components). |
| `EstLoad`, `EstLoadModel`, `PostGam` | Loadings of the model of the selected features; loadings and posterior probabilities of the top models, used for model averaging. |
| `EstSig2`, `EstIntcp` | Residual variances of the features, and intercept of the outcome. |
| `RelabelInfo` | Label switching after burn-in (see below). |

## Interpreting the components

The shared components are identified only up to their order and their signs: swapping two
components, or changing the sign of a component together with all its loadings, gives the same
model. If the component labels change during the MCMC, the component-level summaries would mix
different components. With `relabel = TRUE`, every MCMC draw after burn-in is matched to the
running mean of the components before the summaries are computed, so that each component keeps
the same meaning throughout. `fit$RelabelInfo` reports how many label swaps and sign flips were
found, whether or not `relabel` is used. Chains run separately (different `chainNbr`) are not
aligned with each other.

The sign of a component is arbitrary. In `plotLoadings()` and `plotComponents()`,
`orient = TRUE` orients each component so that its loading on the outcome is positive.

## Plots

```r
plotCompoSel(fit, dataNames = c("Omics 1", "Omics 2", "Outcome"))
plotVarSel(fit, which = 1, dataName = "Omics 1")
plotVarSelGlobal(fit, dataNames = c("Omics 1", "Omics 2", "Outcome"))
plotLoadings(fit, which = 1, dataName = "Omics 1", type = "heatmap")
plotComponents(fit, components = c(2, 3), scale = "diverging", orient = TRUE)
```

Feature names in the plots come from the column names of the data given to `BIP()`.

## License and contact

License: GPL (>= 2). For questions, please contact Thierry Chekouo (tchekouo@umn.edu).
