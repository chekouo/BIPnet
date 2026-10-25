BIP <- function(dataList = dataList, IndicVar = IndicVar, groupList = NULL, Method = Method, nbrcomp = 4, sample = 5000, burnin = 1000, nbrmaxmodels = 50,
                priorcompselv = c(1, 1), priorcompselo = c(1, 1), priorb0 = c(1, 100), priorb = c(1, 1), priorgrpsel = c(1, 1), probvarsel = rep(0.05, length(dataList)), chainNbr = 1, relabel = FALSE) {
    if (sample <= burnin) {
        stop("Argument burnin must be smaller than sample: the number of MCMC draws kept after burn-in.")
    }
    if (sample <= 20) {
        stop("Please specify a larger number of MCMC draws kept after burn-in (sample > 20)")
    }

    ## Input checks: the C code reads these inputs by position and length, so a mismatch would
    ## otherwise give wrong results without an error
    if (!is.list(dataList) || is.data.frame(dataList)) {
        stop("dataList must be a list of data sets (matrices or data frames).")
    }
    Np <- length(dataList)
    if (length(IndicVar) != Np) {
        stop(sprintf("IndicVar must have one entry per data set in dataList (%d), not %d.", Np, length(IndicVar)))
    }
    if (!all(IndicVar %in% c(0, 1, 2))) {
        stop("IndicVar entries must be 0 (omics data), 1 (response) or 2 (clinical covariates).")
    }
    if (sum(IndicVar == 1) != 1) {
        stop("IndicVar must contain exactly one response (1).")
    }
    dataList <- lapply(dataList, as.matrix) # a vector becomes a single column
    FeatureNames <- lapply(dataList, colnames) # used to label the plots; NULL when a data set has no column names
    Yresp <- as.vector(dataList[[which(IndicVar == 1)]]) # the response as given (it is not standardised), for plotComponents()
    nrows <- sapply(dataList, nrow)
    if (any(nrows != nrows[1])) {
        stop(sprintf("All data sets in dataList must have the same number of rows (samples); they have %s.", paste(nrows, collapse = ", ")))
    }
    if (ncol(dataList[[which(IndicVar == 1)]]) != 1) {
        stop("The response must be a single column.")
    }
    for (i in seq_len(Np)) {
        if (!is.numeric(dataList[[i]]) || any(!is.finite(dataList[[i]]))) {
            stop(sprintf("Data set %d of dataList contains missing, infinite or non-numeric values.", i))
        }
        if (IndicVar[i] != 1) {
            sds <- apply(dataList[[i]], 2, sd)
            if (any(sds == 0)) {
                stop(sprintf(
                    "Data set %d of dataList has constant columns, which cannot be standardised: %s.", i,
                    paste(c(head(which(sds == 0), 10), if (sum(sds == 0) > 10) "..."), collapse = ", ")
                ))
            }
        }
    }
    if (length(probvarsel) != Np) {
        stop(sprintf("probvarsel must have one entry per data set in dataList (%d), not %d.", Np, length(probvarsel)))
    }
    priors <- list(priorcompselv = priorcompselv, priorcompselo = priorcompselo, priorb0 = priorb0, priorb = priorb, priorgrpsel = priorgrpsel)
    for (nm in names(priors)) {
        if (!is.numeric(priors[[nm]]) || length(priors[[nm]]) != 2 || any(priors[[nm]] <= 0)) {
            stop(sprintf("%s must be a vector of two positive numbers.", nm))
        }
    }

    if (is.null(nbrcomp)) {
        D <- which(IndicVar == 0)
        mysvd <- lapply(D, function(i) svd(dataList[[i]]))
        mysumsvd <- lapply(1:length(D), function(i) cumsum(mysvd[[i]]$d) / max(cumsum(mysvd[[i]]$d)) * 100)
        KMax <- max(unlist(lapply(1:length(D), function(i) min(which(mysumsvd[[i]] >= 80, arr.ind = TRUE))))) # chooses maximum from D cumulative proportions
        nbrcomp <- min(KMax + 1, 10)
    }


    # Method="GroupInfo"
    if (Method == "BIP") {
        meth <- 0
    } else if (Method == "BIPnet") {
        meth <- 1
    } else {
        stop("You must provide either BIP or BIPnet")
    }

    Np <- length(dataList)
    P <- NULL
    n <- nrow(dataList[[1]])
    P <- NULL
    MeanData <- list()
    SD <- list()
    for (i in 1:Np) {
        dataList[[i]] <- as.matrix(dataList[[i]])
        P[i] <- ncol(dataList[[i]])
        if ((IndicVar[i] == 0) || (IndicVar[i] == 2)) {
            MeanData[[i]] <- apply(dataList[[i]], 2, mean)
            SD[[i]] <- apply(dataList[[i]], 2, sd)
            dataList[[i]] <- scale(as.matrix(dataList[[i]]), T, T)
        }
        dataList[[i]] <- t(dataList[[i]])
    }
    datasets <- unlist(dataList)

    if (is.null(groupList) && (meth == 1)) {
        stop("You must provide a list of grouping information")
    } else if (meth == 1) {
        omics <- which(IndicVar == 0)
        if (length(groupList) != length(omics)) {
            stop(sprintf("groupList must contain one group matrix per omics data set (IndicVar = 0): %d, not %d.", length(omics), length(groupList)))
        }
        for (ll in seq_along(omics)) {
            G <- as.matrix(groupList[[ll]])
            if (nrow(G) != P[omics[ll]]) {
                stop(sprintf("Group matrix %d has %d rows, but data set %d has %d columns (features).", ll, nrow(G), omics[ll], P[omics[ll]]))
            }
            if (!all(G %in% c(0, 1))) {
                stop(sprintf("Group matrix %d must contain only 0 and 1.", ll))
            }
        }
    } else if (is.null(groupList)) {
        for (ll in 1:length(IndicVar)) {
            groupList[[ll]] <- matrix(1, P[ll], 1)
        }
    }


    ll <- 1
    K <- NULL
    for (i in 1:Np) {
        if (IndicVar[i] != 0) {
            K[i] <- 1
        } else {
            K[i] <- ncol(groupList[[ll]])
            ll <- ll + 1
        }
    }
    groupList <- lapply(groupList, t)
    Paths <- unlist(groupList)

    result <- .C("mainfunction",
        Method1 = as.integer(meth), n1 = as.integer(n), P = as.integer(P), r1 = as.integer(nbrcomp), Np1 = as.integer(Np), datasets = as.double(datasets), IndVar = as.integer(IndicVar), K = as.integer(K), Paths = as.integer(Paths), maxmodel1 = as.integer(nbrmaxmodels), nbrsample1 = as.integer(sample), burninsample1 = as.integer(burnin),
        CompoSelMean = as.double(rep(0, Np * nbrcomp)), VarSelMean = as.double(rep(0, nbrcomp * sum(P))),
        VarSelMeanGlobal = as.double(rep(0, sum(P))), GrpSelMean = as.double(rep(0, nbrcomp * sum(K))),
        GrpEffectMean = as.double(rep(0, nbrcomp * sum(K))), IntGrpMean = as.double(rep(0, nbrcomp * Np)),
        EstU = as.double(rep(0, n * nbrcomp)), EstSig2 = as.double(rep(0, sum(P))), InterceptMean = as.double(rep(0, 1)),
        EstLoadMod = as.double(rep(0, nbrmaxmodels * nbrcomp * sum(P))), EstLoad = as.double(rep(0, nbrcomp * sum(P))),
        nbrmodel1 = as.integer(0), postgam = rep(0, nbrmaxmodels), priorcompsel = as.double(priorcompselv),
        priorcompselo = as.double(priorcompselo), priorb0 = as.double(priorb0), priorb = as.double(priorb), priorgrpsel = as.double(priorgrpsel),
        probvarsel = as.double(probvarsel), chainNbr = as.integer(chainNbr),
        AcceptB0Rate = as.double(rep(0, Np * nbrcomp)), relabel1 = as.integer(isTRUE(relabel)), RelabelInfo = as.double(rep(0, 3))
    )


    reseffect <- result$EstLoadMod
    nbrmodel <- result$nbrmodel1
    EstLoadModel <- rep(list(list()), nbrmodel)
    for (mo in 1:nbrmodel) {
        for (m in 1:Np) {
            x <- sum(P[1:(m - 1)])
            y <- sum(P[1:m])
            if (m == 1) {
                x <- 0
            }
            init <- 1 + x * nbrcomp + (mo - 1) * nbrcomp * sum(P)
            final <- y * nbrcomp + (mo - 1) * nbrcomp * sum(P)
            EstLoadModel[[mo]][[m]] <- matrix(reseffect[init:final], nbrcomp, byrow = T)
        }
    }

    resvarsel1 <- result$VarSelMeanGlobal
    resvarsel2 <- result$VarSelMean
    resvarsel3 <- result$GrpSelMean
    resvarsel4 <- result$GrpEffectMean
    resvarsel5 <- result$EstLoad
    EstimateU <- matrix(result$EstU, n, byrow = T)
    CompoSelMean <- matrix(result$CompoSelMean, Np, byrow = T)
    IntGrpMean <- matrix(result$IntGrpMean, Np, byrow = T)
    RelabelInfo <- c(
        relabel = isTRUE(relabel), labelSwaps = result$RelabelInfo[1],
        signFlips = result$RelabelInfo[2], pctRelabelled = result$RelabelInfo[3]
    )
    AcceptB0Rate <- matrix(result$AcceptB0Rate, Np, byrow = T)
    AcceptB0Rate[AcceptB0Rate < 0] <- NA # NA: component never active after burn-in
    VarSelMeanGlobal <- list()
    VarSelMean <- list()
    GrpSelMean <- list()
    GrpEffectMean <- list()
    EstLoad <- list()
    EstSig2 <- list()
    m1 <- m2 <- m3 <- 1
    for (m in 1:Np) {
        VarSelMeanGlobal[[m]] <- resvarsel1[m1:(m1 - 1 + P[m])]
        VarSelMean[[m]] <- matrix(resvarsel2[m2:(m2 - 1 + P[m] * nbrcomp)], P[m], byrow = T)
        GrpSelMean[[m]] <- matrix(resvarsel3[m3:(m3 - 1 + K[m] * nbrcomp)], nbrcomp, byrow = T)
        GrpEffectMean[[m]] <- matrix(resvarsel4[m3:(m3 - 1 + K[m] * nbrcomp)], nbrcomp, byrow = T)
        EstLoad[[m]] <- matrix(resvarsel5[m2:(m2 - 1 + P[m] * nbrcomp)], nbrcomp, byrow = T)
        EstSig2[[m]] <- result$EstSig2[m1:(m1 - 1 + P[m])]
        m1 <- m1 + P[m]
        m2 <- m2 + P[m] * nbrcomp
        m3 <- m3 + K[m] * nbrcomp
    }

    return(list(
        EstU = EstimateU, VarSelMean = VarSelMean, VarSelMeanGlobal = VarSelMeanGlobal,
        CompoSelMean = CompoSelMean, GrpSelMean = GrpSelMean, GrpEffectMean = GrpEffectMean,
        IntGrpMean = IntGrpMean, AcceptB0Rate = AcceptB0Rate, RelabelInfo = RelabelInfo,
        EstLoad = EstLoad, EstLoadModel = EstLoadModel, nbrmodel = result$nbrmodel1,
        EstSig2 = EstSig2, EstIntcp = result$InterceptMean, PostGam = result$postgam,
        IndicVar = IndicVar, nbrcomp = nbrcomp, MeanData = MeanData, SDData = SD, FeatureNames = FeatureNames, Y = Yresp
    ))
}
