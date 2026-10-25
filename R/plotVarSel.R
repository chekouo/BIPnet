## Plot the posterior probabilities of selecting the features: within each component
## (Result$VarSelMean, plotVarSel) or in at least one component (Result$VarSelMeanGlobal,
## plotVarSelGlobal).

## Feature names for data type m: the argument, else the column names kept by BIP(), else numbers
.featureLabels <- function(Result, m, featureNames, P) {
    if (!is.null(featureNames)) {
        if (length(featureNames) != P) {
            stop(sprintf("featureNames must have one name per feature of data type %d (%d).", m, P))
        }
        return(as.character(featureNames))
    }
    nm <- Result$FeatureNames[[m]]
    if (!is.null(nm) && length(nm) == P) {
        return(nm)
    }
    paste("Feature", seq_len(P))
}

## One panel: a thin line per feature from 0 to its probability, blue at or above the threshold
.needlePanel <- function(y, threshold, labels, labelTop, heading, subheading) {
    ink <- "#0b0b0b"
    ink2 <- "#52514e"
    P <- length(y)
    x <- seq_len(P)
    on <- if (is.null(threshold)) rep(TRUE, P) else y >= threshold
    plot(x, y, type = "n", ylim = c(0, 1.18), xlim = c(0.5, P + 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
    abline(h = c(0.25, 0.5, 0.75, 1), col = "#e4e3df", lwd = 1)
    ## segments() fails on empty coordinates, e.g. when no feature reaches the threshold
    if (any(!on)) segments(x[!on], 0, x[!on], y[!on], col = "#b8b7b2", lwd = 1, lend = 1)
    if (any(on)) segments(x[on], 0, x[on], y[on], col = "#2a78d6", lwd = if (P > 300) 1 else 2, lend = 1)
    if (!is.null(threshold)) {
        abline(h = threshold, col = ink2, lty = 2, lwd = 1)
    }
    axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), lwd = 0, lwd.ticks = 0)
    axis(1, lwd = 0, lwd.ticks = 1, col.ticks = "#b8b7b2", tcl = -0.25)
    ## label the highest features: a label is skipped if it would overlap one already placed,
    ## and labels are kept inside the plotting region
    if (labelTop > 0 && any(on)) {
        cand <- order(-y, x)
        cand <- cand[on[cand]]
        usr <- par("usr")
        half <- 0.55 * strwidth(labels, cex = 0.7) # half label width, with some padding, in feature units
        lx <- pmin(pmax(x, usr[1] + half), usr[2] - half)
        keep <- integer(0)
        for (j in cand) {
            if (length(keep) >= labelTop) break
            if (all(abs(lx[j] - lx[keep]) > half[j] + half[keep])) keep <- c(keep, j)
        }
        text(lx[keep], y[keep], labels = labels[keep], pos = 3, offset = 0.3, cex = 0.7, col = ink, xpd = NA)
    }
    mtext(heading, side = 3, line = 0.55, adj = 0, font = 2, col = ink, cex = 0.9)
    mtext(subheading, side = 3, line = 0.55, adj = 1, col = ink2, cex = 0.75)
    mtext("Probability", side = 2, line = 2.6, las = 0, col = ink2, cex = 0.75)
}

.labelWarning <- function(Result) {
    info <- Result$RelabelInfo
    if (!is.null(info) && !isTRUE(as.logical(info["relabel"])) && (info["labelSwaps"] + info["signFlips"]) > 0) {
        warning("Component labels switched during the MCMC (see Result$RelabelInfo); rerun BIP() with relabel = TRUE so that each component keeps the same meaning.")
    }
}

plotVarSel <- function(Result, which = 1, components = NULL, type = c("index", "heatmap"), featureNames = NULL,
                       threshold = 0.5, labelTop = 5, maxFeatures = 50, dataName = NULL, main = NULL) {
    type <- match.arg(type)
    if (is.null(Result$VarSelMean) || is.null(Result$IndicVar)) {
        stop("Result must be the output of BIP().")
    }
    IndicVar <- Result$IndicVar
    if (length(which) != 1 || !(which %in% seq_along(IndicVar)) || IndicVar[which] == 1) {
        stop("which must be the number of one data type that is not the response.")
    }
    probs <- Result$VarSelMean[[which]] # features x components
    P <- nrow(probs)
    r <- ncol(probs)
    if (is.null(components)) {
        components <- seq_len(r)
    }
    if (!all(components %in% seq_len(r))) {
        stop(sprintf("components must be component numbers between 1 and %d.", r))
    }
    .labelWarning(Result)
    labels <- .featureLabels(Result, which, featureNames, P)
    if (is.null(dataName)) {
        dataName <- paste("data type", which)
    }
    if (is.null(main)) {
        main <- sprintf("Posterior probability of selecting each feature of %s, by component", dataName)
    }
    compo <- Result$CompoSelMean[which, ]
    ink <- "#0b0b0b"
    ink2 <- "#52514e"
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    if (type == "index") {
        par(mfrow = c(length(components), 1), mar = c(2.2, 4.2, 2.2, 1), oma = c(1.6, 0, if (nzchar(main)) 2.5 else 0, 0),
            mgp = c(2.6, 0.5, 0), las = 1, col.axis = ink2, fg = ink2, cex.axis = 0.8)
        par(cex = 0.75)
        for (k in components) {
            y <- probs[, k]
            nsel <- if (is.null(threshold)) NA else sum(y >= threshold)
            sub <- sprintf("component selected with probability %.2f%s", compo[k],
                           if (is.na(nsel)) "" else sprintf("  |  %d feature%s at or above %.2f", nsel, if (nsel == 1) "" else "s", threshold))
            .needlePanel(y, threshold, labels, labelTop, paste0("C", k), sub)
        }
        mtext("Feature", side = 1, line = 0.4, outer = TRUE, col = ink2, cex = 0.8)
        if (nzchar(main)) {
            mtext(main, side = 3, line = 0.8, outer = TRUE, font = 2, col = ink, cex = 0.95)
        }
    } else {
        sub <- probs[, components, drop = FALSE]
        best <- apply(sub, 1, max)
        rows <- if (is.null(threshold)) seq_len(P) else which(best >= threshold)
        rows <- rows[order(-best[rows], rows)]
        if (length(rows) > maxFeatures) {
            rows <- rows[seq_len(maxFeatures)]
        }
        if (length(rows) == 0) {
            plot.new()
            text(0.5, 0.5, sprintf("No feature of %s reaches probability %.2f in these components.", dataName, threshold), col = ink2)
            return(invisible(sub[integer(0), , drop = FALSE]))
        }
        z <- sub[rows, , drop = FALSE]
        rownames(z) <- labels[rows]
        colnames(z) <- paste0("C", components)
        nr <- nrow(z)
        nc <- ncol(z)
        pal <- colorRampPalette(c("#f7fafe", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"))(101)
        cexr <- if (nr > 30) 0.55 else 0.75
        par(mar = c(3.5, max(5, max(nchar(rownames(z))) * 0.5 * cexr / 0.75 + 1), if (nzchar(main)) 3.2 else 1, 1),
            mgp = c(2, 0.4, 0), las = 1, col.axis = ink2, fg = ink2)
        image(seq_len(nc), seq_len(nr), t(z[nr:1, , drop = FALSE]), zlim = c(0, 1), col = pal, axes = FALSE, xlab = "", ylab = "")
        abline(v = seq(0.5, nc + 0.5, 1), col = "#fcfcfb", lwd = 2)
        if (nr <= 60) abline(h = seq(0.5, nr + 0.5, 1), col = "#fcfcfb", lwd = 1)
        axis(1, at = seq_len(nc), labels = sprintf("C%d\n(%.2f)", components, compo[components]), lwd = 0, lwd.ticks = 0, cex.axis = 0.8, padj = 0.5)
        axis(2, at = seq_len(nr), labels = rev(rownames(z)), lwd = 0, lwd.ticks = 0, cex.axis = cexr)
        mtext("Component (probability of selecting it)", side = 1, line = 2.4, col = ink2, cex = 0.8)
        if (nzchar(main)) {
            title(main, col.main = ink, font.main = 2, cex.main = 0.95, line = 1.8)
            mtext(sprintf("%d feature%s at or above %.2f in at least one component%s", length(rows), if (length(rows) == 1) "" else "s",
                          if (is.null(threshold)) 0 else threshold, if (length(rows) == maxFeatures) sprintf(" (first %d shown)", maxFeatures) else ""),
                  side = 3, line = 0.5, col = ink2, cex = 0.75)
        }
        return(invisible(z))
    }
    invisible(probs[, components, drop = FALSE])
}

plotVarSelGlobal <- function(Result, which = NULL, featureNames = NULL, threshold = 0.5, labelTop = 5,
                             dataNames = NULL, main = NULL) {
    if (is.null(Result$VarSelMeanGlobal) || is.null(Result$IndicVar)) {
        stop("Result must be the output of BIP().")
    }
    IndicVar <- Result$IndicVar
    Np <- length(IndicVar)
    if (is.null(which)) {
        which <- seq_len(Np)[IndicVar == 0]
    }
    if (!all(which %in% seq_len(Np)) || any(IndicVar[which] == 1)) {
        stop("which must contain the numbers of data types that are not the response.")
    }
    if (is.null(dataNames)) {
        dataNames <- paste("Data type", seq_len(Np))
    } else if (length(dataNames) != Np) {
        stop(sprintf("dataNames must have one name per data type (%d).", Np))
    }
    if (!is.null(featureNames) && !is.list(featureNames)) {
        if (length(which) != 1) stop("With several data types, featureNames must be a list with one vector of names per data type.")
        featureNames <- setNames(list(featureNames), which)
    }
    .labelWarning(Result)
    if (is.null(main)) {
        main <- "Posterior probability of selecting each feature in at least one component"
    }
    ink <- "#0b0b0b"
    ink2 <- "#52514e"
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(length(which), 1), mar = c(2.2, 4.2, 2.2, 1), oma = c(1.6, 0, if (nzchar(main)) 2.5 else 0, 0),
        mgp = c(2.6, 0.5, 0), las = 1, col.axis = ink2, fg = ink2, cex.axis = 0.8)
    par(cex = 0.75) # same text size whatever the number of panels
    out <- list()
    for (m in which) {
        y <- Result$VarSelMeanGlobal[[m]]
        fn <- if (is.null(featureNames)) NULL else if (!is.null(names(featureNames))) featureNames[[as.character(m)]] else featureNames[[match(m, which)]]
        labels <- .featureLabels(Result, m, fn, length(y))
        nsel <- if (is.null(threshold)) NA else sum(y >= threshold)
        sub <- if (is.na(nsel)) "" else sprintf("%d of %d features at or above %.2f", nsel, length(y), threshold)
        .needlePanel(y, threshold, labels, labelTop, dataNames[m], sub)
        out[[dataNames[m]]] <- setNames(y, labels)
    }
    mtext("Feature", side = 1, line = 0.4, outer = TRUE, col = ink2, cex = 0.8)
    if (nzchar(main)) {
        mtext(main, side = 3, line = 0.8, outer = TRUE, font = 2, col = ink, cex = 0.95)
    }
    invisible(out)
}
