## Plot the estimated loadings of one data type, by component: either the loadings of the model
## of the selected features (Result$EstLoad) or their Bayesian model average (Result$EstLoadModel
## weighted by Result$PostGam). Uses .featureLabels() and .labelWarning() from plotVarSel.R.

## Loadings (components x features) of data type m, and of the response, under the chosen model
.loadings <- function(Result, m, model) {
    np <- which(Result$IndicVar == 1)
    if (model == "selected") {
        return(list(L = Result$EstLoad[[m]], y = Result$EstLoad[[np]][, 1]))
    }
    nb <- Result$nbrmodel
    w <- Result$PostGam[seq_len(nb)]
    L <- Reduce(`+`, lapply(seq_len(nb), function(k) w[k] * Result$EstLoadModel[[k]][[m]]))
    y <- Reduce(`+`, lapply(seq_len(nb), function(k) w[k] * Result$EstLoadModel[[k]][[np]][, 1]))
    list(L = L, y = y)
}

## Diverging scale: red (negative) - gray (0) - blue (positive), equal steps on each arm
.divergingPalette <- function(n = 101) {
    colorRampPalette(c("#8f2423", "#e34948", "#f4b6b3", "#f0efec", "#b7d3f6", "#2a78d6", "#0d366b"))(n)
}

## One panel: a thin line per feature from 0 to its loading, blue if positive, red if negative
.signedPanel <- function(y, lim, labels, labelTop, heading, subheading) {
    ink <- "#0b0b0b"
    ink2 <- "#52514e"
    P <- length(y)
    x <- seq_len(P)
    nz <- abs(y) > 1e-10
    pos <- nz & y > 0
    neg <- nz & y < 0
    plot(x, y, type = "n", ylim = c(-lim, lim) * 1.18, xlim = c(0.5, P + 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
    ticks <- pretty(c(-lim, lim), n = 4)
    ticks <- ticks[abs(ticks) <= lim * 1.05]
    abline(h = ticks[ticks != 0], col = "#e4e3df", lwd = 1)
    lw <- if (P > 300) 1 else 2
    if (any(pos)) segments(x[pos], 0, x[pos], y[pos], col = "#2a78d6", lwd = lw, lend = 1)
    if (any(neg)) segments(x[neg], 0, x[neg], y[neg], col = "#e34948", lwd = lw, lend = 1)
    abline(h = 0, col = ink2, lwd = 1)
    axis(2, at = ticks, labels = format(ticks), lwd = 0, lwd.ticks = 0)
    axis(1, lwd = 0, lwd.ticks = 1, col.ticks = "#b8b7b2", tcl = -0.25)
    ## label the largest |loadings|: above the positive ones, below the negative ones; a label is
    ## skipped if it would overlap one already placed on the same side
    if (labelTop > 0 && any(nz)) {
        usr <- par("usr")
        half <- 0.55 * strwidth(labels, cex = 0.7)
        lx <- pmin(pmax(x, usr[1] + half), usr[2] - half)
        cand <- order(-abs(y), x)
        cand <- cand[nz[cand]]
        keep <- integer(0)
        for (j in cand) {
            if (length(keep) >= labelTop) break
            same <- keep[sign(y[keep]) == sign(y[j])]
            if (all(abs(lx[j] - lx[same]) > half[j] + half[same])) keep <- c(keep, j)
        }
        text(lx[keep], y[keep], labels = labels[keep], pos = ifelse(y[keep] > 0, 3, 1), offset = 0.3, cex = 0.7, col = ink, xpd = NA)
    }
    mtext(heading, side = 3, line = 0.55, adj = 0, font = 2, col = ink, cex = 0.9)
    mtext(subheading, side = 3, line = 0.55, adj = 1, col = ink2, cex = 0.75)
    mtext("Loading", side = 2, line = 2.6, las = 0, col = ink2, cex = 0.75)
}

plotLoadings <- function(Result, which = 1, components = NULL, type = c("index", "heatmap"),
                         model = c("selected", "BMA"), featureNames = NULL, labelTop = 5,
                         maxFeatures = 50, orient = FALSE, dataName = NULL, main = NULL) {
    type <- match.arg(type)
    model <- match.arg(model)
    if (is.null(Result$EstLoad) || is.null(Result$IndicVar)) {
        stop("Result must be the output of BIP().")
    }
    IndicVar <- Result$IndicVar
    if (length(which) != 1 || !(which %in% seq_along(IndicVar)) || IndicVar[which] == 1) {
        stop("which must be the number of one data type that is not the response.")
    }
    est <- .loadings(Result, which, model)
    L <- est$L # components x features
    yload <- est$y # loading of each component on the response
    r <- nrow(L)
    P <- ncol(L)
    if (is.null(components)) {
        components <- seq_len(r)
    }
    if (!all(components %in% seq_len(r))) {
        stop(sprintf("components must be component numbers between 1 and %d.", r))
    }
    .labelWarning(Result)
    ## a component's sign is arbitrary (U and A can flip together): optionally orient each one so
    ## that its loading on the response is positive
    if (orient) {
        s <- ifelse(yload < 0, -1, 1)
        L <- L * s
        yload <- yload * s
    }
    labels <- .featureLabels(Result, which, featureNames, P)
    compo <- Result$CompoSelMean[which, ]
    if (is.null(dataName)) {
        dataName <- paste("data type", which)
    }
    if (is.null(main)) {
        main <- sprintf("Loadings of %s, by component", dataName)
    }
    note <- paste0(if (model == "selected") "Model of the selected features" else sprintf("Bayesian model averaging over %d models", Result$nbrmodel),
                   if (orient) "; signs oriented so that each loading on the outcome is positive" else "")
    ink <- "#0b0b0b"
    ink2 <- "#52514e"
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    sub <- L[components, , drop = FALSE]
    lim <- max(abs(sub), 1e-8)

    if (type == "index") {
        par(mfrow = c(length(components), 1), mar = c(2.2, 4.2, 2.2, 1), oma = c(1.6, 0, if (nzchar(main)) 3.4 else 1.2, 0),
            mgp = c(2.6, 0.5, 0), las = 1, col.axis = ink2, fg = ink2, cex.axis = 0.8)
        par(cex = 0.75) # same text size whatever the number of panels
        for (k in components) {
            nnz <- sum(abs(L[k, ]) > 1e-10)
            subhead <- sprintf("component selected with probability %.2f  |  loading on the outcome %+.2f  |  %d non-zero loading%s",
                               compo[k], yload[k], nnz, if (nnz == 1) "" else "s")
            .signedPanel(L[k, ], lim, labels, labelTop, paste0("C", k), subhead)
        }
        mtext("Feature", side = 1, line = 0.4, outer = TRUE, col = ink2, cex = 0.8)
        if (nzchar(main)) {
            mtext(main, side = 3, line = 1.7, outer = TRUE, font = 2, col = ink, cex = 0.95)
        }
        mtext(note, side = 3, line = 0.4, outer = TRUE, col = ink2, cex = 0.75)
        out <- t(sub)
        dimnames(out) <- list(labels, paste0("C", components))
        return(invisible(out))
    }

    ## heatmap: features with a non-zero loading in at least one plotted component, largest first
    best <- apply(abs(sub), 2, max)
    rows <- which(best > 1e-10)
    rows <- rows[order(-best[rows], rows)]
    if (length(rows) > maxFeatures) {
        rows <- rows[seq_len(maxFeatures)]
    }
    if (length(rows) == 0) {
        plot.new()
        text(0.5, 0.5, sprintf("All loadings of %s are 0 in these components.", dataName), col = ink2)
        return(invisible(matrix(numeric(0), 0, length(components))))
    }
    z <- t(sub[, rows, drop = FALSE]) # features x components
    dimnames(z) <- list(labels[rows], paste0("C", components))
    nr <- nrow(z)
    nc <- ncol(z)
    lim <- max(abs(z))
    pal <- .divergingPalette(101)
    cexr <- if (nr > 30) 0.55 else 0.75
    left <- max(5, max(nchar(rownames(z))) * 0.5 * cexr / 0.75 + 1)
    layout(matrix(1:2, 2, 1), heights = c(1, lcm(2.6))) # fixed key height, whatever the size of the device
    par(mar = c(6.2, left, if (nzchar(main)) 4 else 2, 1), mgp = c(2, 0.4, 0), las = 1, col.axis = ink2, fg = ink2)
    image(seq_len(nc), seq_len(nr), t(z[nr:1, , drop = FALSE]), zlim = c(-lim, lim), col = pal, axes = FALSE, xlab = "", ylab = "")
    abline(v = seq(0.5, nc + 0.5, 1), col = "#fcfcfb", lwd = 2)
    if (nr <= 60) abline(h = seq(0.5, nr + 0.5, 1), col = "#fcfcfb", lwd = 1)
    axis(2, at = seq_len(nr), labels = rev(rownames(z)), lwd = 0, lwd.ticks = 0, cex.axis = cexr)
    axis(1, at = seq_len(nc), labels = paste0("C", components), lwd = 0, lwd.ticks = 0, cex.axis = 0.85, font = 2, col.axis = ink)
    axis(1, at = seq_len(nc), labels = sprintf("selected %.2f\noutcome %+.2f", compo[components], yload[components]),
         lwd = 0, lwd.ticks = 0, cex.axis = 0.7, line = 1.6)
    mtext("Component: probability of selecting it, and its loading on the outcome", side = 1, line = 5, col = ink2, cex = 0.75)
    if (nzchar(main)) {
        title(main, col.main = ink, font.main = 2, cex.main = 0.95, line = 2.6)
    }
    mtext(sprintf("%s. %d feature%s with a non-zero loading%s", note, length(rows), if (length(rows) == 1) "" else "s",
                  if (length(rows) == maxFeatures) sprintf(" (the %d largest shown)", maxFeatures) else ""),
          side = 3, line = 0.6, col = ink2, cex = 0.7)
    ## colour key
    par(mar = c(1.6, left, 1.3, 1))
    image(seq(-lim, lim, length.out = 101), 1, matrix(seq(-lim, lim, length.out = 101), ncol = 1), col = pal, axes = FALSE, xlab = "", ylab = "")
    kt <- pretty(c(-lim, lim), n = 4)
    kt <- kt[abs(kt) <= lim]
    axis(1, at = kt, labels = format(kt), lwd = 0, lwd.ticks = 1, tcl = -0.2, cex.axis = 0.7, mgp = c(2, 0.2, 0))
    mtext("Loading (red negative, blue positive)", side = 3, line = 0.25, col = ink2, cex = 0.7)
    invisible(z)
}
