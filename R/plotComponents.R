## Scatter plot of two estimated shared components (columns of Result$EstU, or of other
## component scores such as BIPpredict()$Upredtest), with the samples coloured by the outcome.
## Uses .labelWarning() from plotVarSel.R.
plotComponents <- function(Result, components = c(1, 2), y = NULL, U = NULL,
                           scale = c("sequential", "diverging"), orient = FALSE, arrow = TRUE,
                           outcomeName = "Outcome", main = NULL) {
    scale <- match.arg(scale)
    if (is.null(Result$EstU) || is.null(Result$IndicVar)) {
        stop("Result must be the output of BIP().")
    }
    if (is.null(U)) {
        U <- Result$EstU
    }
    U <- as.matrix(U)
    if (is.null(y)) {
        y <- Result$Y
        if (is.null(y)) {
            stop("Give the outcome values in y (this BIP() output does not contain the response).")
        }
    }
    y <- as.vector(y)
    r <- ncol(U)
    if (length(components) != 2 || components[1] == components[2] || !all(components %in% seq_len(r))) {
        stop(sprintf("components must be two different component numbers between 1 and %d.", r))
    }
    if (length(y) != nrow(U)) {
        stop(sprintf("y must have one value per sample (row of U): %d, not %d.", nrow(U), length(y)))
    }
    if (any(!is.finite(y))) {
        stop("y must not contain missing or infinite values.")
    }
    .labelWarning(Result)

    np <- which(Result$IndicVar == 1)
    ay <- Result$EstLoad[[np]][, 1] # loading of each component on the outcome (model of the selected features)
    pout <- Result$CompoSelMean[np, ] # probability that each component is linked to the outcome
    ## a component's sign is arbitrary (U and A can flip together): optionally orient each one so
    ## that its loading on the outcome is positive, as in plotLoadings()
    s <- if (orient) ifelse(ay < 0, -1, 1) else rep(1, length(ay))
    k1 <- components[1]
    k2 <- components[2]
    u1 <- U[, k1] * s[k1]
    u2 <- U[, k2] * s[k2]
    a1 <- ay[k1] * s[k1]
    a2 <- ay[k2] * s[k2]

    ## colours: one blue hue from light (low) to dark (high), starting at a step that stays visible
    ## on the white surface; or red (below the mean) - gray (mean) - blue (above the mean)
    if (scale == "sequential") {
        pal <- colorRampPalette(c("#86b6ef", "#3987e5", "#1c5cab", "#0d366b"))(100)
        lo <- min(y)
        hi <- max(y)
    } else {
        pal <- colorRampPalette(c("#8f2423", "#e34948", "#f09a97", "#a9a8a3", "#86b6ef", "#2a78d6", "#0d366b"))(100)
        half <- max(abs(y - mean(y)))
        lo <- mean(y) - half
        hi <- mean(y) + half
    }
    if (hi <= lo) { # all outcome values equal: give the colour key a range around that value
        lo <- lo - 0.5
        hi <- hi + 0.5
    }
    idx <- pmin(100, pmax(1, 1 + floor(99.999 * (y - lo) / (hi - lo))))
    fill <- pal[idx]

    ink <- "#0b0b0b"
    ink2 <- "#52514e"
    if (is.null(main)) {
        main <- sprintf("Shared components C%d and C%d, samples coloured by %s", k1, k2, tolower(outcomeName))
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(matrix(1:2, 1, 2), widths = c(1, lcm(3.2))) # fixed key width, whatever the size of the device

    par(mar = c(4.6, 4.6, if (nzchar(main)) 3.6 else 1.5, 0.8), mgp = c(2.6, 0.5, 0), las = 1, col.axis = ink2, fg = ink2, cex.axis = 0.8)
    lim <- max(abs(c(u1, u2))) * 1.08
    plot(u1, u2, type = "n", xlim = c(-lim, lim), ylim = c(-lim, lim), asp = 1, axes = FALSE, xlab = "", ylab = "")
    ticks <- pretty(c(-lim, lim), n = 6)
    abline(v = ticks, h = ticks, col = "#eeede9", lwd = 1)
    abline(v = 0, h = 0, col = "#d6d5d0", lwd = 1)
    axis(1, at = ticks, lwd = 0, lwd.ticks = 0)
    axis(2, at = ticks, lwd = 0, lwd.ticks = 0)
    ## points with a thin surface-coloured ring, so that overlapping samples stay distinct
    points(u1, u2, pch = 21, bg = fill, col = "#fcfcfb", cex = 1.35, lwd = 0.8)
    if (arrow && (abs(a1) + abs(a2)) > 1e-10) {
        len <- 0.75 * lim / sqrt(a1^2 + a2^2)
        arrows(0, 0, a1 * len, a2 * len, length = 0.1, angle = 25, lwd = 2, col = ink)
        text(a1 * len, a2 * len, labels = paste(outcomeName, "increases"), pos = if (a2 >= 0) 3 else 1, cex = 0.75, col = ink, xpd = NA)
    }
    lab <- function(k, a) sprintf("C%d   (outcome loading %+.2f, linked to the outcome with probability %.2f)", k, a, pout[k])
    mtext(lab(k1, a1), side = 1, line = 2.4, col = ink2, cex = 0.75)
    mtext(lab(k2, a2), side = 2, line = 2.6, las = 0, col = ink2, cex = 0.75)
    if (nzchar(main)) {
        title(main, col.main = ink, font.main = 2, cex.main = 0.95, line = 2)
        note <- paste0(if (identical(U, as.matrix(Result$EstU))) "Posterior means of the components" else "Component scores given in U",
                       if (orient) "; signs oriented so that each loading on the outcome is positive" else "")
        mtext(note, side = 3, line = 0.7, col = ink2, cex = 0.72)
    }

    ## colour key
    par(mar = c(4.6, 0.6, if (nzchar(main)) 3.6 else 1.5, 3.2))
    vals <- seq(lo, hi, length.out = 100)
    image(1, vals, matrix(vals, nrow = 1), col = pal, axes = FALSE, xlab = "", ylab = "")
    kt <- pretty(c(lo, hi), n = 5)
    kt <- kt[kt >= lo & kt <= hi]
    axis(4, at = kt, labels = format(kt), lwd = 0, lwd.ticks = 1, tcl = -0.2, cex.axis = 0.75, mgp = c(2, 0.35, 0), las = 1)
    mtext(outcomeName, side = 3, line = 0.5, col = ink2, cex = 0.75)

    out <- data.frame(u1, u2, y)
    names(out) <- c(paste0("C", k1), paste0("C", k2), "y")
    invisible(out)
}
