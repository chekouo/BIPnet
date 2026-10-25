## Plot the posterior probabilities of selecting each component, for each data type
## (Result$CompoSelMean from BIP()).
plotCompoSel <- function(Result, dataNames = NULL, which = NULL, type = c("bar", "heatmap"),
                         threshold = 0.5, showValues = TRUE, main = NULL) {
    type <- match.arg(type)
    probs <- Result$CompoSelMean
    if (is.null(probs) || is.null(Result$IndicVar)) {
        stop("Result must be the output of BIP().")
    }
    Np <- nrow(probs)
    r <- ncol(probs)
    IndicVar <- Result$IndicVar

    ## Names of the data types: "Data type m", "Response" and "Covariates" by default
    if (is.null(dataNames)) {
        dataNames <- paste("Data type", seq_len(Np))
        dataNames[IndicVar == 1] <- "Response"
        dataNames[IndicVar == 2] <- if (sum(IndicVar == 2) > 1) paste("Covariates", seq_len(sum(IndicVar == 2))) else "Covariates"
    } else if (length(dataNames) != Np) {
        stop(sprintf("dataNames must have one name per data type (%d).", Np))
    }
    ## Covariates are always kept (probability 1), so they are left out unless asked for
    if (is.null(which)) {
        which <- seq_len(Np)[IndicVar != 2]
    }
    if (!all(which %in% seq_len(Np))) {
        stop(sprintf("which must contain data type numbers between 1 and %d.", Np))
    }
    probs <- probs[which, , drop = FALSE]
    rownames(probs) <- dataNames[which]
    colnames(probs) <- paste0("C", seq_len(r))

    ## Component labels can mix if they switched during the MCMC and were not realigned
    .labelWarning(Result) # defined in plotVarSel.R

    ink <- "#0b0b0b" # primary text
    ink2 <- "#52514e" # secondary text, axes and the threshold line
    grid <- "#e4e3df" # recessive grid
    bar <- "#2a78d6" # single series: one colour, no legend
    if (is.null(main)) {
        main <- "Posterior probability of selecting each component"
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    if (type == "bar") {
        nr <- nrow(probs)
        par(mfrow = c(nr, 1), mar = c(2.6, 4.2, 2.2, 1), oma = c(1.5, 0, if (nzchar(main)) 2.5 else 0, 0),
            mgp = c(2.6, 0.6, 0), las = 1, col.axis = ink2, fg = ink2, cex.axis = 0.85)
        par(cex = 0.75) # same text size whatever the number of panels
        for (d in seq_len(nr)) {
            p <- probs[d, ]
            mids <- barplot(p, ylim = c(0, 1.12), col = NA, border = NA, axes = FALSE, names.arg = rep("", r), space = 0.5)
            abline(h = c(0.25, 0.5, 0.75, 1), col = grid, lwd = 1)
            barplot(p, ylim = c(0, 1.12), col = bar, border = NA, axes = FALSE, names.arg = colnames(probs),
                    space = 0.5, add = TRUE, col.axis = ink2, cex.names = 0.85)
            axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), lwd = 0, lwd.ticks = 0)
            if (!is.null(threshold)) {
                abline(h = threshold, col = ink2, lty = 2, lwd = 1)
            }
            if (showValues) {
                text(mids, p, labels = formatC(p, format = "f", digits = 2), pos = 3, offset = 0.3, cex = 0.8, col = ink, xpd = NA)
            }
            mtext(rownames(probs)[d], side = 3, line = 0.5, adj = 0, font = 2, col = ink, cex = 0.9)
            mtext("Probability", side = 2, line = 2.6, las = 0, col = ink2, cex = 0.75)
        }
        mtext("Component", side = 1, line = 0.3, outer = TRUE, col = ink2, cex = 0.8)
        if (nzchar(main)) {
            mtext(main, side = 3, line = 0.8, outer = TRUE, font = 2, col = ink, cex = 1)
        }
    } else {
        ## Single-hue sequential scale (light = 0, dark = 1), the printed values carry the numbers
        pal <- colorRampPalette(c("#f7fafe", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"))(101)
        nr <- nrow(probs)
        par(mar = c(3.5, max(6, max(nchar(rownames(probs))) * 0.55 + 1), if (nzchar(main)) 3 else 1, 1),
            mgp = c(2, 0.5, 0), las = 1, col.axis = ink2, fg = ink2, cex.axis = 0.85)
        image(seq_len(r), seq_len(nr), t(probs[nr:1, , drop = FALSE]), zlim = c(0, 1), col = pal,
              axes = FALSE, xlab = "", ylab = "", useRaster = FALSE)
        ## a thin surface-coloured gap between cells
        abline(v = seq(0.5, r + 0.5, 1), h = seq(0.5, nr + 0.5, 1), col = "#fcfcfb", lwd = 2)
        axis(1, at = seq_len(r), labels = colnames(probs), lwd = 0, lwd.ticks = 0)
        axis(2, at = seq_len(nr), labels = rev(rownames(probs)), lwd = 0, lwd.ticks = 0)
        mtext("Component", side = 1, line = 2, col = ink2, cex = 0.8)
        if (showValues) {
            z <- t(probs[nr:1, , drop = FALSE])
            text(rep(seq_len(r), nr), rep(seq_len(nr), each = r), labels = formatC(as.vector(z), format = "f", digits = 2),
                 col = ifelse(as.vector(z) > 0.55, "#ffffff", ink), cex = 0.85)
        }
        if (nzchar(main)) {
            title(main, col.main = ink, font.main = 2, cex.main = 1, line = 1.2)
        }
    }
    invisible(probs)
}
