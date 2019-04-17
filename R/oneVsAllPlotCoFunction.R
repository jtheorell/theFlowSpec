oneVsAllPlotCoFunction <- function(xVar, yVar, color, xlab, ylab, yaxt = "s", dotSize = dotSize) {

    cols <- colorRampPalette(c("black", "grey", color))(256)
    varDf <- data.frame(xVar, yVar)

    ## Use densCols() output to get density at each point. The colors here are only supporting the coming order of the rows further down the script.
    x <- densCols(xVar, yVar, colramp = colorRampPalette(c("black", "white")))
    varDf$dens <- col2rgb(x)[1, ] + 1L
    varDf$col <- cols[varDf$dens]
    plot(yVar ~ xVar, data = varDf[order(varDf$dens), ], main = NULL, pch = 20, cex = dotSize, col = col, xlab = xlab, ylab = ylab, yaxt = yaxt)
}
