#' Plotting all variables against a single variable
#'
#'
#' This function is useful both when setting appropriate gates and when the adjustments of the compensation are done
#' @param dataset This is the full dataset, either a dataframe or a matrix, that should be plotted. If it has more rows than "nrows", a random subsample (without resampling) will be plotted
#' @param yCol Here, the variable to be plotted against all the others is selected. It can be either a number or the column name of interest.
#' @param nRows The number of rows that will be plotted. The fewer, the faster, but the resolution also decreases, naturally. Default is 10000
#' @param color The color for the density peaks in the plots.
#' @param markerName If a specific name for the parameter in the graph name is wanted, it should be added here. Default is name(yCol)_vs_all_others
#' @param saveResult If you do not want to save the results, and view them directly on the screen instead, you set this to false.
#' @return A graph with one sub-graph for each variable that the y-variable should be plotted against.
#'
#' @export oneVsAllPlot
oneVsAllPlot <- function(dataset, yCol, nRows = 10000, markerName = "default", color = "red", saveResult = TRUE) {
    if (class(dataset) == "matrix") {
        dataset <- as.data.frame(dataset)
    }


    if (nrow(dataset) > nRows) {
        datasetInUse <- dataset[sample(1:nrow(dataset), size = nRows), ]
    } else {
        warning(paste0("The number of rows in the dataset was only ", nrow(dataset), ", so the output will be restricted to this"))
        datasetInUse <- dataset
    }

    # If the yCol is given as a character, it is converted to the correct number here
    if (class(yCol) == "character") {
        yCol <- which(colnames(dataset) == yCol)
    }

    if (markerName == "default") {
        markerName <- colnames(dataset)[yCol]
    }


    # Here, the yCol is separated from the rest of the data
    yColData <- datasetInUse[, yCol]
    nonYColData <- datasetInUse[, -yCol]

    # now, dotsize is decided
    dotSize <- 50 / sqrt(nrow(datasetInUse))

    # Now, the plots are created. First, they are divided in sections of 4, as the first one in every quartette will have different axes.


    if (saveResult == TRUE) {
        png(width = 1300, height = 1300, paste0(markerName, "_vs_all_others.png"))
    }
    par(mfrow = c(6, 6), mai = c(0, 0, 0, 0), pty = "m", mar = c(3, 0, 0, 0), oma = c(6, 4, 4, 6), cex = 1.2, mgp = c(1, 0.5, 0))
    for (i in 1:ncol(nonYColData)) {
        if (i == 1 || (i - 1) / 6 == round((i - 1) / 6)) {
            oneVsAllPlotCoFunction(xVar = nonYColData[, i], yVar = yColData, color = color, yaxt = "s", dotSize = dotSize, xlab = "", ylab = "")
            title(xlab = colnames(nonYColData)[i], line = 1.5)
        } else {
            oneVsAllPlotCoFunction(xVar = nonYColData[, i], yVar = yColData, color = color, xlab = "", ylab = "", yaxt = "n", dotSize = dotSize)
            title(xlab = colnames(nonYColData)[i], ylab = "", line = 1.5)
        }
        mtext(paste0(colnames(datasetInUse)[yCol], " vs all other markers"), side = 3, line = 2, outer = TRUE, cex = 2)
    }
    if (saveResult == TRUE) {
        dev.off()
    }

}
