#' Plotting all variables against a single variable
#'
#'
#' This function is useful both when setting appropriate gates and when the
#' adjustments of the compensation are done
#' @param flowData This is the full dataset, either a flowFrame or a flowSet,
#' that should be plotted. If it has more rows than "nRows", a subsample (with
#' equal contributions from each sample if a flowSet) will be plotted.
#' @param yCol Here, the variable to be plotted against all the others is
#' selected. It can be either a number or the column name of interest.
#' @param nRows The number of rows that will be used to construct the plot.
#' The fewer, the faster, but the resolution also decreases, naturally. Default
#' is 100000.
#' @param plotName If a name different from yCol should be used, it can be added
#' here.
#' @param zeroTrim In the case of CyTOF data, the events at zero can often
#' be so dominant, that all other density variation is dwarfed, and thus
#' invisible. With this command, the events that are zero in both y
#' and x are trimmed for each x separately.
#' @param saveResult Should the result be saved as a file?
#' @return A plot with one 2D-graph for each variable that the y-variable
#' should be plotted against.
#' @importFrom flowCore fsApply exprs
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes facet_wrap geom_hex xlab ylab theme
#' element_blank element_line ggsave scale_fill_gradientn
#' @export oneVsAllPlot
oneVsAllPlot <- function(flowData, yCol=1, nRows = 10000, plotName = "default",
                         zeroTrim = TRUE, saveResult = TRUE) {

    plotExprs <- plotDownSample(flowData, nRows)

    # If the yCol is given as a character, it is converted to the correct number
    #here
    if (class(yCol) == "character") {
        yCol <- which(colnames(plotExprs) == yCol)
    }

    if (plotName == "default") {
        plotName <- colnames(plotExprs)[yCol]
    }

    #Here, a separate y-name is created
    yName <- colnames(plotExprs)[yCol]

    #Now, the data is truncated, to minimize influence by individual events
    #on the display
    plotExprs <- apply(plotExprs, 2, function(x) {
        high <- quantile(x, 0.999)
        low <- quantile(x, 0.001)
        x[x > high] <- high
        x[x < low] <- low
        return(x)
    })

    #Here, the dataset is converted into a long format.
    longPlot <- melt(as.data.frame(plotExprs))
    longPlot$Y <- rep(plotExprs[, yCol], times=ncol(plotExprs))

    #Here, the double negative events are made NA for each x separately
    if(zeroTrim == TRUE){
        longPlot[which(longPlot$value == 0 & longPlot$Y == 0), 2:3] <- c(NA, NA)
    }

    #Now, the plots are created.
    ggplot(longPlot, aes(x = value, y = Y)) +
        facet_wrap(~variable, scales = "free") +
        geom_hex(na.rm=TRUE) +
        scale_fill_gradientn(colours=c("#440154FF", "#38598CFF", "#2D708EFF",
                                       "#25858EFF", "#1E9B8AFF", "#2BB07FFF",
                                       "#51C56AFF", "#85D54AFF", "#C2DF23FF",
                                       "#FDE725FF", "#FDE725FF", "#FDE725FF")) +
        xlab("") + ylab(yName) + theme(
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
    if(saveResult == TRUE){
        ggsave(paste0(plotName, "_vs_all_others.pdf"))
    }

}
