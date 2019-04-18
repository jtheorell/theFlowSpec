#' Plotting all variables as histograms
#'
#'
#' This function individually plots all histograms. 
#' @param flowData This is the full dataset, either a flowFrame or a flowSet, 
#' that should be plotted. If it has more rows than "nRows", a subsample (with 
#' equal contributions from each sample if a flowSet) will be plotted.
#' @param nRows The number of rows that will be used to construct the plot. 
#' The fewer, the faster, but the resolution also decreases, naturally. Default
#' is 100000.
#' @param saveResult Should the result be saved as a file? 
#' @return A plot with one histogram for each column in the flowData.
#' @importFrom flowCore fsApply exprs
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes facet_wrap geom_density theme element_blank 
#' element_line ggsave
#' @export
allHist <- function(flowData, plotName = "all_histograms", nRows = 100000, 
                    saveResult = TRUE) {
    
    plotExprs <- plotDownSample(flowData, nRows)
    
    longPlot <- melt(plotExprs)
    ggplot(longPlot, aes(value)) + 
        facet_wrap(~variable, scales = "free") + 
        geom_density() + theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
    ggsave(paste0(name, ".pdf"))
}
