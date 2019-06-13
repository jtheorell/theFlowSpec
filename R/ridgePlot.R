#' Density plots showing gates for all files in a flowSet
#'
#' This function is completely analogous to the the densityplot function in
#' flowViz, but it takes the infor straight from the flowFrame, without the need
#' for any filter functions or objects, and the interface is simpler.
#' @param focusSet The flow set in focus.
#' @param plotVar The variable that should be plotted.
#' @param filtVar The variable where the filter is present.
#' @param sampleKey The name on the keyword from which the names of the samples
#' can be extracted.
#' @param plotName If a name different from plotVar should be used, it can be
#' added here.
#' @param saveResult Should the result be saved as a file?
#' @importFrom viridis scale_fill_viridis
#' @importFrom ggplot2 ggplot aes geom_segment scale_y_discrete
#' @importFrom ggridges geom_density_ridges_gradient theme_ridges
#' @export ridgePlot
ridgePlot <- function(focusSet, plotVar, filtVar, sampleKey = "GUID",
                      plotName = plotVar, saveResult = TRUE){

    #First, the long frame of interest is constructed
    focusExprsList <- fsApply(focusSet, function(x)
        ridgePlotCoFunction1(x[,c(plotVar)], sampleKey = sampleKey))

    longExprs <- do.call("rbind", focusExprsList)

    exprsLinesList <- fsApply(focusSet, function(x){
        focFrame <- x[,c(plotVar, filtVar)]
        x1 <- min(exprs(focFrame[which(exprs(focFrame[,2]) == 1),1]))
        x2 <- max(exprs(focFrame[which(exprs(focFrame[,2]) == 1),1]))
        if(x1 == min(exprs(focFrame[,1]))){
          x1 <- NA
        }
        if(x2 == max(exprs(focFrame[,1]))){
          x2 <- NA
        }
        return(c(x1, x2))
    }, simplify = FALSE)

    exprsLines <- as.data.frame(t(do.call(cbind, exprsLinesList)))
    colnames(exprsLines) <- c("x1", "x2")
    exprsLines$sampleKey <- unique(longExprs$sampleKey)

    ggplot(longExprs, aes(x=longExprs[,1], y=sampleKey, fill=..x..)) +
        geom_density_ridges_gradient() +
        geom_segment(data = exprsLines, aes(x = x1, xend = x1,
                                            y = as.numeric(sampleKey),
                                            yend = as.numeric(sampleKey) + .9),
                     color = "red", na.rm=TRUE) +
        geom_segment(data = exprsLines, aes(x = x2, xend = x2,
                                            y = as.numeric(sampleKey),
                                            yend = as.numeric(sampleKey) + .9),
                     color = "red", na.rm=TRUE) +
        scale_y_discrete(expand = c(0.01, 0)) +
        theme_ridges(grid = FALSE, center = TRUE) +
        xlab(plotVar) + ylab("Sample id")

    if(saveResult == TRUE){
        ggsave(paste0(plotName, "_ridgePlot.pdf"))
    }



}

ridgePlotCoFunction1 <- function(focusFrame, sampleKey){
    if(BiocGenerics::nrow(focusFrame) > 50000){
        focusExprs <-
            as.data.frame(exprs(focusFrame[sample(
                seq_len(BiocGenerics::nrow(focusFrame)), 50000),]))
    } else {
        focusExprs <- as.data.frame(exprs(focusFrame))
    }

    focusExprs$sampleKey <- as.factor(keyword(focusFrame)[sampleKey])
    return(focusExprs)
}

ridgePlotCoFunction2 <- function(focusFrame){

    x1 <- min(exprs(focusFrame[which(exprs(focusFrame[,2]) == 1),1]))
    x2 <- max(exprs(focusFrame[which(exprs(focusFrame[,2]) == 1),1]))
    return(c(x1, x2))
}
