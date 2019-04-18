#' @importFrom BiocGenerics nrow
plotDownSample <- function(flowData, nRows = 10000){
    if(class(flowData)=="flowSet"){
        if(sum(fsApply(flowData, BiocGenerics::nrow)) > nRows){
            numIncluded <- floor(nRows/length(flowData))
            plotExprs <- fsApply(flowData, function(x) 
                return(exprs(x)[sample(seq_len(BiocGenerics::nrow(x)), 
                                       numIncluded),]))
            
        } else {
            plotExprs <- fsApply(flowData, function(x) return(exprs(x)))
        }
    } else if(BiocGenerics::nrow(flowData) > nRows){
        plotExprs <- 
            exprs(flowData)[sample(seq_len(BiocGenerics::nrow(flowData)), 
                                   nRow),]
    } else {
        plotExprs <- exprs(flowData)
        message("The number of rows in the dataset was only ", 
                nrow(plotExprs), 
                ", so the output will be restricted to this")
    }
    return(plotExprs)
}