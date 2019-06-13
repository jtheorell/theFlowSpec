#' Select events based on filter
#' 
#' 
#' @param flowObj The fcs object to be subsetted. Both flowFrames and flowSets
#' are accepted. 
#' @param filterName The name of the filter variable in the flow object to be 
#' used for the subsetting. 
#' @param withinFilter Should the events within or outsitde of the filter be 
#' selected?
#' @importFrom flowCore fsApply
#' @export filterOut
filterOut <- function(flowObj, filterName, withinFilter = TRUE){
    if(withinFilter){
        gateVal <- 1
    } else {
        gateVal <- 0
    }
    if(class(flowObj) == "flowSet"){
        resultObj <- fsApply(flowObj, function(x){
            return(x[which(exprs(x[,filterName])[,1] == gateVal),])
        })
    } else if(class(flowObj) == "flowFrame"){
        resultObj <- flowObj[which(exprs(flowObj[,filterName])[,1] == 
                                       gateVal),]
    } else {
        stop("The flowObj needs to be either a flowSet or a flowFrame")
    }
    return(resultObj)
}