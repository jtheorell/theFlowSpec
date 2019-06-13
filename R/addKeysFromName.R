#' Dissociate the name of flow frames into multiple keys
#' 
#' Adds new keys to flow objects by separating the flowFrame 
#' identifiers/sample names into multiple entities. 
#' 
#' @importFrom flowCore fsApply keyword
#' @param flowObj The fcs object to be transformed. Both flowFrames and flowSets
#' are accepted. 
#' @param keys A list of any number of characteristics that can be derived
#' from the sample names. For each entry, a gsub specification of where to find
#' the information in the file name should be added, such as id="..._|".
#'
#' @return The appended flowObj
#' @export addKeysFromName
addKeysFromName <- function(flowObj, keys){
    if(class(flowObj) == "flowSet"){
        resultObj <- fsApply(flowObj, addKeysFromNameCoFunction, keys=keys)
    } else if(class(flowObj) == "flowFrame") {
        resultObj <- addKeysFromNameCoFunction(focusFrame = flowObj, 
                                               keys = keys)
        } else {
        stop("The flowObj needs to be either a flowSet or a flowFrame")
    }
    
}

addKeysFromNameCoFunction <- function(focusFrame, keys){
    
    newKeys <- lapply(keys, function(x) gsub(x, "\\1", identifier(focusFrame)))
    
    keyword(focusFrame) <- c(keyword(focusFrame), newKeys)
    return(focusFrame)
}