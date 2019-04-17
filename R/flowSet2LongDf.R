#' Convert a flowSet to one long dataframe with all identifiers as separate
#' #columns.
#'
#'
#' This function is mainly used for compatibility with DepecheR.
#'
#' @importFrom flowCore read.flowSet fsApply exprs pData phenoData
#' @importFrom plyr ldply
#' @param frames The flowSet to be converted. Can contain one single flowFrame.
#' @param idInfo A list of any number of characteristics that can be derived
#' from the file names. For each entry, a gsub specification of where to find
#' the information in the file name should be added, such as id=""..._|..."".
#'
#' @return A long data frame with one column per PMT/APD (or fluorochrome,
#' depending on the state of the imported files), one for the acquisition date
#' (for fcs files) and one colum for each specified slot above. If no
#' gsub-pattern is provided, only a single column with the full file name will
#' be used to separate the observations from each file.
#' @export flowSet2longDf
flowSet2LongDf <- function(frames, idInfo) {

    flowSetExprs <- data.frame(fsApply(flowSetRaw, exprs))
    nameVector <- names(pData(phenoData(flowSetRaw))$name)

    if (missing(idInfo)) {
        flowSetExprs$names <-
            unlist(retrieveFlowSetNames(nameVector = nameVector,
                specFlowSet = flowSetRaw, gsubpattern = ""))
    } else {
        for (i in 1:length(idInfo)) {
            flowSetExprs$names <-
                unlist(retrieveFlowSetNames(nameVector = nameVector,
                    specFlowSet = flowSetRaw,
                    gsubpattern = idInfo[[i]]))
            colnames(flowSetExprs)[which(colnames(flowSetExprs) == "names")] <-
                names(idInfo)[i]
        }
    }

    dateVector <- as.vector(fsApply(flowSetRaw,
        function(x) x@description$`$DATE`))
    flowSetExprs$acqDate <-
        unlist(retrieveFlowSetNames(nameVector = dateVector, specFlowSet = flowSetRaw,
            gsubpattern = ""))

    return(flowSetExprs)
}
