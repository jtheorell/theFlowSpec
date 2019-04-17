#' @importFrom flowCore fsApply exprs
retrieveFlowSetNames <- function(nameVector, specFlowSet, gsubpattern) {
    lengthList <- as.list(fsApply(specFlowSet, function(x) nrow(x@exprs)))
    nameVectorShort <- as.list(gsub(pattern = gsubpattern, "\\1", nameVector))
    return(as.vector(mapply(function(x, y) rep(x, times = y), nameVectorShort, lengthList)))
}
