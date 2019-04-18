#' @importFrom flowCore parameters 
#' @importFrom Biobase pData
#' @importFrom BiocGenerics colnames
#' @export heliosChannelCleanup
heliosChannelCleanup <- function(flowFrameRaw, nonMetalChannels =
                                    c("Time", "Event_length", "Center",
                                      "Offset", "Width", "Residual")){

    #First, the phenoData object is retrieved
    pDataFF <- pData(parameters(flowFrameRaw))

    #And an order object is created
    pDataOrder <- seq_len(nrow(pDataFF))

    #Now, all non-metal channels are excluded
    focusNames <- pDataFF$desc[-which(pDataFF$name %in% nonMetalChannels)]

    #And here, rather cheekily, all channels that lack a "_" in their
    #description are excluded.

    includedMetalChannels <- focusNames[grep("_", focusNames)]

    includedNumbers <- c(which(pDataFF$desc %in% includedMetalChannels),
                         which(pDataFF$name %in% nonMetalChannels))

    includedNumbersOrdered <- includedNumbers[order(includedNumbers)]

    #Now, the file and the pheno data are reduced to this subset
    flowFrameReduced <- flowFrameRaw[,includedNumbersOrdered]

    #And here, some adjustments are made to the phenoData, to make it work
    #better with R.

    pDataFFReduced <- pDataFF[includedNumbersOrdered,]

    newColNames <- gsub(".+_", "\\1", unname(pDataFFReduced$desc))

    #Now, all potential lines are exchanged for dots.
    newColNames <- gsub("-", ".", newColNames)

    newColNames[which(pDataFFReduced$name %in% nonMetalChannels)] <-
        nonMetalChannels

    pDataFFReduced$oldName <- pDataFFReduced$name
    pDataFFReduced$name <- newColNames

    BiocGenerics::colnames(flowFrameReduced) <- newColNames

    pData(parameters(flowFrameReduced)) <- pDataFFReduced

    return(flowFrameReduced)

}
