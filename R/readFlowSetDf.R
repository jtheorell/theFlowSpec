#' Import a flow set and convert it to one long dataframe
#'
#'
#' This function is useful to manipulate data within R, as, in most instances, the majority of the meta information stored in the fcs file is not of general use for day.to.day analyses.
#'
#' @importFrom flowCore read.flowSet fsApply exprs pData phenoData
#' @importFrom plyr ldply
#' @param path The path to the folder where the .fcs or .csv files are stored.
#' @param idInfo A list of any number of characteristics that can be derived from the file names. For each entry, a gsub specification of where to find the information in the file name should be added, such as id=""..._|..."".
#'
#' @return A long data frame with one column per PMT/APD (or fluorochrome, depending on the state of the imported files), one for the acquisition date (for fcs files) and one colum for each specified slot above. If no gsub-pattern is provided, the full file name will be used to separate the observations.
#'
#' @export readFlowSetDf
readFlowSetDf <- function(path=".", idInfo, nRows="all"){

  if(grepl(".fcs", dir(path)[1])){
    if(nRows=="all"){
      which.lines=NULL
    } else {
      which.lines <- nRows
    }
    flowSetRaw <- read.flowSet(path=path, alter.names=TRUE, transformation=FALSE,  which.lines=which.lines, phenoData=list(name="FILENAME"))
    flowSetExprs <- data.frame(fsApply(flowSetRaw, exprs))
    nameVector <- names(pData(phenoData(flowSetRaw))$name)

    if(missing(idInfo)){
      flowSetExprs$names <- unlist(retrieveFlowSetNames(nameVector=nameVector, specFlowSet=flowSetRaw, gsubpattern=""))
    } else {
      for(i in 1:length(idInfo)){
        flowSetExprs$names <- unlist(retrieveFlowSetNames(nameVector=nameVector, specFlowSet=flowSetRaw, gsubpattern=idInfo[[i]]))
        colnames(flowSetExprs)[which(colnames(flowSetExprs)=="names")] <- names(idInfo)[i]
      }
    }

    dateVector <- as.vector(fsApply(flowSetRaw, function(x) x@description$`$DATE`))
    flowSetExprs$acqDate  <- unlist(retrieveFlowSetNames(nameVector=dateVector, specFlowSet=flowSetRaw, gsubpattern=""))

  } else {
    if(nRows=="all"){
      nrows=-1
    } else {
      nrows <- nRows
    }
    filenames <- list.files(path = path, pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE)
    filenames <- paste0(path, "/", filenames)
    #Here, the idInfo is prepended with the path information
    idInfo <- lapply(idInfo, function(x) paste0(path, "/", x))

    if(missing(idInfo)==FALSE){
      flowSetExprs <- ldply(filenames, function(x) retrieveCsvWNames(x, idInfo, nrows=nrows))

    } else {
      flowSetExprs <- ldply(filenames, function(x) retrieveCsvWNames(x, nrows=nrows))

    }

  }

  return(flowSetExprs)
}
