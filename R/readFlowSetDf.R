#' Import a flow set and convert it to one long dataframe
#'
#'
#' This function is useful to manipulate data within R, as, in most instances, the majority of the meta information stored in the fcs file is not of general use for day.to.day analyses.
#'
#' @importFrom flowCore read.flowSet fsApply exprs pData phenoData
#' @param path The path to the folder where the .fcs files are stored.
#' @param id If a column specifying which of the original donors that each observation comes from should be included, the gsub-pattern should be specified here.
#' @param group See "ids"
#' @param stim See "ids"
#'
#' @return A long data frame with one column per PMT/APD, one for the acquisition date and one colum for each specified slot above. If the same donor id occurs for two different stimuli, another column creating unique ids is also made. If no gsub-pattern is provided, the full file name will be used to separate the observations.
#'
#' @export readFlowSetDf
readFlowSetDf <- function(path=".", id, group, stim){
  oldWd <- getwd()
  setwd(path)

  flowSetRaw <- read.flowSet(path=path, alter.names=TRUE, transformation=FALSE,  which.lines=NULL, phenoData=list(name="FILENAME"))
  flowSetExprs <- data.frame(fsApply(flowSetRaw, exprs))
  nameVector <- names(pData(phenoData(flowSetRaw))$name)

  if(missing(id) && missing(group) && missing(stim)){
    flowSetExprs$names <- unlist(retrieveFlowSetNames(nameVector=nameVector, specFlowSet=flowSetRaw, gsubpattern=""))
  }
  if(missing(id)==FALSE){
    flowSetExprs$id <- unlist(retrieveFlowSetNames(nameVector=nameVector, specFlowSet=flowSetRaw, gsubpattern=id))
  }
  if(missing(group)==FALSE){
    flowSetExprs$group <- unlist(retrieveFlowSetNames(nameVector=nameVector, specFlowSet=flowSetRaw, gsubpattern=group))

  }
  if(missing(group)==FALSE){
    flowSetExprs$group <- unlist(retrieveFlowSetNames(nameVector=nameVector, specFlowSet=flowSetRaw, gsubpattern=group))

  }

  dateVector <- as.vector(fsApply(flowSetRaw, function(x) x@description$`$DATE`))
  flowSetExprs$acqDate  <- unlist(retrieveFlowSetNames(nameVector=dateVector, specFlowSet=flowSetRaw, gsubpattern=""))

  return(flowSetExprs)
  setwd(oldWd)
}
