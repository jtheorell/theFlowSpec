#' Convert a set of csv files to one long dataframe with all identifiers as
#' separate columns.
#'
#'
#' This function is mainly used for compatibility with DepecheR. It is analogous
#' to flowSet2Df, but takes files that have been pre-processed in e.g. FlowJo
#' instead of a flowSet.
#'
#' @importFrom flowCore read.flowSet fsApply exprs pData phenoData
#' @importFrom plyr ldply
#' @param filePath The path to the files to be imported.
#' @param idInfo A list of any number of characteristics that can be derived
#' from the file names. For each entry, a gsub specification of where to find
#' the information in the file name should be added, such as id=""..._|..."".
#'
#' @return A long data frame with one column per PMT/APD (or fluorochrome,
#' depending on the state of the imported files), one for the acquisition date
#' (for fcs files) and one colum for each specified slot above. If no
#' gsub-pattern is provided, only a single column with the full file name will
#' be used to separate the observations from each file.
#' @export csvSet2LongDf

csvSet2LongDf <- function(path = ".", idInfo, nRows = "all") {
    nrows <- ifelse(nRows == "all", -1, nRows)

    filenames <- list.files(path = path, full.names = TRUE)

    if (missing(idInfo) == FALSE) {
      # Here, the idInfo is prepended with the path information
      idInfo <- file.path(path, idInfo)
      flowSetExprs <- ldply(filenames, function(x)
        retrieveCsvNames(x, idInfo, nrows = nrows))

    } else {

        flowSetExprs <- ldply(filenames, function(x)
          retrieveCsvNames(x, nrows = nrows))

    }
}

retrieveCsvNames <- function(filename, idInfo, nrows = -1) {
  ret <- read.csv(filename, nrows = nrows)
  if (missing(idInfo) == TRUE) {
    ret$names <- filename
  } else
    for (i in 1:length(idInfo)) {
      ret$names <- gsub(pattern = idInfo[[i]], "\\1", filename)
      colnames(ret)[which(colnames(ret) == "names")] <- names(idInfo)[i]
    }
  ret
}
