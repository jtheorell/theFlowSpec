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
#' @export flowSet2Df

csvSet2LongDf <- function(path = ".", idInfo, nRows = "all") {
    if (nRows == "all") {
        nrows = -1
    } else {
        nrows <- nRows
    }
    filenames <- file.path(path, list.files(path = path, pattern = NULL,
        all.files = FALSE,
        full.names = FALSE,
        recursive = FALSE))
    filenames <- paste0(path, "/", filenames)
    # Here, the idInfo is prepended with the path information
    idInfo <- lapply(idInfo, function(x) paste0(path, "/", x))

    if (missing(idInfo) == FALSE) {
        flowSetExprs <- ldply(filenames, function(x) retrieveCsvWNames(x, idInfo, nrows = nrows))

    } else {
        flowSetExprs <- ldply(filenames, function(x) retrieveCsvWNames(x, nrows = nrows))

    }
}
