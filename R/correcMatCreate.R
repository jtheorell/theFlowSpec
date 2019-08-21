#' Generate a correction matrix for cytometry data analysis
#'
#'
#' This function aids the correct unmix function, to create a symmetrical
#' correction matrix that should be used together with a flowframe to correct
#' the errors of unmixing.
#' @param flowObj The dataset that should be corrected.
#' @return A symmetrical matrix of zeros with the right row- and column names.
#' @importFrom BiocGenerics ncol colnames
#' @export correcMatCreate
correcMatCreate <- function(flowObj, corrNames = "default"){
  if(class(flowObj) == "flowSet"){
    flowObj <- flowObj[[1]]
  }
  if(corrNames == "default"){
    allColNames <- BiocGenerics::colnames(flowObj)
    corrNames <- allColNames[-which(grepl("ime", allColNames) |
                                     grepl("SC", allColNames))]
  }
    correcMat <- matrix(0, nrow = length(corrNames),
                        ncol = length(corrNames),
                        dimnames =
                          list(corrNames, corrNames))
    diag(correcMat) <- 1
    return(correcMat)
}
