#' @export creEmptyCorrMat
creEmptyCorrMat <- function(unmixData){
  correcMat <- matrix(0, nrow=ncol(unmixData), ncol=ncol(unmixData), dimnames = list(colnames(unmixData), colnames(unmixData)))
  diag(correcMat) <- 1
  return(correcMat)
}
