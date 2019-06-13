#' @importFrom BiocGenerics colnames
#' @export fixTransNamesAndCoFacs
fixTransNamesAndCoFacs <- function(focusFrame, transNames, transCoFacs,
                                   exclColNmStr){

  # Here, the variables that the normalization should be performed on are
  #selected.
  if(length(transNames) ==1 && transNames == "default"){
    transNames <- BiocGenerics::colnames(focusFrame)
  }

  # Now, the time variable, and possibly other variables are excluded
  excludedTransNames <- grep(paste(exclColNmStr, collapse = "|"),
                             transNames)

  if(missing(excludedTransNames) == FALSE && length(excludedTransNames) > 0){
    transNames <- transNames[-excludedTransNames]
  }

  if(length(transCoFacs) ==1 && transCoFacs == "default"){
    transCoFacs <- massOrFlowTrans(focusFrame = focusFrame,
                                   transNames = transNames)
  } else {
    names(transCoFacs) <- transNames
  }
  return(transCoFacs)
}

