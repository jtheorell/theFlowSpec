#' importFrom stats density
#' @export gating1D
gating1D <- function(dataSet, variable, outputDataset, gateName=variable, color="red", saveResult=TRUE, makeDir=TRUE){

  d <- density(dataSet[,variable])
  plot(d, main=gateName, xlab=variable)
  polygon(d, col="grey", border="black")

  coordinates <- locator(type="n")

  segments(x0=coordinates$x[1], y0=min(coordinates$y), x1=coordinates$x[2], y1=min(coordinates$y), col=color, lwd=3)

  if(saveResult==TRUE){
    dev.copy(pdf,paste0(gateName, "_gate.pdf"))
    dev.off()
  }

  if(missing(outputDataset)){
    outputDataset <- dataSet
  }
  gatedData <- outputDataset[dataSet[,variable]>min(coordinates$x) & dataSet[,variable]<max(coordinates$x),]

  gatedDataPlusRows <- data.frame(gatedData, which(dataSet[,variable]>min(coordinates$x) & dataSet[,variable]<max(coordinates$x)))
  colnames(gatedDataPlusRows)[ncol(gatedDataPlusRows)] <- paste0("Rows_in_", deparse(substitute(dataSet)))

  if(makeDir==TRUE){
    newDir <- paste0("./", gateName)
    dir.create(newDir)
    setwd(newDir)
  }

  return(gatedDataPlusRows)

}
