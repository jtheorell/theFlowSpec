#' @export gating1D
gating1D <- function(dataSet, variable, outputDataset, gateName=variable, color="red", saveResult=TRUE){
  hist(dataSet[,variable], breaks=200, xlab=variable, main=gateName)
  coordinates <- locator(type="n")

  segments(x0=coordinates$x[1], y0=min(coordinates$y), x1=coordinates$x[2], y1=min(coordinates$y), col=color, lwd=3)

  if(saveResult==TRUE){
    dev.copy(png,paste0(gateName, "gate.png"))
    dev.off()
  }

  if(missing(outputDataset)){
    outputDataset <- dataSet
  }
  gatedData <- outputDataset[dataSet[,variable]>min(coordinates$x) & dataSet[,variable]<max(coordinates$x),]

}
