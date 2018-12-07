# @importFrom graphics C_locator
#' @importFrom splancs inpip
#' @export gating2D
gating2D <- function(dataSet, xVar, yVar, gate=TRUE, gateName="default", color="red", sampleSize=10000, dotSize="default", saveResult=TRUE, makeDir=TRUE){

  if(nrow(dataSet)>sampleSize){
    dataSampleRows <- sample(1:nrow(dataSet), size=sampleSize)
    xDataSample <- dataSet[dataSampleRows,xVar]
    yDataSample <- dataSet[dataSampleRows,yVar]
  }
  else {
    xDataSample <- dataSet[,xVar]
    yDataSample <- dataSet[,yVar]
  }
  if(dotSize=="default"){
    dotSize <- 100/sqrt(length(xDataSample))
  }
  if(gateName=="default"){
    gateName <- paste0(xVar, "_vs_", yVar)
  }


  oneVsAllPlotCoFunction(xVar=xDataSample, yVar=yDataSample, color=color, dotSize=dotSize, xlab=xVar, ylab=yVar)

    coordinates <- locator(type="n")
    polygon(x=coordinates$x, y=coordinates$y, col="#0000FF55", border=color, lwd=3)

  if(saveResult==TRUE){
    dev.copy(png,paste0(gateName, "_gate.png"))
    dev.off()
  }

  gateSubset <- inpip(data.frame(x=dataSet[,xVar], y=dataSet[,yVar]), coordinates)

  gatedData <- dataSet[gateSubset,]

  gatedDataPlusRows <- data.frame(gatedData, gateSubset)
  colnames(gatedDataPlusRows)[ncol(gatedDataPlusRows)] <- paste0("Rows_in_", deparse(substitute(dataSet)))

  if(makeDir==TRUE){
    newDir <- paste0("./", gateName)
    dir.create(newDir)
    setwd(newDir)
  }
  return(gatedDataPlusRows)

}
