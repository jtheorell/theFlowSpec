#' @importFrom flowCore flowFrame
#' @importFrom flowVS transFlowVS
# transformList logicleTransform transform
# @importFrom flowTrans flowTrans
#' @export transformData
transformData <- function(rawData, coFactors=rep(400, ncol(rawData))){
#  if(nrow(rawData)>10000){
#    rawDataSubset <- rawData[sample(1:nrow(rawData), size=10000),]
#  } else {
#    rawDataSubset <- rawData
#  }
#  rawSubsetFrame <- flowFrame(as.matrix(rawDataSubset))
#  FCSTransTransform(transformationId = "defaultFCSTransTransform", channelrange = 2^18, channeldecade = 4.5, range = 4096, cutoff = -111, w = NULL, rescale = TRUE)
#  transformParam <- logicleTransform(transformationId="defaultLogicleTransform", w = w, t = t, m = m, a = a)
#  transedData <- transform(flowFrame(as.matrix(rawData)), transformList(colnames(rawData), transformParam))
  rawSet <- flowSet(flowFrame(as.matrix(rawData)))
  transedData <- transFlowVS(rawSet, channels=colnames(rawData), cofactors=coFactors)
  afterExprs <- exprs(transedData[[1]])
}
