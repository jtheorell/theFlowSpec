
madDecision <- function(data, marker, gateVal){
  lowData <- data[data[,marker]<gateVal,]
  highData <- data[data[,marker]>=gateVal,]

  if(class(lowData)=="data.frame" && class(highData)=="data.frame"){
    if(nrow(lowData)>50 & nrow(highData)>50){
      lowCenters <- apply(lowData, 2, median)
      highCenters <- apply(highData, 2, median)

      diffCenters <- abs(lowCenters-highCenters)
      diffCentSquared <- diffCenters^2
      return(sum(diffCentSquared))
    }
    else {
      return(0)
    }
  }
   else {
    return(0)
  }

}
