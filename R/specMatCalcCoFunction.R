specMatCalcCoFunction <- function(sampleX){
  peakValues <- vector()
  for(i in 1:ncol(sampleX)){
    x <- sampleX[,i]
    histData <- hist(x, breaks=100, plot=FALSE)
    #Here, the bin where the most events are present is selected
    peakLow <- histData$breaks[which.max(histData$counts)]
    peakHigh <- histData$breaks[which.max(histData$counts)+1]
    #And now, the median in the defined range is identified
    peakValues[i] <- median(x[which(x>peakLow & x<peakHigh)])

  }

return(peakValues)
}
