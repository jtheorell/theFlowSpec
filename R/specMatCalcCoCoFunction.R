specMatCalcCoCoFunction <- function(sampleX){

  medianList <- list()
  for(i in 1:ncol(sampleX)){
    focusColumn <- sampleX[,i]
    histData <- hist(focusColumn, breaks=100, plot=FALSE)
    #Here, the bin where the most events are present is selected
    peakLow <- histData$breaks[which.max(histData$counts)]
    peakHigh <- histData$breaks[which.max(histData$counts)+1]
    #Now, the dataset in this bin is selected, and the median found for all the columns
    focusPeak <- sampleX[sampleX[,i]>peakLow & sampleX[,i]<peakHigh,]

    #And now, the median in the defined range is identified
    medianList[[i]] <- apply(focusPeak, 2, median)

  }

  #Here, all the results are compiled
  medianDf <- do.call("rbind", medianList)

  #And now the median from all observed peaks is identified
  medianResult <- apply(medianDf, 2, median)

}
