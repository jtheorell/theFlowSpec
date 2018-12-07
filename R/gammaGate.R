#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
# @importFrom GA ga
#' @export gammaGate
gammaGate <- function(DonorList, marker, lowPeak=TRUE, gammaQuantile=0.95, turnOff=FALSE, adjust=2, volRatio=0.05, createDir=TRUE){
  if(createDir==TRUE){
    oldDir <- getwd()
    newDir <- paste0("./Individual_gates_", marker)
    dir.create(newDir)
    setwd(newDir)
  }

  resultList <- list()
  for(i in 1:length(DonorList)){
  markerFocus <- DonorList[[i]][,marker]
  Da = density(markerFocus, adjust=adjust)
  DeltaY = diff(Da$y)
  Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
  #Identify the two highest turnpoints
  topTwo <- Turns[order(Da$y[Turns], decreasing=TRUE)][1:2]

  #And here, it is made clear if it is meaningful to talk about two peaks, or if it is indeed only one.
  #First, it is determined if there are two peaks whatsoever, and if so, by how much these are separated. Currently, they need to be separated by at least 10% of the total range from the first to th 99th percentile to count as two peaks.
  robRange <- quantile(markerFocus, c(0.1, 0.9))
  tenth <- abs(unname((robRange[2]-robRange[1])*0.1))
  if(is.na(topTwo[2])==FALSE && abs(Da$x[topTwo][1]-Da$x[topTwo][2])>tenth){
    #Now, if two peaks are present, identify the lowest turnpoint between the topTwo
    betweenTurns <- Turns[Turns<max(topTwo) & Turns>min(topTwo)]
    minTurn <- betweenTurns[which.min(Da$y[betweenTurns])]
    if(topTwo[1]>topTwo[2]){
      startPeak1 <- minTurn
      stopPeak1 <- Turns[length(Turns)]
      startPeak2 <- Turns[1]
      stopPeak2 <- minTurn
    } else {
      startPeak1 <- Turns[1]
      stopPeak1 <- minTurn
      startPeak2 <- minTurn
      stopPeak2 <- Turns[length(Turns)]
    }

    volumePeak1 <- length(markerFocus[which(markerFocus>Da$x[startPeak1] & markerFocus<Da$x[stopPeak1])])
    volumePeak2 <- length(markerFocus[which(markerFocus>Da$x[startPeak2] & markerFocus<Da$x[stopPeak2])])

    #And here, it is determined if the second peak counts, which it only does if it has a volume that is at least a set percentage of the largest one. Default being 5%.
    if(volumePeak1>0 && volumePeak2>0 && volumePeak2/volumePeak1>volRatio){
      twoPeaks <- TRUE
      if(lowPeak==TRUE){
        peakCenter <- min(topTwo)
      } else {
        peakCenter <- max(topTwo)
      }
    } else {
      peakCenter <- topTwo[1]
          }
  } else {
    peakCenter <- topTwo[1]
  }

  #Now, it is decided whether the peak can be considered symmetrical, and thus, that the full peak, from -3 MAD to +3 MAD (of the more conservative side) can be included, or if a syntehtic peak, only including the more conservative side, is better.
  peakValue <- Da$x[peakCenter]
  if(Turns[1]==peakCenter){
    lowBorder <- Da$x[1]
  } else {
    lowBorder <- Da$x[Turns[which(Turns==peakCenter)-1]]
  }
  if(Turns[length(Turns)]==peakCenter){
    highBorder <- Da$x[length(Da$x)]
  } else {
    highBorder <- Da$x[Turns[which(Turns==peakCenter)+1]]
  }

  #Now, the two "different" peaks are created, representing the mirrored  lower halfof the peak and vice versa
  lowHalfCent <- markerFocus[markerFocus>=lowBorder & markerFocus<peakValue]-peakValue
  lowHalfPeak <- c(lowHalfCent, lowHalfCent*-1)

  highHalfCent <- markerFocus[markerFocus>=peakValue & markerFocus<highBorder]-peakValue
  highHalfPeak <- c(highHalfCent, highHalfCent*-1)

  #Here, the optimal one is selected. It will in most instances be the lower half.
  conservPeak <- if(sd(lowHalfPeak)<sd(highHalfPeak)){lowHalfPeak} else {highHalfPeak}

  if(conservPeak==lowHalfPeak){
    threeMAD <- peakValue+3*sd(conservPeak)
    madPoint <- as.numeric(FNN::knnx.index(Da$x, threeMAD, k=1))
    if(Da$y[madPoint]/Da$y[peakCenter]<0){
      print("Testing, not working")
    }
  }
  optimPeak <-
  optimPeak <- if(sd(lowHalfPeak)<sd(highHalfPeak)){lowHalfPeak} else {highHalfPeak}

  optimSide <- if(sd(lowHalfPeak)<sd(highHalfPeak)){1} else {2}

  #And now, the optimal shape of the gamma distribution is identified
  if(length(optimPeak)>10000){
    optimPeakSample <- sample(optimPeak, 10000)
  } else {
    optimPeakSample <- optimPeak
  }

  optimPeakSample <- scale(optimPeakSample)
  optimPeakDf <- data.frame("Data"=optimPeakSample[,1], "Group"=rep("data", times=length(optimPeakSample[,1])))

  #bestShape <- ga(type = "real-valued", fitness = gammaEstOpt, n=nrow(optimPeakDf), optimSide=optimSide, optimPeakDf=optimPeakDf, adjust=adjust, lower = 1, upper = 100, maxiter=50, run=10, parallel =FALSE)

  nCores <- parallel:::detectCores() - 1

  cl <-  parallel:::makeCluster(nCores, type = "SOCK")
  registerDoSNOW(cl)
  allShapes <- unlist(foreach(i=1:100, .packages="theFlowSpec") %dopar% gammaEstOpt(shape=i, n=nrow(optimPeakDf), optimSide=optimSide, optimPeakDf=optimPeakDf, adjust=adjust))
  parallel:::stopCluster(cl)

  #Now, the optimal solution is generated
  gammaBest <- rgamma(10000, shape=c(1:100)[which.max(allShapes)])
  #Here, the variance is adjusted to the variance of the real population
  DaGammaRaw <- density(gammaBest)
  #Now, define the peak
  DeltaYGamma = diff(DaGammaRaw$y)
  TurnsGamma = which(DeltaYGamma[-1] * DeltaYGamma[-length(DeltaYGamma)] < 0) + 1
  #Identify the two highest turnpoints
  gammaPeak <- TurnsGamma[which.max(DaGammaRaw$y[TurnsGamma])]

  #Now, depending on which side of the original peak that is reliable, define an interval to which it is possible to scale
  if(optimSide==1){
    scaleRangeData <- peakValue-quantile((optimPeak+peakValue), 0.05)
    scaleRangeGamma <- DaGammaRaw$x[gammaPeak]-quantile(gammaBest, 0.05)
  } else {
    scaleRangeData <- quantile((optimPeak+peakValue), 0.95)-peakValue
    scaleRangeGamma <- quantile(gammaBest, 0.95)-DaGammaRaw$x[gammaPeak]
  }
  gammaBestScaled <- (gammaBest/scaleRangeGamma)*scaleRangeData

  #gammaBestScaled <- (gammaBest/sd(gammaBest))*sd(optimPeak)
  #scaledPeak <- (DaGammaRaw$x[gammaPeak]/sd(gammaBest))*sd(optimPeak)

  scaledPeak <- (DaGammaRaw$x[gammaPeak]/scaleRangeGamma)*scaleRangeData
  gammaBestCentered <- (gammaBestScaled-scaledPeak)+peakValue
  DaGamma <- density(gammaBestCentered, adjust=adjust)

  #Here, the height of the density peak is made equal between the two plots
  DaScaled <- Da
  DaScaled$y <- DaScaled$y/DaScaled$y[peakCenter]
  DaGammaScaled <- DaGamma
  DaGammaScaled$y <- DaGamma$y/max(DaGamma$y)

  #Now, it is determined if the gamma data shoulw be used, or if indeed a density-based gate is more useful
  gammaGate <- quantile(gammaBestCentered, gammaQuantile)
  if(twoPeaks==TRUE && turnOff==FALSE){
    if(lowPeak==TRUE && gammaQuantile>0.5 && length(which(Da$x[Turns]>peakValue & Da$x[Turns]<gammaGate))>0){
      print("A turn comes before the gamma level, so it is used instead for gating, for more correct separation")
      gammaGate <- Da$x[Turns[which(Da$x[Turns]>peakValue & Da$x[Turns]<gammaGate)][1]]
      graphName <- paste0("Donor_", i, "_", marker, "_turnGate.pdf")
      } else if(lowPeak==FALSE && gammaQuantile<0.5 && length(which(Da$x[Turns]<peakValue & Da$x[Turns]>gammaGate))>0){
      print("A turn comes before the gamma level, so it is used instead for gating, for more correct separation")
      gammaGate <- Da$x[Turns[which(Da$x[Turns]<peakValue & Da$x[Turns]>gammaGate)][1]]
      graphName <- paste0("Donor_", i, "_", marker, "_turnGate.pdf")
    } else {
      graphName <- paste0("Donor_", i, "_", marker, "_gammaGate.pdf")
    }

  }
  pdf(graphName)
  plot(DaScaled)
  polygon(x=DaGammaScaled$x, y=DaGammaScaled$y, col="#FF000055", border="red", lwd=2, lty=3)
  polygon(x=DaScaled$x, y=DaScaled$y, col="#0000FF55", border="blue", lwd=2)
  abline(v=gammaGate, col="black", lwd=2, lty=2)
  dev.off()
  resultList[[i]] <- quantile(gammaBestCentered, gammaQuantile)
  print(paste0("Done with individual ", i))
  }
  if(createDir==TRUE){
    setwd(oldDir)
  }

return(resultList)
}
