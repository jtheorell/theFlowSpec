#' @export conservGate
conservGate <- function(DonorList, marker, lowPeak=TRUE, madVal=3, turnOff=FALSE, adjust=2, volRatio=0.05, createDir=TRUE){
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
        twoPeaks <- FALSE
      }
    } else {
      peakCenter <- topTwo[1]
      twoPeaks <- FALSE
    }

    #Now, teh most conservative side of the peak is identified
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

    #Here, the tightest one is selected. It will in most instances be the lower half.
    optimPeak <- if(sd(lowHalfPeak)<sd(highHalfPeak)){lowHalfPeak} else {highHalfPeak}

    optimSide <- if(sd(lowHalfPeak)<sd(highHalfPeak)){1} else {2}

    gateVal <- peakValue+mad(optimPeak)*madVal

   #Now, it is determined if the gamma data should be used, or if indeed a density-based gate is more useful
    if(twoPeaks==TRUE && turnOff==FALSE){
      if(lowPeak==TRUE && madVal>0 && length(which(Da$x[Turns]>peakValue & Da$x[Turns]<gateVal))>0){
        print("A turn comes before the gamma level, so it is used instead for gating, for more correct separation")
        gateVal <- Da$x[Turns[which(Da$x[Turns]>peakValue & Da$x[Turns]<gateVal)][1]]
        graphName <- paste0("Donor_", i, "_", marker, "_turnGate.pdf")
      } else if(lowPeak==FALSE && madVal<0 && length(which(Da$x[Turns]<peakValue & Da$x[Turns]>gateVal))>0){
        print("A turn comes before the gamma level, so it is used instead for gating, for more correct separation")
        gateVal <- Da$x[Turns[which(Da$x[Turns]<peakValue & Da$x[Turns]>gateVal)][1]]
        graphName <- paste0("Donor_", i, "_", marker, "_turnGate.pdf")
      } else {
        graphName <- paste0("Donor_", i, "_", marker, "_gammaGate.pdf")
      }

    }
    pdf(graphName)
    plot(Da)
    polygon(x=Da$x, y=Da$y, col="#0000FF55", border="blue", lwd=2)
    abline(v=gateVal, col="black", lwd=2, lty=2)
    dev.off()
    resultList[[i]] <- gateVal
    print(paste0("Done with individual ", i))
  }
  if(createDir==TRUE){
    setwd(oldDir)
  }

  return(resultList)
}
