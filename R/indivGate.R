#' @export indivGate
indivGate <- function(DonorList, marker, madAboveMed=TRUE, selectPos=FALSE, adjust=2, returnThreshold=TRUE, createDir=TRUE){
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
    robRange <- quantile(markerFocus, c(0.1, 0.9))
    tenth <- abs(unname((robRange[2]-robRange[1])*0.1))
    #Identify the two highest turnpoints
    topTwo <- Turns[order(Da$y[Turns], decreasing=TRUE)][1:2]

    #Here, a decision is made regarding in which direction to go.
    #If more than 10% of the events are contained in a second peak
    #(defined as a peak with its maximum separated by more than 10%
    #of the range from the first to the 99th percentile from the
    #highest peak), minDensGate is used.
    #Othewise, madOrTurnGate is used
    if(is.na(topTwo[2])==FALSE && abs(Da$x[topTwo][1]-Da$x[topTwo][2])>tenth){
      #Now, identify the lowest turnpoint between the topTwo
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

      if(volumePeak1>0 && volumePeak2>0 && volumePeak2/volumePeak1>0.05){
        print("minDens")
        resultList[[i]] <- theFlowSpec:::minDensGateCo(DonorList[[i]], markerFocus, Da, Turns, topTwo, minTurn, selectPos, i, marker, returnThreshold)
      }
      else {
        resultList[[i]] <- theFlowSpec:::madOrTurnGateCo(DonorList[[i]], markerFocus, Da, Turns, i, marker, madAboveMed, selectPos, returnThreshold)
      }
    } else {
      resultList[[i]] <- theFlowSpec:::madOrTurnGateCo(DonorList[[i]], markerFocus, Da, Turns, i, marker, madAboveMed, selectPos, returnThreshold)
    }
  }
  if(createDir==TRUE){
    setwd(oldDir)
  }

  return(resultList)
}
