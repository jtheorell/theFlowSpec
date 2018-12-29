#' @export bestGateAngle
bestGateAngle <- function(focusVar, Var2, focusName=deparse(substitute(focusVar)), name2=deparse(substitute(focusVar)), minFocusVarContrib=0.5){

  testData <- cbind(focusVar, Var2)

  focusVarVal <- c(seq(-1, -minFocusVarContrib, length.out=51), seq(minFocusVarContrib, 1, length.out=51))

  var2Val <- 1-abs(focusVarVal)

  focusList <- sapply(1:length(focusVarVal), function(x) testData[,1]*focusVarVal[x]+testData[,2]*var2Val[x], simplify=FALSE)
  resultList <- list()
  for(i in 1:length(focusList)){
    Da = density(focusList[[i]], adjust=2)

    DeltaY = diff(Da$y)
    Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
    #Identify the two highest turnpoints
    topTwo <- Turns[order(Da$y[Turns], decreasing=TRUE)][1:2]
    #Now, identify the lowest turnpoint between them

    robRange <- quantile(focusList[[i]], c(0.1, 0.9))
    tenth <- unname((robRange[2]-robRange[1])*0.1)
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

      volumePeak1 <- length(focusList[[i]][which(focusList[[i]]>Da$x[startPeak1] & focusList[[i]]<Da$x[stopPeak1])])
      volumePeak2 <- length(focusList[[i]][which(focusList[[i]]>Da$x[startPeak2] & focusList[[i]]<Da$x[stopPeak2])])

      if(volumePeak1>0 && volumePeak2>0 && volumePeak2/volumePeak1>0.1){
        resultList[[i]] <- list("lowFraction"=Da$y[minTurn]/Da$y[topTwo[2]], "Da"=Da, "minTurn"=minTurn, "topTwo"=topTwo, "focusVarVal"=focusVarVal[i], "var2Val"=var2Val[i])
      } else {
        resultList[[i]] <- 1
      }

    } else {
      resultList[[i]] <- 1
    }
  }
  lowFraction <- unlist(lapply(resultList, "[[", 1))

  if(min(lowFraction)<1){
    bestOption <- resultList[[which.min(lowFraction)[1]]]
    Da <- bestOption$Da
    topTwo <- bestOption$topTwo
    minTurn <- bestOption$minTurn
    pdf(paste0("Best_Angle_", focusName, "_vs_", name2, ".pdf"))
    plot(Da, xlab="", ylab="", main="")
    points(Da$x[topTwo], Da$y[topTwo], pch=16, cex=1, col="red")
    points(Da$x[minTurn], Da$y[minTurn], pch=16, cex=2, col="blue")
    dev.off()
    varLoadings <- c(bestOption$focusVarVal, bestOption$var2Val)
    names(varLoadings) <- c(focusName, name2)

    return(list("bestValue"=min(lowFraction), varLoadings))
  }
  else {
    return(list("bestValue"=1))
  }

}
