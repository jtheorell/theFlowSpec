#' @importFrom DepecheR dScale
#' @importFrom FNN knnx.index
#' @export turnDerGateCoFunction
turnDerGateCoFunction <- function(euclidFocus, gateMarker, gateVal, gateWeight,
                                  graphName, adjust, n){
  #First, the area closest to the gate is exluded from the populations. This is
  #done by excluding not half of the events, but rather half of the x-axis area
  #that the population takes up.
  allNeg <- euclidFocus[euclidFocus[,gateMarker]<gateVal,]
  allPos <- euclidFocus[euclidFocus[,gateMarker]>=gateVal,]
  clearNegRange <- quantile(allNeg[,gateMarker], c(0.01, 0.99))
  clearPosRange <- quantile(allPos[,gateMarker], c(0.01, 0.99))
  clearNegBorder <- clearNegRange[1]+((clearNegRange[2]-clearNegRange[1])/4)*3
  clearPosBorder <- clearPosRange[1]+((clearPosRange[2]-clearPosRange[1])/4)
  euclidFocus$GateResult <- "Unclear"
  euclidFocus$GateResult[euclidFocus[,gateMarker]<clearNegBorder] <- "neg"
  euclidFocus$GateResult[euclidFocus[,gateMarker]>clearPosBorder] <- "pos"

  if(nrow(euclidFocus[euclidFocus$GateResult=="neg",])>50 &&
     nrow(euclidFocus[euclidFocus$GateResult=="pos",])>50){
    euclidScaled <- dScale(euclidFocus[-ncol(euclidFocus)])
    euclidScaled[,gateMarker] <- euclidScaled[,gateMarker]*gateWeight
    negPop <- euclidScaled[euclidFocus$GateResult=="neg",-ncol(euclidFocus)]
    posPop <- euclidScaled[euclidFocus$GateResult=="pos",-ncol(euclidFocus)]
    lowCent <- apply(negPop, 2, mean)
    highCent <- apply(posPop, 2, mean)
    centMat <- rbind(lowCent, highCent)
    neighbors <- knnx.index(centMat, euclidScaled, k=1)[,1]

  } else {
    neighbors <- rep(1, length.out=nrow(euclidFocus))
    neighbors[euclidFocus[,gateMarker]>=gateVal] <- 2
  }

  negPop <- density(euclidFocus[neighbors==1, gateMarker],
                    from=min(euclidFocus[,gateMarker]),
                    to=max(euclidFocus[,gateMarker]), adjust=adjust, n=n)
  #Here, all points outside the data range are excluded
  minNeg <- min(euclidFocus[neighbors==1, gateMarker])
  maxNeg <- max(euclidFocus[neighbors==1, gateMarker])
  negRangeX <- negPop$x[which(negPop$x>=minNeg & negPop$x<=maxNeg)]
  negY <- negPop$y[which(negPop$x>=minNeg & negPop$x<=maxNeg)]
  #Now, 0-values on the edges are added, to make the graph look right
  negRangeXFull <- c(negRangeX[1], negRangeX, negRangeX[length(negRangeX)])
  negYFull <- c(0, negY, 0)
  #Now, a constant is defined, that makes the size of the population accurate
  negConst <- length(euclidFocus[neighbors==1, gateMarker])/sum(negYFull)
  negYScaled <- negYFull*negConst

  #(negYFull/max(negYFull))*length(euclidFocus[euclidFocus$GateVal=="neg", gateMarker])

  posPop <- density(euclidFocus[neighbors==2, gateMarker],
                    from=min(euclidFocus[,gateMarker]),
                    to=max(euclidFocus[,gateMarker]), adjust=adjust, n=n)
  #Here, all points outside the data range are excluded
  minpos <- min(euclidFocus[neighbors==2, gateMarker])
  maxpos <- max(euclidFocus[neighbors==2, gateMarker])
  posRangeX <- posPop$x[which(posPop$x>=minpos & posPop$x<=maxpos)]
  posY <- posPop$y[which(posPop$x>=minpos & posPop$x<=maxpos)]
  #Now, 0-values on the edges are added, to make the graph look right
  posRangeXFull <- c(posRangeX[1], posRangeX, posRangeX[length(posRangeX)])
  posYFull <- c(0, posY, 0)
  posConst <- length(euclidFocus[neighbors==2, gateMarker])/sum(posYFull)
  posYScaled <- posYFull*posConst

  fullPop <- density(euclidFocus[, gateMarker],
                     from=min(euclidFocus[,gateMarker]),
                     to=max(euclidFocus[,gateMarker]), adjust=adjust, n=n)
  #Here, all points outside the data range are excluded
  minFull <- min(euclidFocus[, gateMarker])
  maxFull <- max(euclidFocus[, gateMarker])
  fullRangeX <- fullPop$x[which(fullPop$x>=minFull & fullPop$x<=maxFull)]
  fullY <- fullPop$y[which(fullPop$x>=minFull & fullPop$x<=maxFull)]
  #Now, 0-values on the edges are added, to make the graph look right
  fullRangeXFull <- c(fullRangeX[1], fullRangeX, fullRangeX[length(fullRangeX)])
  fullYFull <- c(0, fullY, 0)
  fullConst <- length(euclidFocus[, gateMarker])/sum(fullYFull)
  fullYScaled <- fullYFull*fullConst

  pdf(graphName)
  plot(c(0, max(c(negRangeXFull, posRangeXFull, fullRangeXFull))),
       c(0, max(c(negYScaled, posYScaled, fullYScaled))), main="", col="white")
  polygon(x=fullRangeXFull, y=fullYScaled, col="white", border="black", lty=2,
          lwd=2)
  polygon(x=negRangeXFull, y=negYScaled, col="#0000FF55", border="blue", lwd=2)
  polygon(x=posRangeXFull, y=posYScaled, col="#FF000055", border="red", lwd=2)
  abline(v=gateVal, col="black", lwd=2, lty=2)
  dev.off()

  return(neighbors)

}
