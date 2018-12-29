#' @importFrom DepecheR dScale
#' @importFrom FNN knnx.index
#' @export eucSep
eucSep <- function(inDatDonList, gateMarker, gateList, euclidMarkers, gateWeight=2, createDir=TRUE){
  if(createDir==TRUE){
    oldDir <- getwd()
    newDir <- paste0("./Individual_euclid_separation_", gateMarker)
    dir.create(newDir)
    setwd(newDir)
  }
  resultList <- list()
  for(i in 1:length(inDatDonList)){
    euclidFocus <- inDatDonList[[i]][,euclidMarkers]
    if(nrow(euclidFocus[euclidFocus[gateMarker]>=gateList[[i]],])>50){
      euclidScaled <- dScale(euclidFocus)
      euclidScaled[,gateMarker] <- euclidScaled[,gateMarker]*gateWeight
      euclidFocus$GateVal <- "neg"
      euclidFocus$GateVal[euclidFocus[gateMarker]>=gateList[[i]]] <- "pos"
      lowCent <- apply(euclidScaled[euclidFocus$GateVal=="neg",], 2, mean)
      highCent <- apply(euclidScaled[euclidFocus$GateVal=="pos",], 2, mean)
      centMat <- rbind(lowCent, highCent)
      neighbors <- knnx.index(centMat, euclidScaled, k=1)[,1]
      
    } else {
      neighbors <- rep(1, length.out=nrow(euclidFocus))
      neighbors[euclidFocus[,gateMarker]>=gateList[[i]]] <- 2
    }
    
    fullPop <- density(euclidFocus[, gateMarker],from=min(euclidFocus[,gateMarker]), to=max(euclidFocus[,gateMarker]))
    fullPopScaled <- fullPop
    fullPopScaled$y <- fullPop$y*length(euclidFocus[, gateMarker])
    
    negPop <- density(euclidFocus[neighbors==1, gateMarker],from=min(euclidFocus[,gateMarker]), to=max(euclidFocus[,gateMarker]))
    negPopScaled <- negPop
    negPopScaled$y <- negPop$y*length(euclidFocus[neighbors==1, gateMarker])
    
    posPop <- density(euclidFocus[neighbors==2, gateMarker],from=min(euclidFocus[,gateMarker]), to=max(euclidFocus[,gateMarker]))
    posPopScaled <- posPop
    posPopScaled$y <- posPop$y*length(euclidFocus[neighbors==2, gateMarker])
    
    pdf(paste0("Donor_", names(inDatDonList)[i], "_", gateMarker, "_Euclid_sep.pdf"))
    plot(fullPopScaled)
    polygon(x=fullPopScaled$x, y=fullPopScaled$y, col="#99999922", border="black", lwd=2, lty=2)
    polygon(x=negPopScaled$x, y=negPopScaled$y, col="#0000FF55", border="blue", lwd=2)
    polygon(x=posPopScaled$x, y=posPopScaled$y, col="#FF000055", border="red", lwd=2)
    abline(v=gateList[[i]], col="black", lwd=2, lty=2)
    dev.off()
    
    resultList[[i]] <- neighbors
    print(paste0("Done with individual ", names(inDatDonList)[i]))
  }
  if(createDir==TRUE){
    setwd(oldDir)
  }
  return(resultList)
}
