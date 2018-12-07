#' @importFrom FNN knnx.index
#' @export madOrTurnGateCo
madOrTurnGateCo <- function(Donor, markerFocus, Da, Turns, i, marker, madAboveMed, selectPos, returnThreshold){

    if(madAboveMed==FALSE){
      threeMAD <- median(markerFocus)-3*mad(markerFocus)
      twoMAD <- median(markerFocus)-2*mad(markerFocus)
    } else {
      threeMAD <- median(markerFocus)+3*mad(markerFocus)
      twoMAD <- median(markerFocus)+2*mad(markerFocus)
    }
    madPoint <- as.numeric(FNN::knnx.index(Da$x, threeMAD, k=1))

    if((madAboveMed==FALSE && length(which(Da$x[Turns]<twoMAD & Da$x[Turns]>threeMAD)>0)) || (madAboveMed==TRUE && length(which(Da$x[Turns]>twoMAD & Da$x[Turns]<threeMAD))>0)){
      print("Turn")
      if(madAboveMed==FALSE){
        gatePositions <- Turns[which(Da$x[Turns]<twoMAD & Da$x[Turns]>threeMAD)]
        gatePos <- gatePositions[length(gatePositions)]
      } else {
        gatePositions <- Turns[which(Da$x[Turns]>twoMAD & Da$x[Turns]<threeMAD)]
        gatePos <- gatePositions[1]
      }

      pdf(paste0("Donor_", i, "_", marker, "_gate_devTurn.pdf"))
      plot(Da, xlab="", ylab="", main="")
      points(Da$x[madPoint], Da$y[madPoint], pch=16, cex=1, col="blue")
      points(Da$x[gatePos], Da$y[gatePos], pch=16, cex=2, col="red")
      dev.off()
      if(selectPos==TRUE){
        if(returnThreshold==TRUE){
          return(list(Donor[markerFocus>Da$x[gatePos],], "threshVal"=gatePos))
        } else {
          return(Donor[markerFocus>Da$x[gatePos],])
        }

      } else {
        if(returnThreshold==TRUE){
          return(list(Donor[markerFocus<Da$x[gatePos],], "threshVal"=gatePos))
        } else {
          return(Donor[markerFocus<Da$x[gatePos],])
        }

      }

    } else {
      print("MAD")
      pdf(paste0("Donor_", i, "_", marker, "_gate_mad.pdf"))
      plot(Da, xlab="", ylab="", main="")
      points(Da$x[madPoint], Da$y[madPoint], pch=16, cex=2, col="blue")
      dev.off()
      if(selectPos==TRUE){
        if(returnThreshold==TRUE){
          return(list(Donor[markerFocus>threeMAD,], "threshVal"=threeMAD))
        } else {
          return(Donor[markerFocus>threeMAD,])
        }

      } else {
        if(returnThreshold==TRUE){
          return(list(Donor[markerFocus<=threeMAD,], "threshVal"=threeMAD))
        } else {
          return(Donor[markerFocus<=threeMAD,])
        }

      }
    }

}
