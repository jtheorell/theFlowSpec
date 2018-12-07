minDensGateCo <- function(Donor, markerFocus, Da, Turns, topTwo, minTurn, selectPos, i, marker, returnThreshold){

      pdf(paste0("Donor_", i, "_", marker, "_gate_minDens.pdf"))
      plot(Da, xlab="", ylab="", main="")
      points(Da$x[topTwo], Da$y[topTwo], pch=16, cex=1, col="red")
      points(Da$x[minTurn], Da$y[minTurn], pch=16, cex=2, col="blue")
      dev.off()
      gateVal <- Da$x[minTurn]
      if(selectPos==TRUE){
        if(returnThreshold==TRUE){
          return(list(Donor[markerFocus>gateVal,], "threshVal"=gateVal))
        } else {
          return(Donor[markerFocus>gateVal,])
        }

      } else {
        if(returnThreshold==TRUE){
          return(list(Donor[markerFocus<=gateVal,], "threshVal"=gateVal))
        } else {
          return(Donor[markerFocus<=gateVal,])
        }

      }

}
