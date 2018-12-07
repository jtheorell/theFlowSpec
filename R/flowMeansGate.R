#' @export flowMeansGate
flowMeansGate <- kMeansGate <- function(DonorList, primMarker, secondMarkers, posNegOrBoth="both", primWeight=2){
  #First, we enhance the signal for the channel we want to gate on
  resultList <- list()
  for(i in 1:length(DonorList)){
    focusDataSet <- cbind("Focus"=DonorList[[i]][,primMarker]*primWeight, DonorList[[i]][,secondMarkers])
    kmeansResult <- flowMeans(focusDataSet, MaxN = 2, NumC = 2, Standardize=FALSE)

    dispData <- data.frame("Cluster"=as.character(kmeansResult@Label), primMarker=DonorList[[i]][,primMarker])
    ggplot(dispData, aes(primMarker, fill = Cluster, colour = Cluster)) +
      geom_density(alpha = 0.1) + xlab(paste0(primMarker, " k means gate"))
    ggsave(filename=paste0("Donor_", i, "_", primMarker, "_kmeans_gate.pdf"))

    Clust1 <- DonorList[[i]][kmeansResult@Label==1,]
    Clust2 <- DonorList[[i]][kmeansResult@Label==2,]

    posClust <- which.max(c(median(Clust1[,primMarker]), median(Clust2[,primMarker])))

    if(posNegOrBoth=="pos"){
      resultList[[i]] <- DonorList[[i]][kmeansResult@Label==posClust,]
    } else if(posNegOrBoth=="neg"){
      resultList[[i]] <- DonorList[[i]][kmeansResult@Label!=posClust,]
    } else {

      localList <- list(DonorList[[i]][kmeansResult@Label==posClust,], DonorList[[i]][kmeansResult@Label!=posClust,])
      names(localList) <- c(paste0(primMarker, "_pos"), paste0(primMarker, "_neg"))
      resultList[[i]] <- localList
    }


  }
  return(resultList)
}

