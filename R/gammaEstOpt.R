#' @export gammaEstOpt
gammaEstOpt <- function(shape, n, optimSide, optimPeakDf, adjust){
before <- Sys.time()
gammaDist <- rgamma(n, shape)

#Now, half of the distribution is selected, depending on the optimSide
before <- Sys.time()
Da = density(gammaDist, adjust=adjust)
DeltaY = diff(Da$y)
Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
top <- Turns[order(Da$y[Turns], decreasing=TRUE)][1]
after <- Sys.time()-before
compHalf <- if(optimSide==1){gammaDist[gammaDist<=Da$x[top]]-Da$x[top]} else {gammaDist[gammaDist>Da$x[top]]-Da$x[top]}
compPeak <- c(compHalf, compHalf*-1)
#Now, peak is scaled to unit variance
compPeakScaled <- scale(compPeak)

#And here comes the comparison
compPeakDf <- data.frame("Data"=compPeakScaled[,1], "Group"=rep("gamma", times=length(compPeakScaled[,1])))
allData <- rbind(compPeakDf, optimPeakDf)
allDataOrdered <- allData[order(allData[,1]),]
allDataList <- split(allDataOrdered[,2], ceiling(seq_along(allDataOrdered[,2])/100))
result <- 1-sum(sapply(allDataList, function(x) abs(length(which(x=="gamma"))-length(which(x=="data")))^2))

after <- Sys.time()-before

}
