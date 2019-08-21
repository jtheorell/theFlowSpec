#' Function used to generate automated gates on large sets of samples
#'
#'
#' In addition to performing optimized gating based on either the most suitable turning point between two populaitons in the case of a bimodal distritbution,
#' or finding a specific turning point in the third derivative in the case of a single peak distribution, the function goes on to identify which cells that are most
#' similar to the cluster center of each of the two populations created by the gate. It does this not just by dividing along the gate and then calculaitng hte mean of each populaiton in each pre-specified euclidean dimension,
#' but first identifies the median in the gate dimension for each population, and then calculates the mean in each euclidean dimension outside of this blurred middle range. To increase comparability to normal gating, an increased weight is normally
#' put on the gating variable in the final euclidean distance calculation, so that the overlap between the populations is somewhat lowered.
#' @param inDatDonList This is the list of donor dataframes that should be gated. Here, all variables should be included, and the specific gating variables are specified below.
#' @param negDonVec A numeric vector. If some donors should be left ungated, for example because they
#' lack expression of the marker of interest, they should be added here,
#' and all cell will be regarded asnegative for the marker.
#' @param gateMarker Just like it soounds, this specifies which of the columns in each of the donor dataframes that should be used to create the gates. Can be given as a number or a name.
#' @param euclidMarkers This includes all markers that should be used for euclidean neighbor calculations. The gateMarker should also be included here. Can be given as numbers or names.
#' @param lowPeak Logical, defining if the low or the high peak in the data should be the basis for the gates. Currently, only uni- or bimodal distributions are correctly interpreted.
#' @param abovePeak Logical, defining at which side of the peak that the gate should be created.
#' @param adjust A command changing the bandwidth of the density function, both for calculation of the peaks and gate positions and for visualization.
#' @param nDens The number of bins that are used to calculate the density. "default" means the lowest of 512 and the number of events divided by 5. A specific number can also be specified.
#' @param volRatio Defining how small the smaller of the two peaks can be to be considered a true peak. It is a fraction of the volume of the larger peak. Default is 0.05, i.e. if the volume of the second peak is 5 percent or larger than the volume of the first peak, it is considered a peak.
#' @param gateWeight This command defines how similar the selection step is to normal gating. A high value will render the overlap between the populations very small, whereas a low value will have the opposite effect. 2, i.e. a weight of 2x that of any of the other variables, is standard.
#' @param createDir If a new directory should be created.
#'
#' @return A list of vectors with the gate information for each event in each donor.
#'
#' @examples
#' # Here is a simple example
#' library(DepecheR)
#' data(testData)
#' testDataList <- split(testData, f = testData$ids)
#' # This will be the diretory where the folder will be created:
#' getwd()
#'
#' testGates <- turnDerGate(testDataList, 3, c(2:15), adjust = 2)
#' @export turnDerGate
turnDerGate <- function(inDatDonList, negDonVec = 0, gateMarker, euclidMarkers, lowPeak = TRUE, abovePeak = TRUE, adjust = 2, nDens = "default", volRatio = 0.05, gateWeight = 2, createDir = TRUE) {
    if (any(is(gateMarker) == "numeric")) {
        gateMarker <- colnames(inDatDonList[[1]])[gateMarker]
    }
    if (createDir == TRUE) {
        oldDir <- getwd()
        newDir <- paste0("./Individual_", gateMarker, "_gates")
        dir.create(newDir)
        setwd(newDir)
    }

    resultList <- list()
    for (i in 1:length(inDatDonList)) {
        if(i %in% negDonVec){
            resultList[[i]] <- rep(1, nrow(inDatDonList[[i]]))
        } else {
            markerData <- inDatDonList[[i]][, gateMarker]
            if (nDens == "default") {
                n <- min(64, length(markerData) / 5)
            } else {
                n <- nDens
            }
            Da = density(markerData, adjust = adjust, n = n)
            DeltaY = diff(Da$y)
            Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
            # Identify the two highest turnpoints
            topTwo <- Turns[order(Da$y[Turns], decreasing = TRUE)][1:2]

            # And here, it is made clear if it is meaningful to talk about two peaks, or if it is indeed only one.
            # First, it is determined if there are two peaks whatsoever, and if so, by how much these are separated. Currently, they need to be separated by at least 5% of the total range from the first to th 99th percentile to count as two peaks.
            robRange <- quantile(markerData, c(0.1, 0.9))
            tenth <- abs(unname((robRange[2] - robRange[1]) * 0.1))
            if (is.na(topTwo[2]) == FALSE && abs(Da$x[topTwo][1] - Da$x[topTwo][2]) > tenth) {
                # Now, if two peaks are present, identify the lowest turnpoint between the topTwo
                betweenTurns <- Turns[Turns < max(topTwo) & Turns > min(topTwo)]
                minTurn <- betweenTurns[which.min(Da$y[betweenTurns])]
                if (topTwo[1] < topTwo[2]) {
                    startHigh <- Da$x[1]
                    stopHigh <- Da$x[minTurn]
                    startLow <- Da$x[minTurn]
                    stopLow <- Da$x[length(Da$x)]
                    topTurns <- c(topTwo[1], topTwo[2])
                } else {
                    startLow <- Da$x[1]
                    stopLow <- Da$x[minTurn]
                    startHigh <- Da$x[minTurn]
                    stopHigh <- Da$x[length(Da$x)]
                    topTurns <- c(topTwo[2], topTwo[1])
                }

                volHigh <- length(markerData[which(markerData > startHigh & markerData < stopHigh)])
                volLow <- length(markerData[which(markerData > startLow & markerData < stopLow)])

                # And here, it is determined if the second peak counts, which it only does if it has a volume that is at least a set percentage of the largest one. Default being 1%.
                if (volHigh > 0 && volLow > 0 && volLow / volHigh > volRatio) {
                    twoPeaks <- TRUE
                    print("Two peaks are detected.")
                    if (lowPeak == TRUE && abovePeak == TRUE) {
                        print("As the gate should be placed between the peaks, the turn point closest to the lower peak was chosen")
                        gateVal <- Da$x[Turns[which(Turns == topTurns[1]) + 1]]
                        graphName <- paste0("Donor_", names(inDatDonList)[i], "_", gateMarker, "_turnGate.pdf")

                    } else if (lowPeak == FALSE && abovePeak == FALSE) {
                        print("As the gate should be placed between the peaks, the turn point closest to the higher peak was chosen")
                        gateVal <- Da$x[Turns[which(Turns == topTurns[2]) - 1]]
                        graphName <- paste0("Donor_", names(inDatDonList)[i], "_", gateMarker, "_turnGate.pdf")

                    } else {
                        print("It was not suitable to use the turnpoint alternative, as the gate should not be placed between the peaks")
                        peakCenter <- topTwo[1]
                        twoPeaks <- FALSE
                    }

                } else {
                    peakCenter <- topTwo[1]
                    twoPeaks <- FALSE
                }
            } else {
                peakCenter <- topTwo[1]
                twoPeaks <- FALSE
            }
            if (twoPeaks == FALSE) {
                # Only one peak was detected with the current threshold on peak size, and thus, a gate based on curve derivatives is created instead

                graphName <- paste0("Donor_", names(inDatDonList)[i], "_", gateMarker, "_derGate.pdf")

                if (abovePeak == TRUE) {
                    peakRange <- c(peakCenter:length(Da$y))
                    DeltaY2 = diff(diff(Da$y[peakRange]))
                    subDeltaY2 <- unlist(sapply(1:length(DeltaY2) - 1, function(x) DeltaY2[x] * DeltaY2[x + 1]))
                    # Here, the position after the first passing of zero for the second derivative is chosen, i e the first deflection point in the firt derivative after the peak.
                    delta2Val <- (which(subDeltaY2 < 0) + (peakCenter + 1))[1]
                    zoomRange <- c(delta2Val:length(Da$y))
                    DeltaY4 = diff(diff(diff(diff(Da$y[zoomRange]))))
                    subDeltaY4 <- unlist(sapply(1:length(DeltaY4) - 1, function(x) DeltaY4[x] * DeltaY4[x + 1]))
                    deltaPosTurns <- which(subDeltaY4 < 0) + 1
                    Turns4 <- deltaPosTurns + delta2Val
                    gateInterval <- data.frame("x" = c(Da$x[Turns4[1]], Da$x[Turns4[1] - 1]), "y" = c(DeltaY4[deltaPosTurns[1]], DeltaY4[deltaPosTurns[1] - 1]))
                    gateVal <- lm(x ~ y, gateInterval)$coefficients[1]

                } else {
                    peakRange <- c(1:peakCenter)
                    DeltaY2 = diff(diff(Da$y[peakRange]))
                    subDeltaY2 <- unlist(sapply(1:length(DeltaY2) - 1, function(x) DeltaY2[x] * DeltaY2[x + 1]))
                    # Here, the position after the first passing of zero for the second derivative is chosen, i e the first deflection point in the firt derivative after the peak.
                    allDelta2Vals <- which(subDeltaY2 < 0)
                    delta2Val <- allDelta2Vals[length(allDelta2Vals)]
                    zoomRange <- c(1:delta2Val)
                    DeltaY4 = diff(diff(diff(diff(Da$y[zoomRange]))))
                    subDeltaY4 <- unlist(sapply(1:length(DeltaY4) - 1, function(x) DeltaY4[x] * DeltaY4[x + 1]))
                    Turns4 <- which(subDeltaY4 < 0) + 1
                    TurnVal <- Turns4[length(Turns4)]
                    gateInterval <- data.frame("x" = c(Da$x[TurnVal], Da$x[TurnVal - 1]), "y" = c(DeltaY4[TurnVal], DeltaY4[TurnVal - 1]))
                    gateVal <- lm(x ~ y, gateInterval)$coefficients[1]

                }
            }

            resultList[[i]] <-
                turnDerGateCoFunction(euclidFocus =
                                          inDatDonList[[i]][, euclidMarkers],
                                      gateMarker = gateMarker,
                                      gateVal = gateVal,
                                      gateWeight = gateWeight,
                                      graphName = graphName,
                                      adjust = adjust,
                                      n = n)
            }
        print(paste0("Done with individual ", names(inDatDonList)[i]))

        }

    if (createDir == TRUE) {
        setwd(oldDir)
    }
    return(resultList)
}
