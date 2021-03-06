% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/turnDerGate.R
\name{turnDerGate}
\alias{turnDerGate}
\title{Function used to generate automated gates on large sets of samples}
\usage{
turnDerGate(inDatDonList, negDonVec = 0, gateMarker, euclidMarkers,
  lowPeak = TRUE, abovePeak = TRUE, adjust = 2, nDens = "default",
  volRatio = 0.05, gateWeight = 2, createDir = TRUE)
}
\arguments{
\item{inDatDonList}{This is the list of donor dataframes that should be gated. Here, all variables should be included, and the specific gating variables are specified below.}

\item{negDonVec}{A numeric vector. If some donors should be left ungated, for example because they
lack expression of the marker of interest, they should be added here,
and all cell will be regarded asnegative for the marker.}

\item{gateMarker}{Just like it soounds, this specifies which of the columns in each of the donor dataframes that should be used to create the gates. Can be given as a number or a name.}

\item{euclidMarkers}{This includes all markers that should be used for euclidean neighbor calculations. The gateMarker should also be included here. Can be given as numbers or names.}

\item{lowPeak}{Logical, defining if the low or the high peak in the data should be the basis for the gates. Currently, only uni- or bimodal distributions are correctly interpreted.}

\item{abovePeak}{Logical, defining at which side of the peak that the gate should be created.}

\item{adjust}{A command changing the bandwidth of the density function, both for calculation of the peaks and gate positions and for visualization.}

\item{nDens}{The number of bins that are used to calculate the density. "default" means the lowest of 512 and the number of events divided by 5. A specific number can also be specified.}

\item{volRatio}{Defining how small the smaller of the two peaks can be to be considered a true peak. It is a fraction of the volume of the larger peak. Default is 0.05, i.e. if the volume of the second peak is 5 percent or larger than the volume of the first peak, it is considered a peak.}

\item{gateWeight}{This command defines how similar the selection step is to normal gating. A high value will render the overlap between the populations very small, whereas a low value will have the opposite effect. 2, i.e. a weight of 2x that of any of the other variables, is standard.}

\item{createDir}{If a new directory should be created.}
}
\value{
A list of vectors with the gate information for each event in each donor.
}
\description{
In addition to performing optimized gating based on either the most suitable turning point between two populaitons in the case of a bimodal distritbution,
or finding a specific turning point in the third derivative in the case of a single peak distribution, the function goes on to identify which cells that are most
similar to the cluster center of each of the two populations created by the gate. It does this not just by dividing along the gate and then calculaitng hte mean of each populaiton in each pre-specified euclidean dimension,
but first identifies the median in the gate dimension for each population, and then calculates the mean in each euclidean dimension outside of this blurred middle range. To increase comparability to normal gating, an increased weight is normally
put on the gating variable in the final euclidean distance calculation, so that the overlap between the populations is somewhat lowered.
}
\examples{
# Here is a simple example
library(DepecheR)
data(testData)
testDataList <- split(testData, f = testData$ids)
# This will be the diretory where the folder will be created:
getwd()

testGates <- turnDerGate(testDataList, 3, c(2:15), adjust = 2)
}
