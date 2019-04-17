#' Neighbor batch normalisation
#'
#'
#' This function is intended to be used, when files from the same individual has
#' been acquired on multiple time points, and the files from the different
#' time points should be compared. The function works most optimally if all the
#' data came from a sample that was processed and frozen the same day, i.e. a true
#' technical control, but it also works reasonably well with data from the same
#' donor from different time points, as it is robust to outlier populations
#' changing between the dates. It is not recommended to apply this function to
#' files with less than 5000 events, and more than 100 000 events are preferrable.
#' @param normSet A list of files that should be batch normalized.
#' @param normFile The file that contains the events which will be used for the
#' normalization.
#' @param controlFile The data cloud in which the nearest neighbors for the
#' normFile events should be identified.
#' @param normNames The variables that should be normalized. Default is all that
#' do not fall into one of the two excl categories below.
#' @param exclColNmStr A vector of strings that are not accepted in the
#' colnames of any normNames.
#' @param exclCategorical Logical: should normNames with less than 10
#' unique values be excluded? Defaults to TRUE.
#' @param kNeighK The number of nearest neighbors that each of the normalization
#' band events are compared to. The higher the number, the more robust, and the
#' less sensitive to outliers, will the  algorithm be.
#' @param volRatio Defining how small the smaller of the two peaks can be to be
#' considered a true peak. It is a fraction of the volume of the larger peak.
#' Default is 0.05, i.e. if the volume of the second peak is 5 percent or larger
#' than the volume of the first peak, it is considered a peak.
#'
#' @return The normSet normalized. If the normFile was not among the normSet
#' files, it will be added as the final object to the list.
#'
#' @export batchNorm
batchNorm <- function(normSet, normFile, controlFile,
                      normNames = colnames(normFile),
                      exclColNmStr = c("Id", "id", "onor", "ample"),
                      exclCategorical = TRUE,
                      kNeighK = 10, volRatio = 0.1) {

    # First, the nrow of the norm- and cotrolFiles are evaluated, so that they
    # both contain more than 5000 events.
    if (nrow(normFile) < 5000 || nrow(controlFile) < 5000) {
        stop("Either the norm- or the controlFile contain less than 5000
              events, which is the (arbitrary) lower border for any accurate 
              estimation of nearest neighbors and therefore useful 
              normalization")
    }
    # Here, it is evaluated  if the normFile should be added to the normSet, if
    # it is not included in it.
    if (any(as.logical(lapply(seq_along(normSet),
        function(x) identical(normSet[[x]], normFile)))) == FALSE) {
        normSet <- c(normSet, normFile)
    }

    # Now, if the normFile and the controlFile are very large, they are reduced
    # somwehat
    if (nrow(normFile) > 100000) {
        normFile <- normFile[sample(seq_len(nrow(normFile)), 100000), ]
    }
    if (nrow(controlFile) > 1000000) {
        controlFile <- controlFile[sample(seq_len(nrow(controlFile)), 100000), ]
    }

    # Here, the variables that the normalization should be performed on is chosen
    # First, all clearly id-related variables are excluded
    normNames <- normNames[-grep(paste(exclColNmStr, collapse = "|"), normNames)]

    if (exclCategorical) {
        includableVars <-
            colnames(normFile)[which(vapply(normFile,
                function(x) length(unique(x)),
                1) > 10)]
        normNames <- normNames[which(normNames %in% includableVars)]
    }

    # Now, the normFile and the controlFile are reduced to the normNames
    newNormVals <- lapply(seq_along(normNames), function(i)
        batchNormCoFunc1(normVar = normFile[, normNames[i]],
            controlVar = controlFile[, normNames[i]],
            volRatio = volRatio))

    # Now, the analysis takes two paths, depending on if one or two peaks were
    # detected. If two, the data is first normalized to its internal peaks, and
    # then re-organized with the information from the other peaks. If, on the
    # other hand there is only one peak to normalize to, only centering is
    # performed with this peak in mind.

    normSetNormLists <- lapply(normSet, function(x)
        lapply(seq_along(newNormVals), function(y)
            batchNormCoFunc2(newNormValsFocus = newNormVals[[y]],
                focusColNormFile = x[, normNames[y]])))

    normSetNorm <- lapply(normSetNormLists, function(x)
        as.data.frame(do.call("cbind", x)))

    # And now, the colnames and the initially excluded columns are brought back
    normSetComplete <- normSetNorm
    for (i in seq_along(normSetNorm)) {
        colnames(normSetNorm[[i]]) <- normNames
        normSetComplete[[i]] <-
            cbind(normSetNorm[[i]],
                normSet[[i]][, -which(colnames(normSet[[i]]) %in% normNames)])
    }

    return(normSetComplete)
}
