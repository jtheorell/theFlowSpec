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
#' @param normSet A flowSet that should be normalized. If the normFrame should 
#' be normalized, it needs to be included here.
#' @param normFrame The flowFrame that contains the events which will be used 
#' for the normalization. 
#' @param controlFrame The flowFrame containing the control values. 
#' @param normNames The variables that should be normalized. Default is all that
#' do not fall into one of the two excl categories below.
#' @param exclColNmStr A vector of strings that are not accepted in the
#' colnames of any normNames.
#' @param exclCategorical Logical: should normNames with less than 10
#' unique values be excluded? Defaults to TRUE.
#' @param volRatio Defining how small the smaller of the two peaks can be to be
#' considered a true peak. It is a fraction of the volume of the larger peak.
#' Default is 0.05, i.e. if the volume of the second peak is 5 percent or larger
#' than the volume of the first peak, it is considered a peak.
#' @param transCoFac This value defines the standard value for the 
#' transformation during the normalization. This is only applied internally, so
#' transformation needs to be performed afterwards, preferrably with individual
#' values for each channel. In the "default" case, the function defines the file
#' as a CyTOF file, and applies the transformation value 8, if >5% of the values
#' are 0. Otherwise, the value 256 is applied. 
#'
#' @return The normSet normalized. If the normFrame was not among the normSet
#' files, it will be added as the final object to the list.
#'
#' @importFrom BiocGenerics nrow
#' @importFrom flowVS transFlowVS
#' @export batchNorm
batchNorm <- function(normSet, normFrame, controlFrame,
                      normNames = colnames(normFrame),
                      exclColNmStr = c("ime"), volRatio = 0.1, 
                      transCoFac = "default") {

    # First, the nrow of the norm- and cotrolFiles are evaluated, so that they
    # both contain more than 5000 events.
    if (BiocGenerics::nrow(normFrame) < 5000 || 
        BiocGenerics::nrow(controlFrame) < 5000) {
        stop("Either the norm- or the controlFrame contain less than 5000
              events, which is the (arbitrary) lower border for any accurate 
              estimation of nearest neighbors and therefore useful 
              normalization")
    }
    
    # Now, if the normFrame and the controlFrame are very large, they are 
    #reduced somwehat
    if (BiocGenerics::nrow(normFrame) > 100000) {
        normFrame <- normFrame[sample(seq_len(nrow(normFrame)), 100000), ]
    }
    if (BiocGenerics::nrow(controlFrame) > 100000) {
        controlFrame <-
            controlFrame[sample(seq_len(nrow(controlFrame)), 100000), ]
    }
    

    # Here, the variables that the normalization should be performed on is chosen
    # First, the time variable, and possibly other variables are excluded
    normNames <- normNames[-grep(paste(exclColNmStr, collapse = "|"), 
                                 normNames)]

    #Here, the transformation value is defined, depending on if the expression
    #below results in a value above 0.05 or not. 
    if(transCoFac == "default" && sum(exprs(controlFrame) == 0)/
       (dim(controlFrame)[1]*dim(controlFrame)[2]) >0.05){
        transCoFac <- 8
    } else {
        transCoFac <- 256
    }
    
    #Now, the normalization files are transformed, using the transformation 
    #value above, and the transFlowVS function
    
    normFrameTrans <- transFlowVS(flowSet(normFrame), normNames, 
                cofactors=rep(transCoFac, length.out=length(normNames)))[[1]]
    controlFrameTrans <- transFlowVS(flowSet(controlFrame), 
                                     normNames, 
                                     cofactors=
                                         rep(transCoFac, 
                                             length.out=length(normNames)))[[1]]
    
    # Here, the essencial confunction comes into play
    newNormVals <- lapply(seq_along(normNames), function(i)
        batchNormCoFunc1(normVar = exprs(normFrame[, normNames[i]])[,1],
                         controlVar = exprs(controlFrame[, normNames[i]])[,1],
                         volRatio = volRatio))
    
    # Now, the analysis takes two paths, depending on if one or two peaks were
    # detected. If two, the data is first normalized to its internal peaks, and
    # then re-organized with the information from the other peaks. If, on the
    # other hand there is only one peak to normalize to, only centering is
    # performed with this peak in mind.

    normSetNormLists <- fsAapply(normSet, function(x)
        lapply(seq_along(newNormVals), function(y)
            batchNormCoFunc2(newNormValsFocus = newNormVals[[y]],
                focusColnormFrame = x[, normNames[y]])))

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
