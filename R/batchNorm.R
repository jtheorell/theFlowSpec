#' Batch normalisation
#'
#'
#' This function is intended to be used, when files from the same individual has
#' been acquired on multiple time points, and the files from the different
#' time points should be compared. The function works most optimally if all the
#' data came from a sample that was processed and frozen the same day, i.e. a 
#' true technical control, but it also works reasonably well with data from the
#' same donor from different time points, as it is robust to outlier populations
#' changing between the dates. It is not recommended to apply this function to
#' files with less than 5000 events.
#' @param normSet A flowSet that should be normalized. If the intCtrl should
#' be normalized and used downstream, it needs to be included here.
#' @param intCtrl A control flowFrame that was acquired the same day as the
#' normSet.
#' @param extCtrl A control flowFrame that the normalization should be performed
#' against.
#' @param normNames The variables that should be normalized. Default is all that
#' do not fall into one of the two excl categories below.
#' @param exclColNmStr A vector of strings that are not accepted in any
#' normNames.
#' @param exclCategorical Logical: should normNames with less than 10
#' unique values be excluded? Defaults to TRUE.
#' @param volThresh Defining how small the smaller of the two peaks can be to be
#' considered a true peak. It is a fraction of the volume of the larger peak.
#' Default is 0.05, i.e. if the volume of the second peak is 5 percent or larger
#' than the volume of the first peak, it is considered a peak.
#' @param transCoFacs This vector of named values define the values for the
#' transformation during the normalization. This is only applied internally, so
#' transformation needs to be performed afterwards, preferrably with individual
#' values for each channel. In the "default" case, the function defines the file
#' as a CyTOF file, and applies the transformation value 8, if >5 percent of the
#' values are 0. Otherwise, the value 256 is applied. NB. The entries need to be
#' named in the same way as the normNames to secure that the right factor is 
#' added to each variable. 
#' @param zeroTrim In the case of CyTOF data, the events at zero can often
#' be so dominant, that all other density variation is dwarfed, and thus
#' invisible. With this command being TRUE, the events that are confined to the
#' first percentile of the data range for each variable are trimmed
#' to max 10 percent of the data.
#'
#' @return The normSet normalized. If the intCtrl was not among the normSet
#' files, it will be added as the final object to the list.
#'
#' @importFrom BiocGenerics nrow
#' @export batchNorm
batchNorm <- function(normSet, intCtrl, extCtrl,
                      normNames = colnames(intCtrl),
                      exclColNmStr = c("ime"), volThresh = 0.02,
                      transCoFacs = "default", zeroTrim = TRUE) {

    # First, the nrow of the norm- and controlFiles are evaluated, so that they
    # both contain more than 5000 events.
    if (BiocGenerics::nrow(intCtrl) < 5000 ||
        BiocGenerics::nrow(extCtrl) < 5000) {
        stop("Either the norm- or the extCtrl contains less than 5000
              events, which is the (completely arbitrary) lower border for
              useful normalization")
    }

    # Here, zero trimming and downsampling is performed
    if(zeroTrim){
        intCtrl <- oneDZeroTrim(focusFrame = intCtrl, trimFrac = 0.1, 
                                nRowOut = 100000)
        extCtrl <- oneDZeroTrim(focusFrame = extCtrl, trimFrac = 0.1, 
                                nRowOut = 100000)
    } else{
        if(BiocGenerics::nrow(intCtrl) > 100000){
            intCtrl <- intCtrl[sample(seq_len(nrow(intCtrl)), 100000), ]
        } 
        if(BiocGenerics::nrow(extCtrl) > 100000) {
            extCtrl <-
                extCtrl[sample(seq_len(nrow(extCtrl)), 100000), ]
        }
    } 

    # Here, the variables that the normalization should be performed on is chosen
    # First, the time variable, and possibly other variables are excluded
    excludedNormNames <- grep(paste(exclColNmStr, collapse = "|"),
                               normNames)
    
    if(length(excludedNormNames) > 0){
        normNames <- normNames[-excludedNormNames] 
    }

    #Here, the transformation value is defined, depending on if the expression
    #below results in a value above 0.05 or not.
    if(length(transCoFacs) ==1 && transCoFacs == "default"){
        transCoFacs <- massOrFlowTrans(focusFrame = extCtrl, 
                                       transNames = normNames)
        }

    # Here, the essential confunction comes into play
    newNormVals <- lapply(seq_along(normNames), function(i)
        batchNormCoFunc1(intVar = exprs(intCtrl[, normNames[i]])[,1],
                         extVar = exprs(extCtrl[, normNames[i]])[,1],
                         transCoFac = transCoFacs[normNames[i]], 
                         volThresh = volThresh))
    
    names(newNormVals) <- normNames
    
    # Here depending on if the data is from CyTOF or flow cytometry the analysis
    #takes different paths. In the case of CyTOF, the zero values are always the
    #lower peak and they should not be moved.
    exprsInt <- exprs(intCtrl)
    dimExprsInt <- dim(exprsInt)
    if(sum(exprsInt == 0, na.rm = TRUE) >0.05*dimExprsInt[1]*dimExprsInt[2] && 
       sum(exprsInt < 0, na.rm = TRUE) == 0){
        message("The data is supposed to be CyTOF-data, and normalized 
                accordingly")
        newNormVals <- lapply(newNormVals, function(x) {
            if(length(x) == 2){
                x[c(1,2)] <- c(0,0)
            } else {
                x[c(1,3)] <- c(0,0)
            }
            return(x)
        })
    } else {
        message("The data is supposed to be flow cytometry data, and normalized 
                accordingly")
    }
    
    # detected. If two, the data is first normalized to its internal peaks, and
    # then re-organized with the information from the other peaks. If, on the
    # other hand there is only one peak to normalize to, only centering is
    # performed with this peak in mind.

    normSetNorm <- fsApply(normSet, batchNormCoFunc2,
                                newNormVals = newNormVals,
                                normNames = normNames)
    return(normSetNorm)
}
