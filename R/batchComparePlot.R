#' @importFrom flowCore exprs
#' @importFrom BiocGenerics nrow colnames
#' @importFrom DepecheR dColorVector
#' @importFrom beanplot beanplot
#' @export batchComparePlot
#' @export batchComparePlotCoFunc
batchComparePlot <- function(ctrlSet, normSet, pathName = "Batch_comparisons",
                             trimZero = TRUE){

    #First some sanity checks
    if(identical(BiocGenerics::colnames(ctrlSet),
                 BiocGenerics::colnames(normSet))==FALSE){
        stop("The controlset and normalisation set need to have identical
             column names")
    }

    #First, a relevantly sized subset of each file is included in two long
    #frames, one for the control set and one for the normalized set.
    ctrlLong <- batchComparePlotCoFunc(ctrlSet, trimZero = trimZero)
    normLong <- batchComparePlotCoFunc(normSet, trimZero = trimZero)

    message("Done with the pre-processing")

    allLong <- rbind(ctrlLong, normLong)

  areaColors <- as.list(dColorVector(
    c(rep(1, length.out=length(sampleNames(ctrlSet))),
      rep(2, length.out=length(sampleNames(normSet))))))

  plotNames <- BiocGenerics::colnames(ctrlSet)

  #Now, the directory is created, if not already present.
  dir.create(pathName, showWarnings = FALSE)

  #And here, the plots are created.
  for(i in plotNames){
    pdf(file.path(pathName, paste0(i, ".pdf")))
    beanplot(split(allLong[,i], allLong$id),what=c(0,1,1,0), col=areaColors)
    dev.off()
  }
}

batchComparePlotCoFunc <- function(focusSet, trimZero){
    focusList <- lapply(seq_along(focusSet), function(x) {
        focusFrame <- focusSet[[x]]

        #Here, to reduce the negative influence on the plots of the often
        #extreme sparsity of CyTOF data, the data in the lowest
        #percent of the range is is capped to 10% of the total.
        if(trimZero){
            dataExtract <- oneDZeroTrim(focusFrame, trimFrac = 0.1,
                                        nRowOut = 10000,
                                        returnFlowFrame = FALSE)
            } else {
                if(BiocGenerics::nrow(focusFrame)>10000){
                    dataExtract <-
                        as.data.frame(exprs(focusFrame)[sample(seq_len(
                            BiocGenerics::nrow(focusFrame)), 10000),])
                    } else {
                        dataExtract <- as.data.frame(exprs(focusFrame))
                    }
            }

        #Now, an identifier is added to the matrix in question
        dataExtract$id <- rep(sampleNames(focusSet)[x],
                              times = nrow(dataExtract))
        return(dataExtract)
    })
    completeFrame <- do.call("rbind", focusList)
}
