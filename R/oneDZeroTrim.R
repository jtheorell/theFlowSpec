#' importFrom flowCore exprs
#' importFrom BiocGenerics nrow ncol colnames
#' @export oneDZeroTrim
oneDZeroTrim <- function(focusFrame, trimFrac, nRowOut, returnFlowFrame = TRUE,
                         returnPositions = FALSE){
    dataLength <- BiocGenerics::nrow(focusFrame)
    if(dataLength < nRowOut){
        nRowOut <- dataLength
    }

    dataExtract <- as.data.frame(matrix(data = NA, nrow = nRowOut,
                                        ncol = BiocGenerics::ncol(focusFrame)))

    colnames(dataExtract) <- BiocGenerics::colnames(focusFrame)

    for(i in seq_len(BiocGenerics::ncol(focusFrame))){
        focCol <- exprs(focusFrame)[,i]
        focPositions <- seq_along(focCol)
        minimum <- min(focCol)
        maximum <- max(focCol)
        fraction1 <- c(minimum, abs(maximum-minimum)/100)
        minPos <- which(focCol >= fraction1[1] & focCol < fraction1[2])

        #And here, the reduction of the zero fraction comes.
        if(length(minPos) > dataLength*trimFrac){
            dataPositions <- c(focPositions[-minPos],
                               sample(minPos, dataLength*trimFrac))
        } else {
          dataPositions <- focPositions
        }

        if(length(dataPositions) > nRowOut){
          dataPositions <- sample(dataPositions, nRowOut)
            dataExtract[,i] <- focCol[dataPositions]
        } else {

            dataExtract[seq(1, length(dataPositions)),i] <-
              focCol[dataPositions]
        }
    }

    if(returnFlowFrame){
        return(flowFrame(as.matrix(dataExtract)))
    } else if(returnPositions){
      return(dataPositions)
    }
    return(dataExtract)

}
