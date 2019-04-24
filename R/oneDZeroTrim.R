#' importFrom flowCore exprs
#' importFrom BiocGenerics nrow ncol colnames
#' @export oneDZeroTrim
oneDZeroTrim <- function(focusFrame, trimFrac, nRowOut, returnFlowFrame = TRUE){
    dataLength <- BiocGenerics::nrow(focusFrame)
    if(dataLength < nRowOut){
        nRowOut <- dataLength
    }
    
    dataExtract <- as.data.frame(matrix(data = NA, nrow = nRowOut, 
                                        ncol = BiocGenerics::ncol(focusFrame)))
    
    colnames(dataExtract) <- BiocGenerics::colnames(focusFrame)
    
    for(i in seq_len(BiocGenerics::ncol(focusFrame))){
        focCol <- exprs(focusFrame)[,i]
        minimum <- min(focCol)
        maximum <- max(focCol)
        fraction1 <- c(minimum, abs(maximum-minimum)/100)
        minPos <- which(focCol >= fraction1[1] & focCol < fraction1[2])
        
        #And here, the reduction of the zero fraction comes.
        if(length(minPos) > dataLength*trimFrac){
            dataCol <- c(focCol[-minPos], 
                         focCol[minPos[sample(seq_along(minPos), 
                                              dataLength*trimFrac)]])
        } else {
            dataCol <- focCol
        }
        
        if(length(dataCol) >= nRowOut){
            dataExtract[,i] <- dataCol[sample(seq_along(dataCol), 
                                              nRowOut)]
        } else {
            dataExtract[seq(1, length(dataCol)),i] <- dataCol
        }
    }
    
    if(returnFlowFrame){
        return(flowFrame(as.matrix(dataExtract)))
    }
    return(dataExtract)

}