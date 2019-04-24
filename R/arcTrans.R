#' Efficient inverse hyperbolic cosine transformation
#' 
#' This is a simple wrapper function for the base asinh function, that is useful
#' for flowFrames and flowSets. It also allows for reversing the transformation
#' with the argument reverse
#' @param flowObj The fcs object to be transformed. Both flowFrames and flowSets
#' are accepted. 
#' @param transNames The variables that should be normalized. Default is all
#' that do not fall into the excl category below.
#' @param exclColNmStr A vector of strings that are not accepted in the
#' transNames.
#' @param transCoFacs This vector of values define the values for the
#' transformation during the normalization. In the "default" case, the function
#' defines the object as a CyTOF object if >5 percent of the values are 0, and 
#' applies the transformation value 8. Otherwise, the value 256 is applied.
#' @param unTrans If the reverse action should be taken, i.e. if an already 
#' transformed dataset should be un-transformed. NB! It is of great importance
#' that the same transformation factors are used!
#' 
#' @return A flow object containing the transformed data, and with all metadata
#' left untouched. 
#' @importFrom BiocGenerics colnames
#' @export arcTrans
arcTrans <- function(flowObj, transNames = "default",
                     exclColNmStr = c("ime"), transCoFacs = "default", 
                     unTrans = FALSE){
    if(class(flowObj) == "flowSet"){
        focusFrame <- flowObj[[1]]
    } else if(class(flowObj) == "flowFrame"){
        focusFrame <- flowObj
    } else {
        stop("The flowObj needs to be either a flowSet or a flowFrame")
    }
    
    # Here, the variables that the normalization should be performed on are 
    #selected.
    if(length(transNames) ==1 && transNames == "default"){
        transNames <- BiocGenerics::colnames(focusFrame)
    }
    
    # Now, the time variable, and possibly other variables are excluded
    excludedTransNames <- grep(paste(exclColNmStr, collapse = "|"),
                               transNames)
    if(length(excludedTransNames) > 0){
        transNames <- transNames[-excludedTransNames] 
    }
    

    if(length(transCoFacs) ==1 && transCoFacs == "default"){
        transCoFacs <- massOrFlowTrans(focusFrame = focusFrame, 
                                       transNames = transNames)
    } 

    if(class(flowObj) == "flowSet"){
        resultObj <- fsApply(flowObj, theFlowSpec:::arcTransCoFunction, 
                             transNames = transNames, transCoFacs = transCoFacs,
                             unTrans = unTrans)
    } else {
        resultObj <- arcTransCoFunction(flowObj, transNames = transNames, 
                                        transCoFacs = transCoFacs)
    }
    return(resultObj)
        
}

arcTransCoFunction <- function(focusFrame, transNames, transCoFacs, unTrans){
    focusFrameResult <- focusFrame
    if(unTrans){
        for(i in transNames){
            exprs(focusFrameResult)[,i] <- 
                sinh(exprs(focusFrame)[,i])*transCoFacs[i]
        } 
    } else {
        for(i in transNames){
            exprs(focusFrameResult)[,i] <- 
                asinh(exprs(focusFrame)[,i]/transCoFacs[i])
        } 
    }
    
    message("Transformation complete")
    
    return(focusFrameResult)
}