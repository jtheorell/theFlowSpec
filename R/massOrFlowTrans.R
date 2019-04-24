#' @export massOrFlowTrans
massOrFlowTrans <- function(focusFrame, transNames){
    if(unname(sum(exprs(focusFrame) == 0)/(dim(focusFrame)[1]*
                                           dim(focusFrame)[2])) >0.05){
        transCoFacs <- rep(8, times = length(transNames))
    } else {
        transCoFacs <- rep(256, times = length(transNames))
    }
    names(transCoFacs) <- transNames
    
    return(transCoFacs)
}