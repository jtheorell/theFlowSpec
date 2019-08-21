massOrFlowTrans <- function(focusFrame, transNames){
    if(unname(sum(exprs(focusFrame) == 0)/(dim(focusFrame)[1]*
                                           dim(focusFrame)[2])) >0.05){
        transCoFacs <- rep(5, times = length(transNames))
    } else {
        transCoFacs <- rep(200, times = length(transNames))
    }
    names(transCoFacs) <- transNames

    return(transCoFacs)
}
