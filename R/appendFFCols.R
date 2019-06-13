#This function is a slightly adjusted copy of the fr_append_cols 
#function in flowCore. It mainly allows for the presens of the additional pData 
#variable oldName
appendFFCols <- function(focusFrame, newCols){
    pd <- pData(parameters(focusFrame))
    cn <- colnames(newCols)
    new_pid <- max(as.integer(gsub("\\$P", "", rownames(pd)))) + 
        1
    new_pid <- seq(new_pid, length.out = ncol(newCols))
    new_pid <- paste0("$P", new_pid)
    new_pd <- do.call(rbind, lapply(cn, function(i) {
        vec <- newCols[, i]
        rg <- range(vec)
        data.frame(name = i, desc = NA, range = diff(rg) + 1, 
                   minRange = rg[1], maxRange = rg[2], oldName = NA)
    }))
    rownames(new_pd) <- new_pid
    pd <- rbind(pd, new_pd)
    focusFrame@exprs <- cbind(exprs(focusFrame), newCols)
    pData(parameters(focusFrame)) <- pd
    focusFrame
}