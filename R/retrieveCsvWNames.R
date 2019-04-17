retrieveCsvWNames <- function(filename, idInfo, nrows = -1) {
    ret <- read.csv(filename, nrows = nrows)
    if (missing(idInfo) == TRUE) {
        ret$names <- filename
    } else
        for (i in 1:length(idInfo)) {
            ret$names <- gsub(pattern = idInfo[[i]], "\\1", filename)
            colnames(ret)[which(colnames(ret) == "names")] <- names(idInfo)[i]
        }
    ret
}
