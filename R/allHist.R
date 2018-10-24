#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes facet_wrap geom_density theme element_blank element_line ggsave
#' @export
allHist <- function(dataSet, name="dataSet"){
  dataSet %>%
    gather() %>%
    ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_density() + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"))
  ggsave(paste0(name, ".pdf"))
}

