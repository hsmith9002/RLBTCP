#' Summarize transcript content into a tabular format
#' 
#' This takes as input the outlist from content.for.plot()
#' and creates a table that summarized the number of transcripts
#' in each category along with totals and percentages.
#' @param dataset outlist from content.for.plot()
#' @export
#' @examples 
#' summarize.content(content.summary)

summarize.content <- function(dataset) {
  library(plyr)
  alo <- dataset$ALO.ForPlotSS
  all <- dataset$ALL.ForPlotSS
  contentDF <- do.call("rbind", llply(list(alo, all), .fun = function(x) {alo.pib <- c(x$value[1], x$value[2], x$value[3], sum(x$value[c(1:3)]))
  alo.IB <- c(x$value[4], x$value[5], x$value[6], sum(x$value[c(4:6)]))
  alo.RI <- c(x$value[7], x$value[8], x$value[9], sum(x$value[c(7:9)]))
  alo.tot <- c(sum(x$value[c(1,4,7)]), sum(x$value[c(2,5,8)]), sum(x$value[c(3,6,9)]), sum(sum(x$value[c(1:3)]), sum(x$value[c(4:6)]), sum(x$value[c(7:9)])))
  alo.df <- rbind(alo.pib, alo.IB, alo.RI, alo.tot)
  }))
  rownames(contentDF) <- c("Present in at least one of each inbred panel", "Present in at least one classic inbred strain",
                           "Present in at least one recombinant inbred strain", "Total for at least one grouping",
                           "Present in all strains", "Present in all calssic inbred strains", "Present in all recombinant inbred strains",
                           "Total in all grouping")
  colnames(contentDF) <- c("Known transcript of known gene", "Novel transcript of known gene", "Novel transcript of novel gene", "Total")
  return(contentDF)
}