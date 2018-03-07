#' Generate Present/Absent Table
#' 
#' This function identifies transcripts as present/absent, or other based on a user defined
#' criteria, and then assigns a value (1 =Present, 0=Absent, 0.5 =Other) on the 
#' strain level.
#' @param cnts expected count matrix where rows are transcript/gene ids and columns are strains. The first column should be a vector of transcript ids
#' @param datanames Character vector of unique strain names used in cnts data including transcript id col name
#' @param presThresh Upper count threshold to be considered present (default = 10)
#' @param absThresh Lower count threshold to be considered absent (default = 1)
#' @keywords Present/Absent
#' @export
#' @examples 
#' pres.abs.table(cnts, names, 10, 1)


pres.abs.table <- function(cnts, datanames, presThresh = 10, absThresh = 1) {
  df <- sapply(unique(datanames)[-1], function(x){
    cols <- grep(paste0('^', x, '[-_]+', '[1-9]+'), colnames(cnts))
    apply(as.matrix(cnts[, cols]), MARGIN =1, function(y) {
      if(all(y >= presThresh, na.rm = T)) return(1)
      if(all(y <absThresh, na.rm = T)) return(0)
      return(0.5)
    })
  })
  df <- as.data.frame(cbind(cnts[1], df))
  return(df)
}