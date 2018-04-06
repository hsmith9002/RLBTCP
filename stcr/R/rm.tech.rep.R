#' Removing technical replicates
#' 
#' This function identifies technical replicates, finds the replicate with the highest reads count, and then remerges these replicates back into the original dataset. The resulting dataset contains expected counts for all strains minus technical replicates, and if a strain had technical replicates the one with the hiest read count is retained.
#' @param cntsMerged This is the original expected counts matrix generated from using the rsem.read.iso() function
#' @keywords technical replicate
#' @export
#' @examples 
#' rm.tech.rep(cntsMerged)

rm.tech.rep <- function(cntsMerged) {
  df <- cntsMerged
  ## Extract strain name and biological replicate nuber and paste
  names <-  unlist(lapply(strsplit(names(df),split="_",fixed=TRUE),function(a) paste(a[1], a[2], sep = "_")))
  ## generate vector of names that have identical strain and biological replicate numbers. These are technical replicates
  rep.index <- names[duplicated(names)]
  ## Generate dataframe that contains technical replicates only
  df.w.reps <- df[, grepl(paste(rep.index, collapse = "|"), names(df))]
  ## Generate dataframe that does not contain technical replicates
  df.wo.reps <- df[, !grepl(paste(rep.index, collapse = "|"), names(df))]
  ## identify the technical replicates with the highest overall read count
  sums.vect <- sort(apply(df.w.reps, MARGIN = 2, FUN = sum), decreasing = TRUE)
  sums.vect <- as.matrix(sums.vect)
  ## Extract strain name
  strain <- unlist(lapply(strsplit(rownames(sums.vect),split="_",fixed=TRUE),function(a) paste(a[1], a[2], sep = "_")))
  ## Extract batch name
  batch <- unlist(lapply(strsplit(rownames(sums.vect),split="_",fixed=TRUE),function(a) a[3]))
  ## put in dataframe
  tmp <- as.data.frame(cbind(sums.vect, strain, batch))
  ## Extract the sample with the highest read count by selecting its first occurance. This works because the datframe is sorted from highest to lowest, 
  #  so the first occurance of a strain name is its highest read count
  tmp2 <- tmp[!duplicated(tmp$strain), ]
  ## Extract names from dataframe to use as index for data with technical replicates
  new.index <- rownames(tmp2)
  ## Extract columns for technical replicates with highest read counts
  df.w.reps.new <- df.w.reps[, new.index]
  ## recombine these technical replicates with the dataframe that does not contain any technical replicates
  dfFinal <- as.data.frame(cbind(df.wo.reps, df.w.reps.new))
  return(dfFinal)
}
