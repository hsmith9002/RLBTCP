#' Summarize transcript content into a tabular format and perform QC checks
#' 
#' This function summarizes the transcript content and checks to make sure
#' each subcategory of transcripts is a true subset of the complete ALO set of transcripts.
#' @param subtype vector of present options as they are coded in the .results obects
#' @param dataList List .results datasets (i.e. known.results, ntkg.results, ntng.results)
#' @param refDat rbinded data of the 3 datasets listed above with the transcripts labeled "Not present" removed. This should have 114982 transcripts.
#' @param i index indicating whether or not to summarize ALO data or ALL data, i = 1 for ALO, and 2 for ALL
#' @export
#' @examples 
#' contentQC(subtype, dataList, refDat, i)

contentQC <- function(subtype, dataList, refDat, i) {
  st1 <- subtype[1]
  st2 <- subtype[2]
  st3 <- subtype[3]
  
  qclist <- do.call("rbind", llply(dataList, .fun = function(x) {
    known.pib <- x[which(x[,i] %in% c(st1)), ]
    pib.test <- sum(known.pib[,3] %in% refDat$transcript_id) == length(known.pib[,3])
    pib.len <- length(known.pib[,3])
    known.IB <- x[which(x[,i] %in% c(st2)), ]
    IB.test <- sum(known.IB[,3] %in% refDat$transcript_id) == length(known.IB[,3])
    IB.len <- length(known.IB[,3])
    known.RI <- x[which(x[,i] %in% c(st3)), ]
    RI.test <- sum(known.RI[,3] %in% refDat$transcript_id) == length(known.RI[,3])
    RI.len <- length(known.RI[,3])
    known.test <- sum(pib.test, IB.test, RI.test) == 3
    df <- c(pib.len, IB.len, RI.len, sum(pib.len, IB.len, RI.len), known.test)
  }))
  
  qcvect <- do.call("c", llply(dataList, .fun = function(x) {
    known.pib <- x[which(x[,i] %in% c(st1)), ]
    pib.test <- sum(known.pib[,3] %in% refDat$transcript_id) == length(known.pib[,3])
    pib.len <- length(known.pib[,3])
    known.IB <- x[which(x[,i] %in% c(st2)), ]
    IB.test <- sum(known.IB[,3] %in% refDat$transcript_id) == length(known.IB[,3])
    IB.len <- length(known.IB[,3])
    known.RI <- x[which(x[,i] %in% c(st3)), ]
    RI.test <- sum(known.RI[,3] %in% refDat$transcript_id) == length(known.RI[,3])
    RI.len <- length(known.RI[,3])
    known.test <- sum(pib.test, IB.test, RI.test) == 3
    transcriptlist <- c(known.pib[,3], known.IB[,3], known.RI[,3])
  }))
  
  tots <- apply(qclist, 2, function(x) sum(x))
  df.out <- rbind(qclist, tots)
  df.out[,5] <- ifelse(df.out[,5] == c(1), "Yes", NA)
  rownames(df.out) <- c("Known transcript of known gene", "Novel transcript of known gene", "Novel transcript of novel gene", "Totals")
  colnames(df.out) <- c("Present in both", "Present in classic inbreds", "Present in recombinant inbreds", "Totals", "Subset of full dataset")
  transcriptList <- c()
  #names(qclist) <- c("KTKG", "NTKG", "NTNG")
  return(list(qcvect, df.out))
}