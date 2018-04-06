#' Find Chromosome position for use in UCSC track visualizer
#' 
#' Brief Description: This function takes as input the output of exon.tpt.table(), and creates a table of chromosome positions to be used
#' with the UCSC genome browser custome track vizualizer. You can look at 3 different transcripts based on number of exons chosen by the user.
#' @param eptOut output object from exon.tpt.table()
#' @param ec1 Number of exons associated with transcript 1
#' @param ec2 Number of exons associated with transcript 2
#' @param ec3 Number of exons associated with transcript 3
#' @keywords UCSC Genome browser custome track
#' @export
#' @examples 
#' track.position.table(eptOut, 25, 10, 5)

track.position.table <- function(eptOut, ec1 = 25, ec2 = 10, ec3 = 5) {
  ept.ntng <- eptOut$EPT.ntng
  ept.ntkg <- eptOut$EPT.ntkg
  
  ## Assign transcript summary tables to objects 
  tbl.ntng <- ept$table.ntng
  tbl.ntkg <- ept$table.ntkg
  
  ## Create transcript name index for subsetting reference dataset
  
  fortrack.ntkg.index <- ept.ntkg[which(ept.ntkg$Exon_Count %in% c(max(ept.ntkg$Exon_Count), ec1, ec2, ec3)), ]
  fortrack.ntkg.index <- fortrack.ntkg.index[!duplicated(fortrack.ntkg.index$Exon_Count), 1]
  
  fortrack.ntng.index <- ept.ntng[which(ept.ntng$Exon_Count %in% c(max(ept.ntng$Exon_Count), ec1, ec2, ec3)), ]
  fortrack.ntng.index <- fortrack.ntng.index[!duplicated(fortrack.ntng.index$Exon_Count), 1]
  
  ## extract transcripts in index above from reference table
  fortrack.ntng <- tbl.ntng[which(tbl.ntng$transcript_id %in% fortrack.ntng.index), ]
  fortrack.ntkg <- tbl.ntkg[which(tbl.ntkg$transcript_id %in% fortrack.ntkg.index), ]
  
  ## Separate by transcript ID and build transcript range dataframe to use as reference when looking at tracks
  fortrack.ntng.list <- split(fortrack.ntng, f = fortrack.ntng$transcript_id)
  ntng.chrnum <- sapply(fortrack.ntng.list, function(x) paste("chr", x[1, 4], sep = ""))
  ntng.maxpos <- sapply(fortrack.ntng.list, function(x) max(x$start, na.rm=TRUE))
  ntng.minpos <- sapply(fortrack.ntng.list, function(x) min(x$end, na.rm=TRUE))
  ntng.range <- rbind(ntng.chrnum, ntng.maxpos, ntng.minpos)
  
  fortrack.ntkg.list <- split(fortrack.ntkg, f = fortrack.ntkg$transcript_id)
  ntkg.chrnum <- sapply(fortrack.ntkg.list, function(x) paste("chr", x[1, 4], sep = ""))
  ntkg.maxpos <- sapply(fortrack.ntkg.list, function(x) max(x$start, na.rm=TRUE))
  ntkg.minpos <- sapply(fortrack.ntkg.list, function(x) min(x$end, na.rm=TRUE))
  ntkg.range <- rbind(ntkg.chrnum, ntkg.maxpos, ntkg.minpos)
  return(list("ntng.range" = ntng.range, "ntkg.range" = ntkg.range))
}