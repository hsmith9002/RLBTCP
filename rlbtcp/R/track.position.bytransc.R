#' Find Chromosome position for use in UCSC track visualizer with geneID
#' 
#' Brief Description: This function does the same thing as track.position.table() but allows the user to search for transcripts from a
#' gene ID chosen by the user a priori. These transcript lists go into the the function in place of the number of exons for the track.position.table function.
#' @param eptOut output object from exon.tpt.table()
#' @param tlist1 transcript list 1
#' @param tlist2 transcript list 2
#' @param tlist3 transcript list 3
#' @keywords UCSC Genome browser custome track
#' @export
#' @examples 
#' track.position.bytransc(eptOut, tlist1, tlist2, tlist2)

track.position.bytransc <- function(eptOut, tlist1, tlist2, tlist3) {
  ept.known <- eptOut$EPT.known
  ept.ntng <- eptOut$EPT.ntng
  ept.ntkg <- eptOut$EPT.ntkg
  
  ## Assign transcript summary tables to objects 
  tbl.known <- ept$table.known
  tbl.ntng <- ept$table.ntng
  tbl.ntkg <- ept$table.ntkg
  
  ## Create transcript name index for subsetting reference dataset
  fortrack.known.index <- ept.known[which(ept.known$transcript_id %in% tlist1), 1]
  #fortrack.known.index <- fortrack.known.index[!duplicated(fortrack.known.index$Exon_Count), 1]
  
  fortrack.ntkg.index <- ept.ntkg[which(ept.ntkg$transcript_id %in% tlist2), 1]
  #fortrack.ntkg.index <- fortrack.ntkg.index[!duplicated(fortrack.ntkg.index$Exon_Count), 1]
  
  fortrack.ntng.index <- ept.ntng[which(ept.ntng$transcript_id %in% tlist3), 1]
  #fortrack.ntng.index <- fortrack.ntng.index[!duplicated(fortrack.ntng.index$Exon_Count), 1]
  
  ## extract transcripts in index above from reference table
  fortrack.known <- tbl.known[which(tbl.known$transcript_id %in% fortrack.known.index), ]
  fortrack.ntng <- tbl.ntng[which(tbl.ntng$transcript_id %in% fortrack.ntng.index), ]
  fortrack.ntkg <- tbl.ntkg[which(tbl.ntkg$transcript_id %in% fortrack.ntkg.index), ]
  
  ## Separate by transcript ID and build transcript range dataframe to use as reference when looking at tracks
  fortrack.known.list <- split(fortrack.known, f = fortrack.known$transcript_id)
  known.chrnum <- sapply(fortrack.known.list, function(x) paste("chr", x[1, 4], sep = ""))
  known.maxpos <- sapply(fortrack.known.list, function(x) max(x$start, na.rm=TRUE))
  known.minpos <- sapply(fortrack.known.list, function(x) min(x$end, na.rm=TRUE))
  known.range <- rbind(known.chrnum, known.maxpos, known.minpos)
  fr.known <- c(min(known.range[3, ]), max(known.range[2, ]))
  
  fortrack.ntng.list <- split(fortrack.ntng, f = fortrack.ntng$transcript_id)
  ntng.chrnum <- sapply(fortrack.ntng.list, function(x) paste("chr", x[1, 4], sep = ""))
  ntng.maxpos <- sapply(fortrack.ntng.list, function(x) max(x$start, na.rm=TRUE))
  ntng.minpos <- sapply(fortrack.ntng.list, function(x) min(x$end, na.rm=TRUE))
  ntng.range <- rbind(ntng.chrnum, ntng.maxpos, ntng.minpos)
  fr.ntng <- c(min(ntng.range[3, ]), max(ntng.range[2, ]))
  
  fortrack.ntkg.list <- split(fortrack.ntkg, f = fortrack.ntkg$transcript_id)
  ntkg.chrnum <- sapply(fortrack.ntkg.list, function(x) paste("chr", x[1, 4], sep = ""))
  ntkg.maxpos <- sapply(fortrack.ntkg.list, function(x) max(x$start, na.rm=TRUE))
  ntkg.minpos <- sapply(fortrack.ntkg.list, function(x) min(x$end, na.rm=TRUE))
  ntkg.range <- rbind(ntkg.chrnum, ntkg.maxpos, ntkg.minpos)
  fr.ntkg <- c(min(ntkg.range[3, ]), max(ntkg.range[2, ]))
  return(list("known.range" = known.range, "R.known" = fr.known, "ntng.range" = ntng.range, "R.ntng" = fr.ntng, "ntkg.range" = ntkg.range, "R.ntkg" = fr.ntkg))
}