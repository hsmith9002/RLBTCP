#' Calculate Exons Per Transcript
#' 
#' Brief Description: This function takes as input the cleaned
#' chromosome file and the data frames created from the merged
#' liver gtf, which have been subsetted based on transcripts 
#' present in the known, ntng or ntkg data. It then creates a reference
#' table of transcript names and exons, and calculates the 
#' number of exons/transcript.
#' @param ST_df.known R version of gtf subsetted for known transcripts
#' @param ST_df.ntkg R version of gtf subsetted for novel transcripts of known genes
#' @param ST_df.ntng R version of gtf subsetted for novel transcripts of novel genes
#' @param directory directory pointing to chromosome sizes file
#' @param chrSizes name of chromosome sizes file
#' @param tissue name of tissue (this has to do with the way my files are structured)
#' @keywords gtf chromosome sizes
#' @export
#' @examples 
#' exon.tpt.table(ST_df.known = NULL, ST_df.ntng = st.ntkg, ST_df.ntkg = st.ntng, directory = "~/genomedat/, chrSize, tissue = "liver")

exon.tpt.table <- function(ST_df.known = NULL, ST_df.ntng, ST_df.ntkg, directory, chrSize) {
  ###Create workspace and load packages###
  options(stringsAsFactors = FALSE)
  library(rtracklayer) 
  library(GenomicRanges) 
  library(dplyr)
  
  ##Import cleaned Chromosome sizes dataset
  dir <- directory #define file directory
  chrSize <- chrSize
  chro.clean <- read.table(paste(dir, tissue, "/quantitation/", chrSize, sep = ""), header = FALSE, sep = "\t")
  colnames(chro.clean) <- c("Chromosome", "Length")
  known_df <- ST_df.known
  ntng_df <- ST_df.ntng
  ntkg_df <- ST_df.ntkg
  
  #Create a table for counting exons per transcript - NTNG
  tctexonTable.known <- as.data.frame.matrix(table(ST_df.known$gene_id, ST_df.known$type))
  tctexonTable.known$gene_id <- rownames(tctexonTable.known)
  tbldfMerge.known <- merge(tctexonTable.known, ST_df.known, by = "gene_id")
  #Create a table for counting exons per transcript - NTNG
  tctexonTable.ntng <- as.data.frame.matrix(table(ST_df.ntng$gene_id, ST_df.ntng$type))
  tctexonTable.ntng$gene_id <- rownames(tctexonTable.ntng)
  tbldfMerge.ntng <- merge(tctexonTable.ntng, ST_df.ntng, by = "gene_id")
  #Create a table for counting exons per transcript - NTKG
  tctexonTable.ntkg <- as.data.frame.matrix(table(ST_df.ntkg$gene_id, ST_df.ntkg$type))
  tctexonTable.ntkg$gene_id <- rownames(tctexonTable.ntkg)
  tbldfMerge.ntkg <- merge(tctexonTable.ntkg, ST_df.ntkg, by = "gene_id")
  
  ##Calculate the number of exons per transcript - NTNG
  exnPERtpt.known <- as.data.frame(table(tbldfMerge.known$transcript_id)) # this counts all transcript IDs assigned to both exons and the original transcript which is not an exon
  exnPERtpt.known$Freq <- exnPERtpt.known$Freq - 1 # this removes the count attributed to the transcript
  exnPERtpt.known$Var1 <- as.character(exnPERtpt.known$Var1)
  exnPERtpt.known <- exnPERtpt.known[order(exnPERtpt.known$Freq, decreasing = T), ]
  colnames(exnPERtpt.known) <- c("transcript_id", "Exon_Count")
  
  ##Calculate the number of exons per transcript - NTNG
  exnPERtpt.ntng <- as.data.frame(table(tbldfMerge.ntng$transcript_id)) # this counts all transcript IDs assigned to both exons and the original transcript which is not an exon
  exnPERtpt.ntng$Freq <- exnPERtpt.ntng$Freq - 1 # this removes the count attributed to the transcript
  exnPERtpt.ntng$Var1 <- as.character(exnPERtpt.ntng$Var1)
  exnPERtpt.ntng <- exnPERtpt.ntng[order(exnPERtpt.ntng$Freq, decreasing = T), ]
  colnames(exnPERtpt.ntng) <- c("transcript_id", "Exon_Count")
  
  ##Calculate the number of exons per transcript - NTKG
  exnPERtpt.ntkg <- as.data.frame(table(tbldfMerge.ntkg$transcript_id)) # this counts all transcript IDs assigned to both exons and the original transcript which is not an exon
  exnPERtpt.ntkg$Freq <- exnPERtpt.ntkg$Freq - 1 # this removes the count attributed to the transcript
  exnPERtpt.ntkg$Var1 <- as.character(exnPERtpt.ntkg$Var1)
  exnPERtpt.ntkg <- exnPERtpt.ntkg[order(exnPERtpt.ntkg$Freq, decreasing = T), ]
  colnames(exnPERtpt.ntkg) <- c("transcript_id", "Exon_Count")
  return(list("EPT.known" = exnPERtpt.known, "EPT.ntng" = exnPERtpt.ntng, "EPT.ntkg" = exnPERtpt.ntkg, "table.known" = tbldfMerge.known, "table.ntng" = tbldfMerge.ntng, "table.ntkg" = tbldfMerge.ntkg))
}