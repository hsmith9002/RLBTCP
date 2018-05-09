#' Reading in RSEM output - Isoform
#' 
#' This function reads in the RSEM output for each sample and merges the data into one count matrix with rows as transcript IDs and columns as sample names
#' @param dir Global directory. Often points to project start folder. This will be the Parent folder for the folder containing lists of file names for each rsem batch output
#' @param dirs2 Directory pointing parent folder where RSEM output is
#' @param tissue Character value defining name of tissue in file path
#' @param N Numerical vector of batch numbers
#' @param qtype What quantitation type should the function read in. Can read counts or TPM
#' @param level Should the function generate expression data at the gene level or transcript level
#' @param filename file name for output as RData file
#' @keywords rsem
#' @export
#' @examples 
#' rsem.read.iso()

rsem.read.iso <- function(dir, dirs2, tissue, N, qtype, level, filename) {
  for(h in unique(N)) {
    ## Define "not in" function
    '%!in%' <- function(x,y)!('%in%'(x,y))
    
    ## Read in batch file list and extract unique strain names with replicate number attached
    flist <- as.character(read.table(file = paste(dir, "metadata/", tissue, "/", "RSEM.", "batch", h, ".filelist.txt", sep = ""), header = F)$V1)
    sampleID <- unique(unlist(lapply(strsplit(flist,split=".",fixed=TRUE),function(a) a[1])))
    sampleID <- sampleID[which(sampleID %!in% c("prep", paste("prepL", h, sep = ""), "rsem", "rsemb1r1", "RSEM", paste("rsemp", c(1:10), sep = ""), paste("rsemb", c(1:10), sep = ""), paste("rsem", c(1:10), sep = ""), "rsem_redo", "new_rsem", paste("rsemL", N, sep = "")))]
    fcheck <- rep(NA, length(sampleID))
    for(v in 1:length(sampleID)){
      if(c(paste(sampleID[v], ".", "isoforms.results", sep = "")) %in% flist) {fcheck[v] = 1.0} else {fcheck[v] = 0.0}
    }
    fcheck <- ifelse(fcheck == 1, TRUE, FALSE)
    if(sum(fcheck) != length(sampleID)) {stop(paste('There is no .results file present in directory for one of the strains: ', sampleID[!fcheck], "\n", sep = ""))}
  }
  
  for(h in unique(N)) {
    ## Define "not in" function
    '%!in%' <- function(x,y)!('%in%'(x,y))
    
    ## Read in batch file list and extract unique strain names with replicate number attached
    flist <- as.character(read.table(file = paste(dir, "metadata/", tissue, "/", "RSEM.", "batch", h, ".filelist.txt", sep = ""), header = F)$V1)
    sampleID <- unique(unlist(lapply(strsplit(flist,split=".",fixed=TRUE),function(a) a[1])))
    sampleID <- sampleID[which(sampleID %!in% c("prep", paste("prepL", h, sep = ""), "rsem", "rsemb1r1", "RSEM", paste("rsemp", c(1:10), sep = ""), paste("rsem", c(1:10), sep = ""), "rsem_redo", "new_rsem", paste("rsemL", N, sep = "")))]
    if(level == "transcript"){
      if(qtype == "counts") {
        ## Read in transcript level data
        for(i in sampleID){
          x = read.table(file=paste(dirs2, "batch", h, "/RI.reconst.v1/",i,".isoforms.results",sep=""),sep="\t",header=TRUE)
          y = data.frame(transcript_id = x$transcript_id, i=x$expected_count)
          colnames(y)[2] = i
          if(i==sampleID[1]) cntsT = y
          if(i!=sampleID[1]) cntsT = merge(cntsT,y,by="transcript_id")
        }
        colnames(cntsT)[-1] <- paste(colnames(cntsT)[-1], "_batch", h, sep = "")
        
        ## Merge datasets
        if(h==N[1]) cntsMerged = cntsT
        if(h!=N[1]) cntsMerged = merge(cntsMerged,cntsT,by="transcript_id")
        
        ## save cnts
        save(list = c("cntsMerged"), file = paste(dir, "Transcript_Reconstruction/", tissue, "/quantitation/",  filename, ".Rdata", sep = ""))
      } else{
        ## Read in transcript level data
        for(i in sampleID){
          x = read.table(file=paste(dirs2, "batch", h, "/RI.reconst.v1/",i,".isoforms.results",sep=""),sep="\t",header=TRUE)
          y = data.frame(transcript_id = x$transcript_id, i=x$TPM)
          colnames(y)[2] = i
          if(i==sampleID[1]) cntsT = y
          if(i!=sampleID[1]) cntsT = merge(cntsT,y,by="transcript_id")
        }
        colnames(cntsT)[-1] <- paste(colnames(cntsT)[-1], "_batch", h, sep = "")
        
        ## Merge datasets
        if(h==N[1]) cntsMerged = cntsT
        if(h!=N[1]) cntsMerged = merge(cntsMerged,cntsT,by="transcript_id")
        
        ## save cnts
        save(list = c("cntsMerged"), file = paste(dir, "Transcript_Reconstruction/", tissue, "/quantitation/",  filename, ".Rdata", sep = ""))
      }
    }
    else {
      if(qtype == "counts") {
        ## Read in transcript level data
        for(i in sampleID){
          x = read.table(file=paste(dirs2, "batch", h, "/RI.reconst.v1/",i,".genes.results",sep=""),sep="\t",header=TRUE)
          y = data.frame(gene_id = x$gene_id, i=x$expected_count)
          colnames(y)[2] = i
          if(i==sampleID[1]) cntsG = y
          if(i!=sampleID[1]) cntsG = merge(cntsG,y,by="gene_id")
        }
        colnames(cntsG)[-1] <- paste(colnames(cntsG)[-1], "_batch", h, sep = "")
        
        ## Merge datasets
        if(h==N[1]) cntsMerged.gene = cntsG
        if(h!=N[1]) cntsMerged.gene = merge(cntsMerged.gene,cntsG,by="gene_id")
        
        ## save cnts
        save(list = c("cntsMerged.gene"), file = paste(dir, "Transcript_Reconstruction/", tissue, "/quantitation/",  filename, ".Rdata", sep = ""))
      } else{
        ## Read in transcript level data
        for(i in sampleID){
          x = read.table(file=paste(dirs2, "batch", h, "/RI.reconst.v1/",i,".genes.results",sep=""),sep="\t",header=TRUE)
          y = data.frame(gene_id = x$gene_id, i=x$TPM)
          colnames(y)[2] = i
          if(i==sampleID[1]) cntsG = y
          if(i!=sampleID[1]) cntsG = merge(cntsG,y,by="gene_id")
        }
        colnames(cntsG)[-1] <- paste(colnames(cntsG)[-1], "_batch", h, sep = "")
        
        ## Merge datasets
        if(h==N[1]) cntsMerged.gene = cntsG
        if(h!=N[1]) cntsMerged.gene = merge(cntsMerged.gene,cntsG,by="gene_id")
        
        ## save cnts
        save(list = c("cntsMerged.gene"), file = paste(dir, "Transcript_Reconstruction/", tissue, "/quantitation/",  filename, ".Rdata", sep = ""))
      }
    }
  }
}