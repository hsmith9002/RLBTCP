#' Transcript Categorization With Sampling
#' 
#' This function first randomly samples a given number of strains from each panel, then assigns codes to transcripts based on wether they are in 1 of 4 categories.
#' If in at least 1 strain in RI (Known and novel datasets
#' analyzed separatly), then return 1, else return 5.
#' If in at least 1 strain in IB (Known and novel datasets
#' analyzed separatly), then return 1, else return 10. This
#' will cretae 2 columns to be summed. The sum of the columns
#' has meaning:
#' 2 = present in both data sets
#' 6 = present in IB and NOT in RI
#' 11 = present in RI and NOT in IB
#' 15 = not in either data set
#' and for the in all subgroup
#' 2 = present in all strains
#' 6 = present in all IB
#' 11 = present in all RI
#' 15 = not no strains
#' @param data The present absent tabel created by pres.abs.table() function
#' @param RI.num Number of strains to sample from recombinant inbred panel
#' @param IB.num Number of strains to sample from calsic inbred panel
#' @export
#' @examples 
#' sep.panel.sample(PA.table, 5, 5)

sep.panel.sample <- function(data, RI.num, IB.num) {
  x <- c("HXB", "BXH","BN", "SHR") # Define RI names
  a <- c("ENSRNOT") # Define prefix that denotes known transcripts
  
  ## Create RI only dataset
  y <- grepl(paste(x, collapse = "|"), colnames(data))
  RI.tmp <- data[, y]
  RI.tmp <- as.data.frame(cbind(RI.tmp, data$type))
  RI <- RI.tmp[, colnames(RI.tmp) != "SHRSP"]
  rownames(RI) <- data$transcript_id
  colnames(RI)[length(colnames(RI))] <- "type"
  sample_index.RI <- c(sample(c(1:32), RI.num), 33)
  RI.sample <- RI[, sample_index.RI]
  
  ## Create Classic IB only dataset
  z <- !grepl(paste(x, collapse = "|"), colnames(data))
  IB.temp <- data[, z]
  IB <- as.data.frame(cbind(IB.temp[, -1], RI.tmp[, colnames(RI.tmp) == "SHRSP"]))
  colnames(IB)[colnames(IB) == "RI.tmp[, colnames(RI.tmp) == \"SHRSP\"]"] <- "SHRSP"
  rownames(IB) <- data$transcript_id
  sample_index.IB <- c(sample(c(1:10, 12), IB.num), 11)
  IB.sampled <- IB[, sample_index.IB]
  ## This removes the column labeled type
  RI.ktkg <- RI.sample[which(RI.sample$type %in% "KTKG"), -which(names(RI.sample) %in% "type")]
  RI.ntkg <- RI.sample[which(RI.sample$type %in% "NTKG"), -which(names(RI.sample) %in% "type")]
  RI.ntng <- RI.sample[which(RI.sample$type %in% "NTNG"), -which(names(RI.sample) %in% "type")]
  IB.ktkg <- IB.sampled[which(IB.sampled$type %in% "KTKG"), -which(names(IB.sampled) %in% "type")] 
  IB.ntkg <- IB.sampled[which(IB.sampled$type %in% "NTKG"), -which(names(IB.sampled) %in% "type")]
  IB.ntng <- IB.sampled[which(IB.sampled$type %in% "NTNG"), -which(names(IB.sampled) %in% "type")]
  
  ## separate each dataset into a dataset of known transcripts, and a dataset of novel transcripts
  #b <- grepl(a, data$transcript_id)
  #RI.KN <- RI[b, ]
  #RI.NV <- RI[!b, ]
  #IB.KN <- IB[b, ]
  #IB.NV <- IB[!b, ]
  
  ## Recode genes based on there location in each dataset
  ###############################################################
  # If in at least 1 strain in RI (Known and novel datasets
  # analyzed separatly), then return 1, else return 5.
  # If in at least 1 strain in IB (Known and novel datasets
  # analyzed separatly), then return 1, else return 10. This
  # will cretae 2 columns to be summed. The sum of the columns
  # has meaning:
  #     2 = present in both data sets
  #     6 = present in IB and NOT in RI
  #     11 = present in RI and NOT in IB
  #     15 = not in either data set
  ###############################################################
  RIKN.result <- apply(RI.ktkg, MARGIN = 1, FUN = function(x){
    if(any(x == 1.0, na.rm = T)) return(1)
    else return(5)
  })
  IBKN.result <- apply(IB.ktkg, MARGIN = 1, FUN = function(x){
    if(any(x == 1.0, na.rm = T)) return(1)
    else return(10)
  })
  RINTKG.result <- apply(RI.ntkg, MARGIN = 1, FUN = function(x){
    if(any(x == 1.0, na.rm = T)) return(1)
    else return(5)
  })
  RINTNG.result <- apply(RI.ntng, MARGIN = 1, FUN = function(x){
    if(any(x == 1.0, na.rm = T)) return(1)
    else return(5)
  })
  IBNTKG.result <- apply(IB.ntkg, MARGIN = 1, FUN = function(x){
    if(any(x == 1.0, na.rm = T)) return(1)
    else return(10)
  })
  IBNTNG.result <- apply(IB.ntng, MARGIN = 1, FUN = function(x){
    if(any(x == 1.0, na.rm = T)) return(1)
    else return(10)
  })
  
  ## Recode genes based on there location in each dataset
  ###############################################################
  # If in all strains in RI (Known and novel datasets
  # analyzed separatly), then return 1, else return 5.
  # If in all strains in IB (Known and novel datasets
  # analyzed separatly), then return 1, else return 10. This
  # will cretae 2 columns to be summed. The sum of the columns
  # has meaning:
  #     2 = present in all strains
  #     6 = present in all IB
  #     11 = present in all RI
  #     15 = not no strains
  ###############################################################
  RIKN.result2 <- apply(RI.ktkg, MARGIN = 1, FUN = function(x){
    if(all(x == 1.0, na.rm = T)) return(1)
    else return(5)
  })
  IBKN.result2 <- apply(IB.ktkg, MARGIN = 1, FUN = function(x){
    if(all(x == 1.0, na.rm = T)) return(1)
    else return(10)
  })
  RINTKG.result2 <- apply(RI.ntkg, MARGIN = 1, FUN = function(x){
    if(all(x == 1.0, na.rm = T)) return(1)
    else return(5)
  })
  RINTNG.result2 <- apply(RI.ntng, MARGIN = 1, FUN = function(x){
    if(all(x == 1.0, na.rm = T)) return(1)
    else return(5)
  })
  IBNTKG.result2 <- apply(IB.ntkg, MARGIN = 1, FUN = function(x){
    if(all(x == 1.0, na.rm = T)) return(1)
    else return(10)
  })
  IBNTNG.result2 <- apply(IB.ntng, MARGIN = 1, FUN = function(x){
    if(all(x == 1.0, na.rm = T)) return(1)
    else return(10)
  })
  
  ## Put new coded data in new dataframes
  knrslt.tmp <- as.data.frame(cbind(RIKN.result, IBKN.result))
  ntkgrslt.tmp <- as.data.frame(cbind(RINTKG.result, IBNTKG.result))
  ntngrslt.tmp <- as.data.frame(cbind(RINTNG.result, IBNTNG.result))
  knrslt.tmp2 <- as.data.frame(cbind(RIKN.result2, IBKN.result2))
  ntkgrslt.tmp2 <- as.data.frame(cbind(RINTKG.result2, IBNTKG.result2))
  ntngrslt.tmp2 <- as.data.frame(cbind(RINTNG.result2, IBNTNG.result2))
  kn.atleastsum <- rowSums(knrslt.tmp)
  kn.allsum <- rowSums(knrslt.tmp2)
  ntkg.atleastsum <- rowSums(ntkgrslt.tmp)
  ntng.atleastsum <- rowSums(ntngrslt.tmp)
  ntkg.allsum <- rowSums(ntkgrslt.tmp2)
  ntng.allsum <- rowSums(ntngrslt.tmp2)
  knrslt <- as.data.frame(cbind(knrslt.tmp, knrslt.tmp2, kn.atleastsum, kn.allsum))
  colnames(knrslt)[c(5,6)] <- c("RS.atleast1", "RS.all")
  ntkgrslt <- as.data.frame(cbind(ntkgrslt.tmp, ntkgrslt.tmp2, ntkg.atleastsum, ntkg.allsum))
  colnames(ntkgrslt)[c(5,6)] <- c("RS.atleast1", "RS.all")
  ntngrslt <- as.data.frame(cbind(ntngrslt.tmp, ntngrslt.tmp2, ntng.atleastsum, ntng.allsum))
  colnames(ntngrslt)[c(5,6)] <- c("RS.atleast1", "RS.all")
  return(list(RI_DF = RI.sample, InbredDF = IB.sampled, KN_Results =knrslt, NTKG_Results = ntkgrslt, NTNG_Results = ntngrslt))
}