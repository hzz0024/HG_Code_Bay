#########################
# process the mafs file #
#########################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/07_mafs_process/original_mafs/")
process_mafs <- function(fname){
  #fname = "CHR19_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs"
  dat = read.delim(fname, header = TRUE, sep='\t')
  dat$id <- paste0(dat$chromo,'_',dat$position)
  k = round(dat$knownEM*dat$nInd*2)
  message(paste0("number of SNPs in ", fname, " is ", length(dat$chromo)))
  message("maximum allele frequency is ", max(dat$knownEM))
  message("allele frequency larger than 0.5 is ", length(dat$knownEM[which(dat$knownEM > 0.5)]))
  message(paste0("number of counts range from ", min(k), " to ", max(k)))
  hist(dat$knownEM, main=paste0("Allele frequency in ", fname),
       xlab="Allele frequencies" )
  return(dat$id)
  #dat[which(dat$knownEM <0.95),]
}

ARN <- process_mafs("ARN_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
CHR19 <- process_mafs("CHR19_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
CHR20 <- process_mafs("CHR20_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
COH <- process_mafs("COH_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
HC <- process_mafs("HC_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
NB <- process_mafs("NB_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
REF19 <- process_mafs("REF19_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
REF20 <- process_mafs("REF20_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
SR <- process_mafs("SR_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")

#########################
# extract common shared #
#########################

common_id <- Reduce(intersect, list(ARN, CHR19, CHR20, COH, HC, NB, REF19, REF20, SR))
length(common_id)
# [1] 893419

#########################
# process the ps file   #
#########################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/07_mafs_process/indep_SGS_angular_ps_shared_snp")
process_ps <- function(fname){
  #fname = "ps_Del19_challenge.txt"
  dat = read.delim(fname, header = FALSE, sep='\t')
  dat$id <- paste0(dat$V1,'_',dat$V2)
  message(paste0("number of SNPs in ", fname, " is ", length(dat$V1)))
  message("number of shared SNP is ", length(common_id))
  hist(dat$V5, main=paste0("Transformed D in ", fname),
       xlab="Transformed D" )
  extract <- dat[dat$id %in% common_id, ]
  message(length(extract$id))
  write.table(extract[,1:7], file = paste0("./shared/",fname), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

process_ps("ps_ARN_COH.txt")
process_ps("ps_Del19_challenge.txt")
process_ps("ps_Del19_HC_NB.txt")
process_ps("ps_Del20_challenge.txt")
process_ps("ps_SR_REF19.txt")
