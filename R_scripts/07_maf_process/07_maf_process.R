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

extract_id <- function(fname){
  dat = read.delim(fname, header = TRUE, sep='\t')
  dat$id <- paste0(dat$chromo,'_',dat$position)
  extract <- dat[dat$id %in% common_id, ]
  message(length(extract$id))
  write.table(extract, file = paste0("../shared/",fname), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  #return(extract)
} 

ARN <- extract_id("ARN_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
CHR19 <- extract_id("CHR19_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
CHR20 <- extract_id("CHR20_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
COH <- extract_id("COH_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
HC <- extract_id("HC_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
NB <- extract_id("NB_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
REF19 <- extract_id("REF19_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
REF20 <- extract_id("REF20_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
SR <- extract_id("SR_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
# 
# > CHR19 <- process_mafs("./shared/CHR19_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
# number of SNPs in ./shared/CHR19_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs is 893419
# maximum allele frequency is 0.951551
# allele frequency larger than 0.5 is 36231
# number of counts range from 1 to 77
# > CHR20 <- process_mafs("./shared/CHR20_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs")
# number of SNPs in ./shared/CHR20_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs is 893419
# maximum allele frequency is 0.695183
# allele frequency larger than 0.5 is 20103
# number of counts range from 1 to 65
# 
# source("manhattan.R")
# dat = read.delim("CHR19_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs", header = TRUE, sep='\t')
# dat$SNP <- paste0(dat$chromo,'_',dat$position)
# manhattan(chr="chromo",bp="position",p="knownEM", snp = "SNP", dat, logp=FALSE, cex.axis = 1.2, ylim = c(0, 1),
#           col=c("grey50","black"),genomewideline=F, suggestiveline=F,
#           ylab=expression(log~FDR), cex.lab=1.5) 
# 
