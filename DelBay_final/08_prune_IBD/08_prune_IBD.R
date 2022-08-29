########################################################
# Step 1: convert the snp list for Angsd
########################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/08_prune_IBD")

format_snp_list <- function(pname){
  #pname = "WILD.unlinked.id"
  snp_list = read.delim("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.txt", header = FALSE, sep='\t')
  dat = read.delim(pname, header = FALSE, sep=':')
  dat = dat[with(dat, order(V1, V2)),]
  angsd_list <- snp_list[which(paste0(snp_list$V1, "_",snp_list$V2) %in% paste0(dat$V1, "_", dat$V2)),]
  angsd_list = angsd_list[with(angsd_list, order(V1, V2)),]
  write.table(angsd_list, paste0(strsplit(pname, split = ".txt")[[1]][1], ".pruned.site.txt"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_snp_list("WILD.unlinked.id")
