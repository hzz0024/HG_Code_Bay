# set up the location for data running
setwd("~/Documents/Ryan_workplace/DelBay_adult/13_env_gen_association/key_snps")
# extract 3006 SNPs from the whole snp list
fname1 = 'Del19_SGS_outlier_3006.list'  
dat_sites = read.delim(fname1, header = FALSE, sep='\t')
dat_sites$id = paste0(dat_sites$V1,'_',dat_sites$V2)
colnames(dat_sites) <- c("chr","pos","major","minor", "id")
# specify the number of decimal point
options("scipen"=100, "digits"=6)

for (pop in c("HC", "ARN", "COH", "NB", "SR")){
  print(paste("Work on population:", pop))
  fname = paste0(pop,"_maf0.05_pctind0.7_maxdepth3.mafs")  
  df = read.delim(fname, header = TRUE, sep='\t')
  df$id = paste0(df$chromo,'_',df$position)
  outlier_list <- as.data.frame(cbind(df$chromo[df$id %in% dat_sites$id], 
                                      df$pos[df$id %in% dat_sites$id], 
                                      df$major[df$id %in% dat_sites$id], 
                                      df$minor[df$id %in% dat_sites$id],
                                      df$anc[df$id %in% dat_sites$id],
                                      df$knownEM[df$id %in% dat_sites$id],
                                      df$nInd[df$id %in% dat_sites$id]))
  colnames(outlier_list) <- c("chromo","position","major","minor", "anc", "knownEM","nInd")
  write.table(outlier_list,paste0("./outlier_mafs/", pop,"_maf0.05_pctind0.7_maxdepth3.mafs"), row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
}


