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

# extract 3006 SNPs and format it into a region file
fname1 = 'Del19_SGS_outlier_3006.list'  
dat_sites = read.delim(fname1, header = FALSE, sep='\t')
range = paste0(dat_sites$V1, ":",dat_sites$V2-1, "-", dat_sites$V2)
write.table(range, "./outlier_mafs/outlier_3006_list.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)

############################################################
# extract LD pruning SNPs and format it into angsd site file
setwd("~/Documents/Ryan_workplace/DelBay_adult/15_ngsLD")
fname1 = 'WILD_pruned_450776.snplist.txt'  
dat_sites = read.delim(fname1, header = FALSE, sep='_')
dat_sites$id = paste0(dat_sites$V1,'_',"0", dat_sites$V2,"_",dat_sites$V3)
colnames(dat_sites) <- c("chr1","chr2","pos", "id")
# specify the number of decimal point
options("scipen"=100, "digits"=6)


fname = "WILD_all_1680985.snplist.txt"  
df = read.delim(fname, header = FALSE, sep='\t')
df$id = paste0(df$V1,'_',df$V2)
pruned_list <- as.data.frame(cbind(df$V1[df$id %in% dat_sites$id], 
                                    df$V2[df$id %in% dat_sites$id], 
                                    df$V3[df$id %in% dat_sites$id], 
                                    df$V4[df$id %in% dat_sites$id]))
write.table(pruned_list,"WILD_pruned_723274.snplist.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)


############################################################
# 04132021 extract LD pruning SNPs from ngsLD outputs and format it into BayPass format
setwd("~/Documents/Ryan_workplace/DelBay_adult/15_ngsLD")
fname1 = 'WILD_pruned_450776.snplist.txt'  
dat_sites = read.delim(fname1, header = FALSE, sep=':')
dat_sites$id = paste0(dat_sites$V1,'_',dat_sites$V2)
colnames(dat_sites) <- c("chr","pos", "id")
# specify the number of decimal point
options("scipen"=100, "digits"=6)


fname = "Del19_20_1x_share_snps.list"  
df = read.delim(fname, header = FALSE, sep='\t')
df$id = paste0(df$V1,'_',df$V2)
pruned_list <- as.data.frame(cbind(df$V1[df$id %in% dat_sites$id], 
                                   df$V2[df$id %in% dat_sites$id], 
                                   df$V3[df$id %in% dat_sites$id], 
                                   df$V4[df$id %in% dat_sites$id]))
write.table(pruned_list,"WILD_pruned_450776.snplist.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)


