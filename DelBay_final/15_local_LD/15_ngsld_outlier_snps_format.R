
#######################################################################
##########  format SGS and permutation shared SNP for local LD ########
#######################################################################
distance = 250

# SGS FDR0.1 outliers
shared_format = function(SGS_outliers, contrast){
  #contrast = "19_Sur_Ref"
  SGS_outlier = read.delim(SGS_outliers, header = TRUE, sep='\t') 
  SGS_outlier_id <- SGS_outlier$id
  # permutation raw p-values 
  df = read.delim("out_selection_index.txt", header = TRUE, sep='\t')
  df$SNP <- paste0(df$chromo, "_", df$position)
  hist(df$pvalue)
  #outlier <- df
  outlier <- df[which(df$pvalue < 0.005),]
  message("number of SNPs from permutation is ", dim(outlier)[1])
  # No.shared outliers
  shared <- SGS_outlier[which(SGS_outlier$id %in% outlier$SNP), ]
  message("number of SNPs shared is ", dim(shared)[1])
  message("number of SNPs with positive delta_p is ", length(shared$deltap[which(shared$deltap>=0)]))
  message("number of SNPs with negative delta_p is ", length(shared$deltap[which(shared$deltap< 0)]))
  shared_sort = shared[with(shared, order(chr, pos)),]
  #shared_sort[which(shared_sort$id == 'NC_035784.1_16552716'),] #check if chr5_16552716 is there
  system("mkdir Shared")
  write.table(shared_sort, paste0("./Shared/", contrast, "_shared_outlier.txt"), row.names=F, col.names = T, quote=F, sep="\t")
  
  # Step 1: format the shared outliers and convert it to bed format (must sort the position for ngsLD running)
  bed_list <- paste0(shared_sort$chr, "\t" , shared_sort$pos-distance, "\t", shared_sort$pos+distance)
  write.table(bed_list, paste0("./Shared/", contrast, "_shared_outlier.bed"), row.names=F, col.names = F, quote=F, sep="\t")
  
  # Step 2: merge the bed files
  bedtools = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/bedtools2/bin/bedtools";
  system(paste(bedtools, " merge -i ", "./Shared/", contrast, "_shared_outlier.bed > ", "./Shared/", contrast, "_shared_outlier.bed.merged.txt",  sep=""))

  # Step 3: convert the bed format to rf input for Angsd
  dat = read.delim(paste0("./Shared/", contrast, "_shared_outlier.bed.merged.txt"), header = FALSE, sep='\t')
  angsd_list <- paste0(dat$V1, ":" , dat$V2, "-", dat$V3)
  write.table(angsd_list, paste0("./Shared/", contrast, ".rf.txt"), row.names=F, col.names = F, quote=F, sep="\t")
}

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/test/")
shared_format("test_outlier.list","test")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Challenge_19")
shared_format("19_SGS_Sur_Ref_ps_outlier.list","19_Sur_Ref")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Challenge_20/")
shared_format("20_SGS_Sur_Ref_ps_outlier.list","20_Sur_Ref")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Wild_18/")
shared_format("18_SGS_HC_NB_ps_outlier.list","18_HC_NB")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Wild_19/")
shared_format("19_SGS_HC_NB_ps_outlier.list","19_HC_NB")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Wild_21/")
shared_format("21_SGS_HC_NB_ps_outlier.list","21_HC_NB")

# format for GEA outliers
bed_format = function(bed_file_name, contrast){
  #contrast = "CD5"
  bed_file = read.delim(bed_file_name, header = F, sep='\t') 
  bed_file_id <- paste0(bed_file$V1, "_", bed_file$V2)
  message("number of SNPs in ", contrast, " is ", length(bed_file_id))
  bed_file_sort = bed_file[with(bed_file, order(V1, V2)),]
  #shared_sort[which(shared_sort$id == 'NC_035784.1_16552716'),] #check if chr5_16552716 is there
  write.table(bed_file_sort, paste0(contrast, "_outlier.txt"), row.names=F, col.names = F, quote=F, sep="\t")
  
  # Step 1: format the shared outliers and convert it to bed format (must sort the position for ngsLD running)
  bed_list <- paste0(bed_file_sort$V1, "\t" , bed_file_sort$V2-distance, "\t", bed_file_sort$V3+distance)
  write.table(bed_list, paste0(contrast, "_outlier.bed"), row.names=F, col.names = F, quote=F, sep="\t")
  
  # Step 2: merge the bed files
  bedtools = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/bedtools2/bin/bedtools";
  system(paste(bedtools, " merge -i ", contrast, "_outlier.bed > ", contrast, "_outlier.bed.merged.txt",  sep=""))
  
  # Step 3: convert the bed format to rf input for Angsd
  dat = read.delim(paste0(contrast, "_outlier.bed.merged.txt"), header = FALSE, sep='\t')
  angsd_list <- paste0(dat$V1, ":" , dat$V2, "-", dat$V3)
  write.table(angsd_list, paste0(contrast, ".outlier.rf.txt"), row.names=F, col.names = F, quote=F, sep="\t")
}

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/CD5/")
bed_format("CD5_shared_lfmm_rda.bed","CD5")

########################################################
# Step 4: extract the rf part based on chromosomes
########################################################
# path /Users/HG/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/shared/ngsld_outlier_block
# path /Users/HG/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/shared/random_theta_distribution
# for i in {0..9}; do
# #for i in NC_035780.1 NC_035781.1 NC_035782.1 NC_035783.1 NC_035784.1 NC_035785.1 NC_035786.1 NC_035787.1 NC_035788.1 NC_035789.1; do
# j=$[i+1]
# grep 'NC_03578'$i'.1' ps_Del20_challenge.rf.txt > 'ps_Del20_challenge.rf.chr'$j'.txt'
# done
# 
# for i in {0..9}; do
# #for i in NC_035780.1 NC_035781.1 NC_035782.1 NC_035783.1 NC_035784.1 NC_035785.1 NC_035786.1 NC_035787.1 NC_035788.1 NC_035789.1; do
# j=$[i+1]
# grep 'NC_03578'$i'.1' ps_Del19_HC_NB.rf.txt > 'ps_Del19_HC_NB.rf.chr'$j'.txt'
# done
# 
# for i in {0..9}; do
# #for i in NC_035780.1 NC_035781.1 NC_035782.1 NC_035783.1 NC_035784.1 NC_035785.1 NC_035786.1 NC_035787.1 NC_035788.1 NC_035789.1; do
# j=$[i+1]
# grep 'NC_03578'$i'.1' ps_Del19_challenge.rf.txt > 'ps_Del19_challenge.rf.chr'$j'.txt'
# done
# 
# for i in {0..9}; do
# #for i in NC_035780.1 NC_035781.1 NC_035782.1 NC_035783.1 NC_035784.1 NC_035785.1 NC_035786.1 NC_035787.1 NC_035788.1 NC_035789.1; do
# j=$[i+1]
# grep 'NC_03578'$i'.1' ps_Del19_ARN_COH.rf.txt > 'ps_Del19_ARN_COH.rf.chr'$j'.txt'
# done
