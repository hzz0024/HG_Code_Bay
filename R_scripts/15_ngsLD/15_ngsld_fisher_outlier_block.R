
################################################
#######  process fisher outlier for ngsLD ######
################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/fisher_outlier")
########################################################
# Step 1: format the ps file and convert it to bed format
########################################################
# must sort the position for ngsLD running
format_bed <- function(pname, distance){
  #pname = "format/REF19_CHR19_NB_HC_out_0.05_fish.txt"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat = dat[with(dat, order(V1, V2)),]
  bed_list <- paste0(dat$V1, "\t" , dat$V2, "\t", dat$V2+distance)
  print(head(bed_list))
  write.table(bed_list, paste0(pname, ".bed"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_bed("./format/REF19_CHR19_NB_HC_out_0.05_fish.txt", 250)

########################################################
# Step 2: using bedtools to merge intervals
########################################################
# path in /workdir/hz269/DelBay_all_angsd_final/15_LD_prunning/process_bed_file

# for f in *.bed; do
# echo $f
# cat $f | wc -l
# bedtools merge -i $f > $f'.merged.txt'
# cat $f'.merged.txt' | wc -l
# done

# ps_Del19_ARN_COH.txt.bed
# 1113
# 922
# ps_Del19_challenge.txt.bed
# 2185
# 1854
# ps_Del19_HC_NB.txt.bed
# 2073
# 1696
# ps_Del20_challenge.txt.bed
# 3117
# 2525

########################################################
# Step 3: convert the bed format to rf input for Angsd
########################################################

format_rf <- function(pname){
  dat = read.delim(pname, header = FALSE, sep='\t')
  angsd_list <- paste0(dat$V1, ":" , dat$V2, "-", dat$V3)
  write.table(angsd_list, paste0(strsplit(pname, split = ".txt")[[1]][1], ".rf.txt"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_rf("format/REF19_CHR19_NB_HC_out_0.05_fish.txt.bed")

########################################################
# Step 4: extract the rf part based on chromosomes
########################################################

# using the script called /Users/HG/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/fisher_outlier/format/split_chr.sh