#######################
#  Adjust the p-value #
#######################
setwd("/Volumes/cornell/Fisher_exact")
#------------------------------------------------------------
# result by domajorminor 3 and domaf 1
library(qvalue)
write_table <- function(res_name, output_fdr1, output_fdr5){
  #res_name = deparse(substitute(dat))
  res <- read.delim(res_name, header = FALSE, sep='\t')
  p1 <- res$V3
  p0 <- res$V4
  p_value <- res$V7
  length(which(p_value<0.05))
  hist(p_value)
  qobj <- qvalue(p = p_value)
  fdr <- qobj$qvalues
  length(fdr[fdr<0.1])
  outlier <- as.data.frame(cbind(res$V1[fdr<0.1], res$V2[fdr<0.1], res$V4[fdr<0.1], res$V5[fdr<0.1], p_value[fdr<0.1], fdr[fdr<0.1]))
  colnames(outlier) <- c("chr","pos", "p0", "delta_p", "p_value","fdr")
  # export the dataframe to txt file
  write.table(outlier,output_fdr1, quote = FALSE, row.names = FALSE)
  length(fdr[fdr<0.05])
  outlier <- as.data.frame(cbind(res$V1[fdr<0.05], res$V2[fdr<0.05], res$V4[fdr<0.05], res$V5[fdr<0.05], p_value[fdr<0.05], fdr[fdr<0.05]))
  colnames(outlier) <- c("chr","pos", "p0", "delta_p", "p_value","fdr")
  # export the dataframe to txt file
  write.table(outlier,output_fdr5, quote = FALSE, row.names = FALSE)
  print('Done...')
  return(length(fdr[fdr<0.05]))
}

#------------------------------------------------------------
# result by domajorminor 3 and domaf 1
#write_table("fish_CH_REF_HCCOH.txt",  "./results/Fish_CH_REF_HCCOH_fdr1.txt" , "./results/Fish_CH_REF_HCCOH_fdr5.txt")
#write_table("fish_CH_REF_HCNB.txt",  "./results/Fish_CH_REF_HCNB_fdr1.txt" , "./results/Fish_CH_REF_HCNB_fdr5.txt")
#write_table("fish_HC_COH.txt",  "./results/Fish_HC_COH_fdr1.txt" , "./results/Fish_HC_COH_fdr5.txt")
#write_table("fish_HC_NB.txt", "./results/Fish_HC_NB_fdr1.txt" , "./results/Fish_HC_NB_fdr5.txt")

# result by domajorminor 5 and domaf 2
#write_table("fish_CH_REF_HCCOH_anc.txt",  "./results/Fish_CH_REF_HCCOH_anc_fdr1.txt" , "./results/Fish_CH_REF_HCCOH_anc_fdr5.txt")
#write_table("fish_CH_REF_HCNB_anc.txt", "./results/Fish_CH_REF_HCNB_anc_fdr1.txt" , "./results/Fish_CH_REF_HCNB_anc_fdr5.txt")
#write_table("fish_HC_COH_anc.txt",  "./results/Fish_HC_COH_anc_fdr1.txt" , "./results/Fish_HC_COH_anc_fdr5.txt")
#write_table("fish_HC_NB_anc.txt", "./results/Fish_HC_NB_anc_fdr1.txt" , "./results/Fish_HC_NB_anc_fdr5.txt")

#######################
#  check shared SNPs  #
#######################
setwd("~/Documents/Ryan_workplace/DelBay/Selection_Fish/Fish_repeat_results")

make_id <- function(dt_name){
  dat <- read.delim(dt_name, header = FALSE, sep='\t')
  dat_ID = paste0(dat$V1,'_',dat$V2)
  return(dat_ID)
}

file1 = "REF-CH-NB-HC_out_0.1.txt"
file2 = "REF-CH-SR-HC_out_0.1.txt"

a <- make_id(file1)
b <- make_id(file2)
c = intersect(a, b)

message('First file \"',file1,'\": ',length(a),', second file \"',file2, '\": ',length(b),'. Intersection: ',length(c))


#################################
# Check major and minor alleles #
#################################

file1 = 'CH_maf0.05_pctind0.7_cv30_anc.mafs'
file2 = 'CH_maf0.05_pctind0.7_cv30.mafs'
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat2 = read.delim(file2, header = TRUE, sep='\t')
names1 = paste0(dat1$chromo,'_',dat1$position)
names2 = paste0(dat2$chromo,'_',dat2$position)
datt1 = dat1[names1 %in% c,]
datt2 = dat2[names2 %in% c,]
datt1 = dat1
datt2 = dat2

index1 =  datt1$major == datt2$minor & datt1$minor == datt2$major
index = index1
diff_list1 = paste(datt1$chromo[index], datt1$position[index], datt1$major[index], datt1$minor[index], datt1$anc[index], datt1$unknownEM[index], datt2$major[index], datt2$minor[index], datt2$anc[index], datt2$knownEM[index], sep = ' ')
length(diff_list1)

index2 =  datt1$major == datt2$major & datt1$minor == datt2$minor
index = index2
diff_list2 = paste(datt1$chromo[index], datt1$position[index], datt1$major[index], datt1$minor[index], datt1$anc[index], datt1$knownEM[index], datt2$major[index], datt2$minor[index], datt2$anc[index],datt2$knownEM[index], sep = ' ')
length(diff_list2)

index =   !(index1 | index2)
diff_list = paste(datt1$chromo[index], datt1$position[index], datt1$major[index], datt1$minor[index], datt1$anc[index], datt1$unknownEM[index], datt2$major[index], datt2$minor[index], datt2$anc[index], datt2$knownEM[index], sep = ' ')
length(diff_list)

length(diff_list1)+length(diff_list2)+length(diff_list)
write.table(diff_list, 'diff.txt', quote=FALSE, sep = "\t", row.names=FALSE, col.names=FALSE)

