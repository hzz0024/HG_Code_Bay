######################
### outlier format ###
######################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/24_deeplearning")

format <- function(input){
  #input = "18_HC_NB_shared_outlier.txt"
  dat = read.delim(input, header = TRUE, sep='\t')
  angsd_list <- data.frame(dat$chr , dat$pos,  dat$pos)
  return(angsd_list)
}

d1 <- format("18_HC_NB_shared_outlier.txt")
d2 <- format("19_HC_NB_shared_outlier.txt")
d3 <- format("21_HC_NB_shared_outlier.txt")
d4 <- format("19_Sur_Ref_shared_outlier.txt")
d5 <- format("20_Sur_Ref_shared_outlier.txt")

df <- rbind(d1,d2,d3)
df1 <- df[!duplicated(df), ]
df2 = df1[with(df1, order(dat.chr, dat.pos)),]
df_final <- paste0(df2$dat.chr, ":" , df2$dat.pos, "-", df2$dat.pos)
write.table(df_final, "DB.wild.outlier.6329.rf.txt", row.names=F, col.names = F, quote=F, sep="\t")


