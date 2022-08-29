#!/usr/bin/env Rscript
# Honggang Zhao 03252022
# usage: Rscript --vanilla ngsparalog.R input.txt out.txt[default:snps_after_paralog_filtering_p05.txt]
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Input missing! The samtools pileup input file is required", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "snps_after_paralog_filtering_p05.txt"
}

#lr <- read.table("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.output")
lr <- read.table(args[1]) # read in ngsParalog calcLR output
lr$pval <- 0.5*pchisq(lr$V5,df=1,lower.tail=FALSE) # append column of p-values
lr$pval.adj <- p.adjust(lr$pval, method="BH") # p-values adjusted for number of tested sites
# The 7th column of the lr data.frame is the adjusted p-value for rejecting the null hypothesis that reads
# covering the site derive from a single locus. Of course you can use any p-value adjustment of your
# choosing, e.g. "fdr".

# generate list of sites that don't show evidence of mismapping at FDR 0.05 significance level:
qc.sites <- lr[-which(lr$pval.adj < 0.05),1:2]
# generate list of sites that don't show evidence of mismapping at 0.05 significance level:
#qc.sites <- lr[-which(lr$pval < 0.05),1:2]
write.table(qc.sites,file=args[2],quote=F,row.names=F,col.names=F)

