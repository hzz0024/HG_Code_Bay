setwd("/Volumes/cornell/CVreseq_density")
library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)

# for i in 100 500 1000 5000; do
# vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.CS.recode.vcf --SNPdensity $i --out 'CS.'$i'bp'
# vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.DEBY.recode.vcf --SNPdensity $i --out 'DEBY.'$i'bp'
# vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.NEH.recode.vcf --SNPdensity $i --out 'NEH.'$i'bp'
# vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.OBOYS2.recode.vcf --SNPdensity $i --out 'OBOYS2.'$i'bp'
# vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.CS.recode.vcf --SNPdensity $i --out 'SL.'$i'bp'
# done
par(mfrow=c(2,2))
jpeg("SNP_density_maf01.jpg", width = 16, height = 9, units = 'in', res = 300)
for(win in c(500, 1000, 5000, 10000)){
# load the original density results
  name1 = toString(paste0("CS.maf01.",win,"bp.snpden"))
  DT1 = read.delim(name1, header = TRUE, sep='\t')
  dat1 <- data.frame(chr=DT1$CHROM, pos=DT1$BIN_START, count=DT1$SNP_COUNT/win)
  name2 = toString(paste0("DEBY.maf01.",win,"bp.snpden"))
  DT2 = read.delim(name2, header = TRUE, sep='\t')
  dat2 <- data.frame(chr=DT2$CHROM, pos=DT2$BIN_START, count=DT2$SNP_COUNT/win)
  name3 = toString(paste0("NEH.maf01.",win,"bp.snpden"))
  DT3 = read.delim(name3, header = TRUE, sep='\t')
  dat3 <- data.frame(chr=DT3$CHROM, pos=DT3$BIN_START, count=DT3$SNP_COUNT/win)
  name4 = toString(paste0("SL.maf01.",win,"bp.snpden"))
  DT4 = read.delim(name4, header = TRUE, sep='\t')
  dat4 <- data.frame(chr=DT4$CHROM, pos=DT4$BIN_START, count=DT4$SNP_COUNT/win)
  name5 = toString(paste0("OBOYS2.maf01.",win,"bp.snpden"))
  DT5 = read.delim(name5, header = TRUE, sep='\t')
  dat5 <- data.frame(chr=DT5$CHROM, pos=DT5$BIN_START, count=DT5$SNP_COUNT/win)
  plotdat = data.frame(SNP_density=c(dat1$count,dat2$count,dat3$count, dat4$count, dat5$count), 
                       Population=c(rep('CS',length(dat1$count)), rep('DEBY',length(dat2$count)), rep('NEH',length(dat3$count)), rep('SL',length(dat4$count)), rep('OBOYS2',length(dat5$count))))
  p <- ggplot(plotdat, aes(x=Population, y=SNP_density, fill=Population)) +
       geom_violin(trim=FALSE) + geom_boxplot(width = 0.05) + theme_classic() +
       labs(title=paste0("Frequency of segregation sites at window ", win, "bp, maf 0.01"),x="Population", y = "SNP density") 
  p +  scale_color_manual(values=c("#FFAEBC", "#A0E7E5", "#B4F8C8", "#FBE7C6", "#FFA384"))
} #graph2ppt(file="pi_CS_NEH_1K.pptx", width=9, height=6.5)
dev.off()

par(mfrow=c(2,2))
jpeg("SNP_density_maf05.jpg", width = 16, height = 9, units = 'in', res = 300)
for(win in c(500, 1000, 5000, 10000)){
  # load the original density results
  name1 = toString(paste0("CS.",win,"bp.snpden"))
  DT1 = read.delim(name1, header = TRUE, sep='\t')
  dat1 <- data.frame(chr=DT1$CHROM, pos=DT1$BIN_START, count=DT1$SNP_COUNT/win)
  name2 = toString(paste0("DEBY.",win,"bp.snpden"))
  DT2 = read.delim(name2, header = TRUE, sep='\t')
  dat2 <- data.frame(chr=DT2$CHROM, pos=DT2$BIN_START, count=DT2$SNP_COUNT/win)
  name3 = toString(paste0("NEH.",win,"bp.snpden"))
  DT3 = read.delim(name3, header = TRUE, sep='\t')
  dat3 <- data.frame(chr=DT3$CHROM, pos=DT3$BIN_START, count=DT3$SNP_COUNT/win)
  name4 = toString(paste0("SL.",win,"bp.snpden"))
  DT4 = read.delim(name4, header = TRUE, sep='\t')
  dat4 <- data.frame(chr=DT4$CHROM, pos=DT4$BIN_START, count=DT4$SNP_COUNT/win)
  name5 = toString(paste0("OBOYS2.",win,"bp.snpden"))
  DT5 = read.delim(name5, header = TRUE, sep='\t')
  dat5 <- data.frame(chr=DT5$CHROM, pos=DT5$BIN_START, count=DT5$SNP_COUNT/win)
  plotdat = data.frame(SNP_density=c(dat1$count,dat2$count,dat3$count, dat4$count, dat5$count), 
                       Population=c(rep('CS',length(dat1$count)), rep('DEBY',length(dat2$count)), rep('NEH',length(dat3$count)), rep('SL',length(dat4$count)), rep('OBOYS2',length(dat5$count))))
  p <- ggplot(plotdat, aes(x=Population, y=SNP_density, fill=Population)) +
    geom_violin(trim=FALSE) + geom_boxplot(width = 0.05) + theme_classic() +
    labs(title=paste0("Frequency of segregation sites at window ", win, "bp, maf 0.05"),x="Population", y = "SNP density") 
  p +  scale_color_manual(values=c("#FFAEBC", "#A0E7E5", "#B4F8C8", "#FBE7C6", "#FFA384"))
} #graph2ppt(file="pi_CS_NEH_1K.pptx", width=9, height=6.5)
dev.off()


