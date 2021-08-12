library(KRIS)
options(scipen=999)

prefix = "CS_NEH"
bed = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".noinvers.bed")
bim = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".noinvers.bim")
fam = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".noinvers.fam")
snp <- read.bed(bed, bim, fam )
sample_labels = c(rep('pop1',6), rep('pop2',6))
idx1 <- which(sample_labels == 'pop1')
idx2 <- which(sample_labels == 'pop2')
fst.pairwise <- fst.each.snp.hudson(snp$snp, idx1, idx2)
print(fst.pairwise[1:10])

dat.bedgraph <- data.frame(Chromosome=snp$snp.info$chr, Start=as.numeric(as.character(snp$snp.info$position))-1, End=snp$snp.info$position, Fst=fst.pairwise)
dat.bedgraph <- dat.bedgraph[complete.cases(dat.bedgraph), ]
dat.bedgraph$Chromosome <- as.numeric(as.character(dat.bedgraph$Chromosome))
dat.bedgraph$Start <- as.numeric(as.character(dat.bedgraph$Start))
dat.bedgraph$End <- as.numeric(as.character(dat.bedgraph$End))
dat.bedgraph$Fst <- as.numeric(as.character(dat.bedgraph$Fst))
dat.sort.bedgraph = dat.bedgraph[order(dat.bedgraph$Chromosome, dat.bedgraph$Start),]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(10)) 
  dat.sort.bedgraph$Chromosome[dat.sort.bedgraph$Chromosome==i] = chr_str_list[i]
write.table(dat.sort.bedgraph, file = paste0(prefix, ".sort.bedgraph"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)



prefix = "CS_DEBY"
bed = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".noinvers.bed")
bim = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".noinvers.bim")
fam = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".noinvers.fam")
snp <- read.bed(bed, bim, fam )
sample_labels = c(rep('pop1',6), rep('pop2',6))
idx1 <- which(sample_labels == 'pop1')
idx2 <- which(sample_labels == 'pop2')
fst.pairwise <- fst.each.snp.hudson(snp$snp, idx1, idx2)
print(fst.pairwise[1:10])

dat.bedgraph <- data.frame(Chromosome=snp$snp.info$chr, Start=as.numeric(as.character(snp$snp.info$position))-1, End=snp$snp.info$position, Fst=fst.pairwise)
dat.bedgraph <- dat.bedgraph[complete.cases(dat.bedgraph), ]
dat.bedgraph$Chromosome <- as.numeric(as.character(dat.bedgraph$Chromosome))
dat.bedgraph$Start <- as.numeric(as.character(dat.bedgraph$Start))
dat.bedgraph$End <- as.numeric(as.character(dat.bedgraph$End))
dat.bedgraph$Fst <- as.numeric(as.character(dat.bedgraph$Fst))
dat.sort.bedgraph = dat.bedgraph[order(dat.bedgraph$Chromosome, dat.bedgraph$Start),]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(10)) 
  dat.sort.bedgraph$Chromosome[dat.sort.bedgraph$Chromosome==i] = chr_str_list[i]
write.table(dat.sort.bedgraph, file = paste0(prefix, ".sort.bedgraph"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)



prefix = "SL_OBOYS2"
bed = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".noinvers.bed")
bim = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".noinvers.bim")
fam = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".noinvers.fam")
snp <- read.bed(bed, bim, fam )
sample_labels = c(rep('pop1',6), rep('pop2',6))
idx1 <- which(sample_labels == 'pop1')
idx2 <- which(sample_labels == 'pop2')
fst.pairwise <- fst.each.snp.hudson(snp$snp, idx1, idx2)
print(fst.pairwise[1:10])

dat.bedgraph <- data.frame(Chromosome=snp$snp.info$chr, Start=as.numeric(as.character(snp$snp.info$position))-1, End=snp$snp.info$position, Fst=fst.pairwise)
dat.bedgraph <- dat.bedgraph[complete.cases(dat.bedgraph), ]
dat.bedgraph$Chromosome <- as.numeric(as.character(dat.bedgraph$Chromosome))
dat.bedgraph$Start <- as.numeric(as.character(dat.bedgraph$Start))
dat.bedgraph$End <- as.numeric(as.character(dat.bedgraph$End))
dat.bedgraph$Fst <- as.numeric(as.character(dat.bedgraph$Fst))
dat.sort.bedgraph = dat.bedgraph[order(dat.bedgraph$Chromosome, dat.bedgraph$Start),]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(10)) 
  dat.sort.bedgraph$Chromosome[dat.sort.bedgraph$Chromosome==i] = chr_str_list[i]
write.table(dat.sort.bedgraph, file = paste0(prefix, ".sort.bedgraph"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
# /workdir/hz269/CVreseq_fst_SM_hudson/genome
# cat haplotig_masked_genome.size
# NC_035780.1	65668440
# NC_035781.1	61752955
# NC_035782.1	77061148
# NC_035783.1	59691872
# NC_035784.1	98698416
# NC_035785.1	51258098
# NC_035786.1	57830854
# NC_035787.1	75944018
# NC_035788.1	104168038
# NC_035789.1 32650045

# ./bedGraphToBigWig CS_NEH.sort.bedgraph haplotig_masked_genome.size CS_NEH.singleSNP.Hudson.fst.bw
# ./bedGraphToBigWig CS_DEBY.sort.bedgraph haplotig_masked_genome.size CS_DEBY.singleSNP.Hudson.fst.bw
# ./bedGraphToBigWig SL_OBOYS2.sort.bedgraph haplotig_masked_genome.size SL_OBOYS2.singleSNP.Hudson.fst.bw

prefix = "CS_HC"
bed = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.bed")
bim = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.bim")
fam = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.fam")
snp <- read.bed(bed, bim, fam )
sample_labels = c(rep('pop1',6), rep('pop2',6))
idx1 <- which(sample_labels == 'pop1')
idx2 <- which(sample_labels == 'pop2')
fst.pairwise <- fst.each.snp.hudson(snp$snp, idx1, idx2)
print(fst.pairwise[1:10])

dat.bedgraph <- data.frame(Chromosome=snp$snp.info$chr, Start=as.numeric(as.character(snp$snp.info$position))-1, End=snp$snp.info$position, Feature="Fst", Fst=fst.pairwise)
dat.bedgraph <- dat.bedgraph[complete.cases(dat.bedgraph), ]
dat.bedgraph$Chromosome <- as.numeric(as.character(dat.bedgraph$Chromosome))
dat.bedgraph$Start <- as.numeric(as.character(dat.bedgraph$Start))
dat.bedgraph$End <- as.numeric(as.character(dat.bedgraph$End))
dat.bedgraph$Fst <- as.numeric(as.character(dat.bedgraph$Fst))
dat.sort.bedgraph = dat.bedgraph[order(dat.bedgraph$Chromosome, dat.bedgraph$Start),]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(10)) 
  dat.sort.bedgraph$Chromosome[dat.sort.bedgraph$Chromosome==i] = chr_str_list[i]
write.table(dat.sort.bedgraph, file = paste0(prefix, ".sort.igv"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

prefix = "HCVA_CLP"
bed = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.bed")
bim = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.bim")
fam = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.fam")
snp <- read.bed(bed, bim, fam )
sample_labels = c(rep('pop1',6), rep('pop2',6))
idx1 <- which(sample_labels == 'pop1')
idx2 <- which(sample_labels == 'pop2')
fst.pairwise <- fst.each.snp.hudson(snp$snp, idx1, idx2)
print(fst.pairwise[1:10])

dat.bedgraph <- data.frame(Chromosome=snp$snp.info$chr, Start=as.numeric(as.character(snp$snp.info$position))-1, End=snp$snp.info$position, Fst=fst.pairwise)
dat.bedgraph <- dat.bedgraph[complete.cases(dat.bedgraph), ]
dat.bedgraph$Chromosome <- as.numeric(as.character(dat.bedgraph$Chromosome))
dat.bedgraph$Start <- as.numeric(as.character(dat.bedgraph$Start))
dat.bedgraph$End <- as.numeric(as.character(dat.bedgraph$End))
dat.bedgraph$Fst <- as.numeric(as.character(dat.bedgraph$Fst))
dat.sort.bedgraph = dat.bedgraph[order(dat.bedgraph$Chromosome, dat.bedgraph$Start),]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(10)) 
  dat.sort.bedgraph$Chromosome[dat.sort.bedgraph$Chromosome==i] = chr_str_list[i]
write.table(dat.sort.bedgraph, file = paste0(prefix, ".sort.bedgraph"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

prefix = "CS_HCVA"
bed = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.bed")
bim = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.bim")
fam = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.fam")
snp <- read.bed(bed, bim, fam )
sample_labels = c(rep('pop1',6), rep('pop2',6))
idx1 <- which(sample_labels == 'pop1')
idx2 <- which(sample_labels == 'pop2')
fst.pairwise <- fst.each.snp.hudson(snp$snp, idx1, idx2)
print(fst.pairwise[1:10])

dat.bedgraph <- data.frame(Chromosome=snp$snp.info$chr, Start=as.numeric(as.character(snp$snp.info$position))-1, End=snp$snp.info$position, Fst=fst.pairwise)
dat.bedgraph <- dat.bedgraph[complete.cases(dat.bedgraph), ]
dat.bedgraph$Chromosome <- as.numeric(as.character(dat.bedgraph$Chromosome))
dat.bedgraph$Start <- as.numeric(as.character(dat.bedgraph$Start))
dat.bedgraph$End <- as.numeric(as.character(dat.bedgraph$End))
dat.bedgraph$Fst <- as.numeric(as.character(dat.bedgraph$Fst))
dat.sort.bedgraph = dat.bedgraph[order(dat.bedgraph$Chromosome, dat.bedgraph$Start),]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(10)) 
  dat.sort.bedgraph$Chromosome[dat.sort.bedgraph$Chromosome==i] = chr_str_list[i]
write.table(dat.sort.bedgraph, file = paste0(prefix, ".sort.bedgraph"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

prefix = "HC_CLP"
bed = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.bed")
bim = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.bim")
fam = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.fam")
snp <- read.bed(bed, bim, fam )
sample_labels = c(rep('pop1',6), rep('pop2',6))
idx1 <- which(sample_labels == 'pop1')
idx2 <- which(sample_labels == 'pop2')
fst.pairwise <- fst.each.snp.hudson(snp$snp, idx1, idx2)
print(fst.pairwise[1:10])

dat.bedgraph <- data.frame(Chromosome=snp$snp.info$chr, Start=as.numeric(as.character(snp$snp.info$position))-1, End=snp$snp.info$position, Fst=fst.pairwise)
dat.bedgraph <- dat.bedgraph[complete.cases(dat.bedgraph), ]
dat.bedgraph$Chromosome <- as.numeric(as.character(dat.bedgraph$Chromosome))
dat.bedgraph$Start <- as.numeric(as.character(dat.bedgraph$Start))
dat.bedgraph$End <- as.numeric(as.character(dat.bedgraph$End))
dat.bedgraph$Fst <- as.numeric(as.character(dat.bedgraph$Fst))
dat.sort.bedgraph = dat.bedgraph[order(dat.bedgraph$Chromosome, dat.bedgraph$Start),]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(10)) 
  dat.sort.bedgraph$Chromosome[dat.sort.bedgraph$Chromosome==i] = chr_str_list[i]
write.table(dat.sort.bedgraph, file = paste0(prefix, ".sort.bedgraph"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

prefix = "CS_UMFS"
bed = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.bed")
bim = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.bim")
fam = paste0("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.sort.",prefix,".nobiginvers.fam")
snp <- read.bed(bed, bim, fam )
sample_labels = c(rep('pop1',6), rep('pop2',6))
idx1 <- which(sample_labels == 'pop1')
idx2 <- which(sample_labels == 'pop2')
fst.pairwise <- fst.each.snp.hudson(snp$snp, idx1, idx2)
print(fst.pairwise[1:10])
dat.bedgraph <- data.frame(Chromosome=snp$snp.info$chr, Start=as.numeric(as.character(snp$snp.info$position))-1, End=snp$snp.info$position, Fst=fst.pairwise)
dat.bedgraph <- dat.bedgraph[complete.cases(dat.bedgraph), ]
dat.bedgraph$Chromosome <- as.numeric(as.character(dat.bedgraph$Chromosome))
dat.bedgraph$Start <- as.numeric(as.character(dat.bedgraph$Start))
dat.bedgraph$End <- as.numeric(as.character(dat.bedgraph$End))
dat.bedgraph$Fst <- as.numeric(as.character(dat.bedgraph$Fst))
dat.sort.bedgraph = dat.bedgraph[order(dat.bedgraph$Chromosome, dat.bedgraph$Start),]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(10)) 
  dat.sort.bedgraph$Chromosome[dat.sort.bedgraph$Chromosome==i] = chr_str_list[i]
write.table(dat.sort.bedgraph, file = paste0(prefix, ".sort.bedgraph"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

##################################### vcftools fixed SNPs ##################################
##################################### vcftools fixed SNPs ##################################
##################################### vcftools fixed SNPs ##################################
name1 = "CS_NEH.sort.bedgraph"
DT1 = read.delim(name1, header = FALSE, sep='\t')
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
fixed <- DT1[DT1$V4 == max(DT1$V4),]
for(chr in chr_str_list){
  cnt = nrow(fixed[fixed$V1 == chr,])
  print(paste0("Chr ", chr, " has ", cnt, " fixed SNPs"))
}
dat.fixed <- data.frame(Chromosome=fixed$V1, Start=fixed$V2, End=fixed$V3, Feature = "fixed", Fst=fixed$V4)
write.table(dat.fixed, file = "CS_NEH_fixed.bed", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

name1 = "CS_DEBY.sort.bedgraph"
DT1 = read.delim(name1, header = FALSE, sep='\t')
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
fixed <- DT1[DT1$V4 == max(DT1$V4),]
for(chr in chr_str_list){
  cnt = nrow(fixed[fixed$V1 == chr,])
  print(paste0("Chr ", chr, " has ", cnt, " fixed SNPs"))
}
dat.fixed <- data.frame(Chromosome=fixed$V1, Start=fixed$V2, End=fixed$V3, Feature = "fixed", Fst=fixed$V4)
write.table(dat.fixed, file = "CS_DEBY_fixed.bed", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

name1 = "SL_OBOYS2.sort.bedgraph"
DT1 = read.delim(name1, header = FALSE, sep='\t')
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
fixed <- DT1[DT1$V4 == max(DT1$V4),]
for(chr in chr_str_list){
  cnt = nrow(fixed[fixed$V1 == chr,])
  print(paste0("Chr ", chr, " has ", cnt, " fixed SNPs"))
}
dat.fixed <- data.frame(Chromosome=fixed$V1, Start=fixed$V2, End=fixed$V3, Feature = "fixed", Fst=fixed$V4)
write.table(dat.fixed, file = "L_OBOYS2_fixed.bed", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

jpeg("Pairwise_weir_fst_fixed_SNP.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(dat1, chr="chr",bp="pos",p="fst", highlight1 = fixed$SNP, logp=FALSE, cex.axis = 1, ylim = c(min(dat1$fst), 1.02), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = "CS DEBY single SNP Fst", cex.main=1.5)
#graph2ppt(file="pi_CS_DEBY_500K.pptx", width=9, height=6.5)
dev.off()
print("Pi Plotting done")

