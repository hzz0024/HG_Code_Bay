setwd("/Volumes/cornell/CVreseq_fst")
source("manhattan.R")
library(export)

####################### loop over different windows ####################
####################### loop over different windows ####################
####################### loop over different windows ####################
jpeg("CS_DEBY_fst_chr8_vcftools.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
name1 = "CS_DEBY_single_SNP.weir.fst"
DT1 = read.delim(name1, header = TRUE, sep='\t')
id1 = paste0(DT1$CHROM,'_',DT1$POS)
DT1 <- as.data.frame(cbind(DT1,id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]
dat1 <- data.frame(chr=DT1$CHROM, pos=DT1$POS, SNP=DT1$id1, fst=DT1$WEIR_AND_COCKERHAM_FST)
dat1$chr <- as.numeric(dat1$chr)
dat1$pos <- as.numeric(dat1$pos)
dat1$fst <- as.numeric(dat1$fst)
manhattan(chr="chr",bp="pos",p="fst", subset(dat1, chr == 8), logp=FALSE, cex.axis = 1, ylim = c(min(dat1$fst), 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = "CS vs DEBY single SNP Fst", cex.main=1.5)
for(win in c(50, 100, 500, 1000, 5000)){
  #name = paste0("CS_DEBY.chr8.", win, "bp.csv") # hudson
  #DT = read.delim(name, header = TRUE, sep=',') # hudson
  #mid_pos <- DT$mid # hudson
  #id = paste0(DT$scaffold,'_',mid_pos) # hudson
  #DT <- as.data.frame(cbind(DT,mid_pos, id)) # hudson
  # remove NA rows in the dxy column
  name = paste0("CS_DEBY_", win, "bp.windowed.weir.fst")
  DT = read.delim(name, header = TRUE, sep='\t')
  mid_pos <- round((DT$BIN_START + DT$BIN_END)/2)
  id = paste0(DT$CHROM,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(chr=DT$CHROM, start=DT$BIN_START, end=DT$BIN_END, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT$WEIGHTED_FST)
  dat$chr <- as.numeric(dat$chr)
  dat$mid_pos <- as.numeric(dat$mid_pos)
  dat$fst <- as.numeric(dat$fst)
  #dat <- data.frame(chr=DT$scaffold, mid_pos=DT$mid, SNP=DT$id, fst=DT$Fst_CS_DEBY) # hudson
  manhattan(chr="chr",bp="mid_pos",p="fst", subset(dat, chr == 8), logp=FALSE, cex.axis = 1, ylim = c(0, 1), 
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = paste0("CS vs DEBY Fst window size:", win, ' bp'), cex.main=1.5)
}
dev.off()

jpeg("CS_NEH_fst_chr8_vcftools.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
name1 = "CS_NEH_single_SNP.weir.fst"
DT1 = read.delim(name1, header = TRUE, sep='\t')
id1 = paste0(DT1$CHROM,'_',DT1$POS)
DT1 <- as.data.frame(cbind(DT1,id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]
dat1 <- data.frame(chr=DT1$CHROM, pos=DT1$POS, SNP=DT1$id1, fst=DT1$WEIR_AND_COCKERHAM_FST)
dat1$chr <- as.numeric(dat1$chr)
dat1$pos <- as.numeric(dat1$pos)
dat1$fst <- as.numeric(dat1$fst)
manhattan(chr="chr",bp="pos",p="fst", subset(dat1, chr == 8), logp=FALSE, cex.axis = 1, ylim = c(min(dat1$fst), 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = "CS vs NEH single SNP Fst", cex.main=1.5)
for(win in c(50, 100, 500, 1000, 5000)){
  #name = paste0("CS_NEH.chr8.", win, "bp.csv") # hudson
  #DT = read.delim(name, header = TRUE, sep=',') # hudson
  #mid_pos <- DT$mid # hudson
  #id = paste0(DT$scaffold,'_',mid_pos) # hudson
  #DT <- as.data.frame(cbind(DT,mid_pos, id)) # hudson
  # remove NA rows in the dxy column
  name = paste0("CS_NEH_", win, "bp.windowed.weir.fst")
  DT = read.delim(name, header = TRUE, sep='\t')
  mid_pos <- round((DT$BIN_START + DT$BIN_END)/2)
  id = paste0(DT$CHROM,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(chr=DT$CHROM, start=DT$BIN_START, end=DT$BIN_END, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT$WEIGHTED_FST)
  dat$chr <- as.numeric(dat$chr)
  dat$mid_pos <- as.numeric(dat$mid_pos)
  dat$fst <- as.numeric(dat$fst)
  #dat <- data.frame(chr=DT$scaffold, mid_pos=DT$mid, SNP=DT$id, fst=DT$Fst_CS_NEH) # hudson
  manhattan(chr="chr",bp="mid_pos",p="fst", subset(dat, chr == 8), logp=FALSE, cex.axis = 1, ylim = c(0, 1), 
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="Fst", cex.lab=1.5, main = paste0("CS vs NEH Fst window size:", win, ' bp'), cex.main=1.5)
}
dev.off()

jpeg("SL_OBOYS2_fst_chr8_vcftools.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
name1 = "SL_OBOYS2_single_SNP.weir.fst"
DT1 = read.delim(name1, header = TRUE, sep='\t')
id1 = paste0(DT1$CHROM,'_',DT1$POS)
DT1 <- as.data.frame(cbind(DT1,id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]
dat1 <- data.frame(chr=DT1$CHROM, pos=DT1$POS, SNP=DT1$id1, fst=DT1$WEIR_AND_COCKERHAM_FST)
dat1$chr <- as.numeric(dat1$chr)
dat1$pos <- as.numeric(dat1$pos)
dat1$fst <- as.numeric(dat1$fst)
manhattan(chr="chr",bp="pos",p="fst", subset(dat1, chr == 8), logp=FALSE, cex.axis = 1, ylim = c(min(dat1$fst), 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = "SL vs OBOYS2 single SNP Fst", cex.main=1.5)
for(win in c(50, 100, 500, 1000, 5000)){
  #name = paste0("SL_OBOYS2.chr8.", win, "bp.csv") # hudson
  #DT = read.delim(name, header = TRUE, sep=',') # hudson
  #mid_pos <- DT$mid # hudson
  #id = paste0(DT$scaffold,'_',mid_pos) # hudson
  #DT <- as.data.frame(cbind(DT,mid_pos, id)) # hudson
  # remove NA rows in the dxy column
  name = paste0("SL_OBOYS2_", win, "bp.windowed.weir.fst")
  DT = read.delim(name, header = TRUE, sep='\t')
  mid_pos <- round((DT$BIN_START + DT$BIN_END)/2)
  id = paste0(DT$CHROM,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(chr=DT$CHROM, start=DT$BIN_START, end=DT$BIN_END, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT$WEIGHTED_FST)
  dat$chr <- as.numeric(dat$chr)
  dat$mid_pos <- as.numeric(dat$mid_pos)
  dat$fst <- as.numeric(dat$fst)
  #dat <- data.frame(chr=DT$scaffold, mid_pos=DT$mid, SNP=DT$id, fst=DT$Fst_SL_OBOYS2) # hudson
  manhattan(chr="chr",bp="mid_pos",p="fst", subset(dat, chr == 8), logp=FALSE, cex.axis = 1, ylim = c(0, 1), 
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="Fst", cex.lab=1.5, main = paste0("SL vs OBOYS2 Fst window size:", win, ' bp'), cex.main=1.5)
}
dev.off()

############################# from SM output ############################# 
jpeg("CS_DEBY_fst_chr8.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
for(win in c(50, 100, 500, 1000, 5000, 10000)){
  name = paste0("CS_DEBY.chr8.", win, "bp.csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- DT$mid
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  # remove NA rows in the dxy column
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(chr=DT$scaffold, mid_pos=DT$mid, SNP=DT$id, fst=DT$Fst_CS_DEBY)
  manhattan(chr="chr",bp="mid_pos",p="fst", subset(dat, chr == 8), logp=FALSE, cex.axis = 1, ylim = c(0, 1), 
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="Fst", cex.lab=1.5, main = paste0("CS vs DEBY Fst window size:", win, ' bp'), cex.main=1.5)
}
dev.off()

jpeg("CS_NEH_fst_chr8.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
for(win in c(50, 100, 500, 1000, 5000, 10000)){
  name = paste0("CS_NEH.chr8.", win, "bp.csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- DT$mid
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  # remove NA rows in the dxy column
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(chr=DT$scaffold, mid_pos=DT$mid, SNP=DT$id, fst=DT$Fst_CS_NEH)
  manhattan(chr="chr",bp="mid_pos",p="fst", subset(dat, chr == 8), logp=FALSE, cex.axis = 1, ylim = c(0, 1), 
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="Fst", cex.lab=1.5, main = paste0("CS vs NEH Fst window size:", win, ' bp'), cex.main=1.5)
}
dev.off()

jpeg("SLOBOYS2_fst_chr8.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
for(win in c(50, 100, 500, 1000, 5000, 10000)){
  name = paste0("SL_OBOYS2.chr8.", win, "bp.csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- DT$mid
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  # remove NA rows in the dxy column
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(chr=DT$scaffold, mid_pos=DT$mid, SNP=DT$id, fst=DT$Fst_SL_OBOYS2)
  manhattan(chr="chr",bp="mid_pos",p="fst", subset(dat, chr == 8), logp=FALSE, cex.axis = 1, ylim = c(0, 1), 
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="Fst", cex.lab=1.5, main = paste0("SL vs OBOYS2 Fst window size:", win, ' bp'), cex.main=1.5)
}
dev.off()


####################### windows plot ####################
####################### windows plot ####################
####################### windows plot ####################

name1 = "CS_DEBY_chr8.csv"
DT1 = read.delim(name1, header = TRUE, sep=',')
mid_pos <- DT1$mid
id1 = paste0(DT1$scaffold,'_',mid_pos)
DT1 <- as.data.frame(cbind(DT1,mid_pos, id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]
dat1 <- data.frame(chr=DT1$scaffold, mid_pos=DT1$mid, SNP=DT1$id1, fst=DT1$Fst_CS_DEBY)

name2 = "CS_NEH_chr8.csv"
DT2 = read.delim(name2, header = TRUE, sep=',')
mid_pos <-  DT2$mid
id2 = paste0(DT2$scaffold,'_',mid_pos)
DT2 <- as.data.frame(cbind(DT2,mid_pos,id2))
# remove NA rows in the dxy column
DT2 <- DT2[complete.cases(DT2), ]
dat2 <- data.frame(chr=DT2$scaffold, mid_pos=DT2$mid, SNP=DT2$id2, fst=DT2$Fst_CS_NEH)

name3 = "SL_OBOYS2_chr8.csv"
DT3 = read.delim(name3, header = TRUE, sep=',')
mid_pos <-  DT3$mid
id3 = paste0(DT3$scaffold,'_',mid_pos)
DT3 <- as.data.frame(cbind(DT3,mid_pos,id3))
# remove NA rows in the dxy column
DT3 <- DT3[complete.cases(DT3), ]
dat3 <- data.frame(chr=DT3$scaffold, mid_pos=DT3$mid, SNP=DT3$id3, fst=DT3$Fst_SL_OBOYS2)

file1 = 'snpsOfInterest_chr8.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)

jpeg("fst_chr8.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chr",bp="mid_pos",p="fst", subset(dat1, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = "CS vs DEBY Fst (1000 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="fst", subset(dat2, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = "CS vs NEH Fst (1000 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="fst", subset(dat3, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = "SL vs OBOYS2 Fst (1000 bp/window)", cex.main=1.5)
dev.off()
#  print("Pi Plotting done")

##################################### vcftools single SNP fst  plot ##################################
##################################### vcftools single SNP fst  plot ##################################
##################################### vcftools single SNP fst  plot ##################################

# import the highlight snps list
file1 = 'snpsOfInterest_chr8.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
print('Prepare Data Done')

name1 = "CS_DEBY_chr8_fst.weir.fst"
DT1 = read.delim(name1, header = TRUE, sep='\t')
id1 = paste0(DT1$CHROM,'_',DT1$POS)
DT1 <- as.data.frame(cbind(DT1,id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]

dat1 <- data.frame(chr=DT1$CHROM, pos=DT1$POS, SNP=DT1$id1, fst=DT1$WEIR_AND_COCKERHAM_FST)
dat1$chr <- as.numeric(dat1$chr)
dat1$pos <- as.numeric(dat1$pos)
dat1$fst <- as.numeric(dat1$fst)

h1$V1 %in% dat1$SNP
dat1[dat1$SNP %in% h1$V1,]

name2 = "CS_NEH_chr8_fst.weir.fst"
DT2 = read.delim(name2, header = TRUE, sep='\t')
id2 = paste0(DT2$CHROM,'_',DT2$POS)
DT2 <- as.data.frame(cbind(DT2,id2))
# remove NA rows in the dxy column
DT2 <- DT2[complete.cases(DT2), ]

dat2 <- data.frame(chr=DT2$CHROM, pos=DT2$POS, SNP=DT2$id2, fst=DT2$WEIR_AND_COCKERHAM_FST)
dat2$chr <- as.numeric(dat2$chr)
dat2$pos <- as.numeric(dat2$pos)
dat2$fst <- as.numeric(dat2$fst)

h1$V1 %in% dat2$SNP
dat2[dat2$SNP %in% h1$V1,]

name3 = "SL_OBOYS2_chr8_fst.weir.fst"
DT3 = read.delim(name3, header = TRUE, sep='\t')
id3 = paste0(DT3$CHROM,'_',DT3$POS)
DT3 <- as.data.frame(cbind(DT3,id3))
# remove NA rows in the dxy column
DT3 <- DT3[complete.cases(DT3), ]

dat3 <- data.frame(chr=DT3$CHROM, pos=DT3$POS, SNP=DT3$id3, fst=DT3$WEIR_AND_COCKERHAM_FST)
dat3$chr <- as.numeric(dat3$chr)
dat3$pos <- as.numeric(dat3$pos)
dat3$fst <- as.numeric(dat3$fst)

h1$V1 %in% dat3$SNP
dat3[dat3$SNP %in% h1$V1,]

jpeg("Pairwise_weir_fst_single_snp.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chr",bp="pos",p="fst", subset(dat1, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(min(dat1$fst), 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = "CS DEBY single SNP Fst", cex.main=1.5)
manhattan(chr="chr",bp="pos",p="fst", subset(dat2, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(-0.2, 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = "CS NEH single SNP Fst", cex.main=1.5)
manhattan(chr="chr",bp="pos",p="fst", subset(dat3, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(-0.2, 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.5, main = "SL OBOYS2 single SNP Fst", cex.main=1.5)
#graph2ppt(file="pi_CS_DEBY_500K.pptx", width=9, height=6.5)
dev.off()
print("Pi Plotting done")
