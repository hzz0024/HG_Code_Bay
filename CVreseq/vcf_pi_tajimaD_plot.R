# load the R function
setwd("~/Documents/Ryan_workplace/CVreseq_evaluation/plot")
source("manhattan.R")

###################### start plotting pi ######################
plot_pi <- function(name1, name2, plot_name, pop){
  # load the original windowed pi
  # name1 = "SL_chr2_original_5000.windowed.pi"
  DT1 = read.delim(name1, header = TRUE, sep='\t')
  id1 = paste0(DT1$CHROM,'_',DT1$BIN_START)
  DT1 <- as.data.frame(cbind(DT1,id1))
  # load the mask windowed pi
  #name2 = "SL_chr2_mask_5000.windowed.pi"
  DT2 = read.delim(name2, header = TRUE, sep='\t')
  id2 = paste0(DT2$CHROM,'_',DT2$BIN_START)
  DT2 <- as.data.frame(cbind(DT2,id2))
  common_id = intersect(id1,id2)
  idx1 = DT1$id1 %in% common_id
  idx2 = DT2$id2 %in% common_id
  diff = DT2$PI[idx2]-DT1$PI[idx1]
  dat <- as.data.frame(cbind(DT1$CHROM[idx1], DT1$BIN_START[idx1], common_id, DT1$PI[idx1], DT2$PI[idx2], diff))
  colnames(dat) <- c("chr","Bin_start", "SNP", "pi1", "pi2","diff")
  dat$chr <- as.numeric(dat$chr)
  dat$Bin_start <- as.numeric(dat$Bin_start)
  dat$pi1 <- as.numeric(dat$pi1)
  dat$pi2 <- as.numeric(dat$pi2)
  dat$diff <- as.numeric(dat$diff)
  #return(dat)
  print('Prepare Data Done')
  ###################### start pi plot ######################
  jpeg(plot_name, width = 16, height = 9, units = 'in', res = 300)
  par(mfrow=c(3,1))
  manhattan(chr="chr",bp="Bin_start",p="pi1", subset(dat, chr == 2), logp=FALSE, cex.axis = 1, ylim = c(0, 0.005),
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab=expression(paste(pi)), cex.lab=2, main = paste(pop, "_original (5000 bp/window)"), cex.main=2)
  manhattan(chr="chr",bp="Bin_start",p="pi2", subset(dat, chr == 2), logp=FALSE, cex.axis = 1, ylim = c(0, 0.005),
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab=expression(paste(pi)), cex.lab=2, main = paste(pop, "_mask (5000 bp/window)"), cex.main=2)
  manhattan(chr="chr",bp="Bin_start",p="diff", subset(dat, chr == 2), logp=FALSE, cex.axis = 1, ylim = c(-5e-4, 5e-4),
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab=expression(paste(pi)), cex.lab=2, main = "Mask minus Original (5000 bp/window)", cex.main=2)
  dev.off()
  print("Pi Plotting done")
}


plot_pi2 <- function(name1, name2, plot_name, pop1, pop2){
  #  load the original windowed pi
  #name1 = "SL_chr2_original_5000.windowed.pi"
  DT1 = read.delim(name1, header = TRUE, sep='\t')
  id1 = paste0(DT1$CHROM,'_',DT1$BIN_START)
  DT1 <- as.data.frame(cbind(DT1,id1))
  # load the mask windowed pi
  #name2 = "SL_chr2_mask_5000.windowed.pi"
  DT2 = read.delim(name2, header = TRUE, sep='\t')
  id2 = paste0(DT2$CHROM,'_',DT2$BIN_START)
  DT2 <- as.data.frame(cbind(DT2,id2))
  common_id = intersect(id1,id2)
  idx1 = DT1$id1 %in% common_id
  idx2 = DT2$id2 %in% common_id
  diff = DT2$PI[idx2]-DT1$PI[idx1]
  dat <- as.data.frame(cbind(DT1$CHROM[idx1], DT1$BIN_START[idx1], common_id, DT1$PI[idx1], DT2$PI[idx2], diff))
  colnames(dat) <- c("chr","Bin_start", "SNP", "pi1", "pi2","diff")
  dat$chr <- as.numeric(dat$chr)
  dat$Bin_start <- as.numeric(dat$Bin_start)
  dat$pi1 <- as.numeric(dat$pi1)
  dat$pi2 <- as.numeric(dat$pi2)
  dat$diff <- as.numeric(dat$diff)
  #return(dat)
  print('Prepare Data Done')
  ###################### start pi plot ######################
  jpeg(plot_name, width = 16, height = 9, units = 'in', res = 300)
  par(mfrow=c(1,1))
  #manhattan(chr="chr",bp="Bin_start",p="pi1", subset(dat, chr == 2), logp=FALSE, cex.axis = 1, ylim = c(0, 0.005),
  #          col=c("grey","black"),genomewideline=F, suggestiveline=F,
  #          ylab=expression(paste(pi)), cex.lab=2, main = paste(pop1, "_original (5000 bp/window)"), cex.main=2)
  #manhattan(chr="chr",bp="Bin_start",p="pi2", subset(dat, chr == 2), logp=FALSE, cex.axis = 1, ylim = c(0, 0.005),
  #          col=c("grey","black"),genomewideline=F, suggestiveline=F,
  #          ylab=expression(paste(pi)), cex.lab=2, main = paste(pop2, "_original (5000 bp/window)"), cex.main=2)
  manhattan(chr="chr",bp="Bin_start",p="diff", subset(dat, chr == 2), logp=FALSE, cex.axis = 1, ylim = c(-8e-4, 8e-4),
            col=c("grey","black"),genomewideline=T, suggestiveline=F,
            ylab=expression(paste(pi)), cex.lab=2, main = "SL-Louisiana wild line vs. OBOYS2-selected line (5000 bp/window)", cex.main=2)
  abline(h=0,col='red')
  dev.off()
  print("Pi Plotting done")
}
SL_chr2 <- plot_pi("SL_chr2_original_5000.windowed.pi","SL_chr2_mask_5000.windowed.pi", "SL_pi_chr2_w5000.jpg", "SL" )
OBOYS2_chr2 <- plot_pi("OBOYS2_chr2_original_5000.windowed.pi","OBOYS2_chr2_mask_5000.windowed.pi", "OBOYS2_pi_chr2_w5000.jpg", "OBOYS2" )
OBOYS2_chr2 <- plot_pi2("OBOYS2_chr2_original_5000.windowed.pi","SL_chr2_original_5000.windowed.pi", "SL_OBOYS2_pi_original_chr2_w5000.jpg", "OBOYS2", "SL")
OBOYS2_chr2 <- plot_pi2("OBOYS2_chr2_mask_5000.windowed.pi","SL_chr2_mask_5000.windowed.pi", "SL_OBOYS2_pi_mask_chr2_w5000.jpg", "OBOYS2", "SL")

###################### start plotting Tajima's D ######################
plot_TajimaD <- function(name1, name2, plot_name, pop){
  #name1 = "OBOYS2_chr2_original_D5000.Tajima.D"
  DT1 = read.delim(name1, header = TRUE, sep='\t')
  DT1 <- DT1[complete.cases(DT1), ]
  id1 = paste0(DT1$CHROM,'_',DT1$BIN_START)
  DT1 <- as.data.frame(cbind(DT1,id1))
  #name2 = "OBOYS2_chr2_mask_D5000.Tajima.D"
  DT2 = read.delim(name2, header = TRUE, sep='\t')
  DT2 <- DT2[complete.cases(DT2), ]
  id2 = paste0(DT2$CHROM,'_',DT2$BIN_START)
  DT2 <- as.data.frame(cbind(DT2,id2))
  common_id = intersect(id1,id2)
  idx1 = DT1$id1 %in% common_id
  idx2 = DT2$id2 %in% common_id
  diff = DT2$TajimaD[idx2]-DT1$TajimaD[idx1]
  dat <- as.data.frame(cbind(DT1$CHROM[idx1], DT1$BIN_START[idx1], common_id, DT1$TajimaD[idx1], DT2$TajimaD[idx2], diff))
  colnames(dat) <- c("chr","Bin_start", "SNP", "TD1", "TD2","diff")
  dat$chr <- as.numeric(dat$chr)
  dat$Bin_start <- as.numeric(dat$Bin_start)
  dat$TD1 <- as.numeric(dat$TD1)
  dat$TD2 <- as.numeric(dat$TD2)
  dat$diff <- as.numeric(dat$diff)
  print('Prepare Data Done')
  jpeg(plot_name, width = 16, height = 9, units = 'in', res = 300)
  par(mfrow=c(3,1))
  manhattan(chr="chr",bp="Bin_start",p="TD1", subset(dat, chr == 2), logp=FALSE, cex.axis = 1, ylim = c(-3, 3),
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="Tajima's D", cex.lab=2, main = paste(pop, "_original (5000 bp/window)"), cex.main=2)
  manhattan(chr="chr",bp="Bin_start",p="TD2", subset(dat, chr == 2), logp=FALSE, cex.axis = 1, ylim = c(-3, 3),
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="Tajima's D", cex.lab=2, main = paste(pop, "_mask (5000 bp/window)"), cex.main=2)
  manhattan(chr="chr",bp="Bin_start",p="diff", subset(dat, chr == 2), logp=FALSE, cex.axis = 1, ylim = c(-3, 3),
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="Tajima's D", cex.lab=2, main = "Mask minus Original (5000 bp/window)", cex.main=2)
  dev.off()
  print("TajimaD Plotting done")
}

SL_chr2 <- plot_TajimaD("SL_chr2_original_D5000.Tajima.D","SL_chr2_mask_D5000.Tajima.D", "SL_TajimaD_chr2_w5000.jpg", "SL" )
OBOYS2_chr2 <- plot_TajimaD("OBOYS2_chr2_original_D5000.Tajima.D","OBOYS2_chr2_mask_D5000.Tajima.D", "OBOYS2_TajimaD_chr2_w5000.jpg", "OBOYS2" )

