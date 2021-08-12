setwd("~/Documents/Ryan_workplace/CVreseq_pi")
source("manhattan.R")
library(export)

##################################### vcftools output plot ##################################
name1 = "CS_DEBY_chr8_CS_pi.sites.pi"
DT1 = read.delim(name1, header = TRUE, sep='\t')
id1 = paste0(DT1$CHROM,'_',DT1$POS)
DT1 <- as.data.frame(cbind(DT1,id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]

name2 = "CS_DEBY_chr8_DEBY_pi.sites.pi"
DT2 = read.delim(name2, header = TRUE, sep='\t')
id2 = paste0(DT2$CHROM,'_',DT2$POS)
DT2 <- as.data.frame(cbind(DT2,id2))
# remove NA rows in the dxy column
DT2 <- DT2[complete.cases(DT2), ]

name3 = "CS_NEH_chr8_CS_pi.sites.pi"
DT3 = read.delim(name3, header = TRUE, sep='\t')
id3 = paste0(DT3$CHROM,'_',DT3$POS)
DT3 <- as.data.frame(cbind(DT3,id3))
# remove NA rows in the dxy column
DT3 <- DT3[complete.cases(DT3), ]

name4 = "CS_NEH_chr8_NEH_pi.sites.pi"
DT4 = read.delim(name4, header = TRUE, sep='\t')
id4 = paste0(DT4$CHROM,'_',DT4$POS)
DT4 <- as.data.frame(cbind(DT4,id4))
# remove NA rows in the dxy column
DT4 <- DT4[complete.cases(DT4), ]

name5 = "SL_OBOYS2_chr8_SL_pi.sites.pi"
DT5 = read.delim(name5, header = TRUE, sep='\t')
id5 = paste0(DT5$CHROM,'_',DT5$POS)
DT5 <- as.data.frame(cbind(DT5,id5))
# remove NA rows in the dxy column
DT5 <- DT5[complete.cases(DT5), ]

name6 = "SL_OBOYS2_chr8_OBOYS2_pi.sites.pi"
DT6 = read.delim(name6, header = TRUE, sep='\t')
id6 = paste0(DT6$CHROM,'_',DT6$POS)
DT6 <- as.data.frame(cbind(DT6,id6))
# remove NA rows in the dxy column
DT6 <- DT6[complete.cases(DT6), ]

common_id = intersect(DT1$id1,DT2$id2)
common_id <- as.vector(common_id)
sub_DT1 = DT1[DT1$id1 %in% common_id,]
sub_DT2 = DT2[DT2$id2 %in% common_id,]

sub_DT1 = sub_DT1[match(common_id, sub_DT1$id1),]
sub_DT2 = sub_DT2[match(common_id, sub_DT2$id2),]

dat <- data.frame(chr=sub_DT1$CHROM, pos=sub_DT1$POS, SNP=common_id, D1_pi=sub_DT1$PI, D2_pi=sub_DT2$PI)
dat$chr <- as.numeric(dat$chr)
dat$pos <- as.numeric(dat$pos)
dat$D1_pi <- as.numeric(dat$D1_pi)
dat$D2_pi <- as.numeric(dat$D2_pi)
write.table(dat, file = "single_chr8_pi.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
# import the highlight snps list
file1 = 'snpsOfInterest_chr8.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
print('Prepare Data Done')
###################### start CS_DEBY plot ######################
#jpeg("test.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="pos",p="D1_pi", subset(dat, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "CS single SNP pi", cex.main=1.5)
manhattan(chr="chr",bp="pos",p="D2_pi", subset(dat, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 1), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "DEBY single SNP pi", cex.main=1.5)
graph2ppt(file="pi_CS_DEBY_500K.pptx", width=9, height=6.5)
#dev.off()
print("Pi Plotting done")

###################### start CS_NEH plot ######################
common_id = intersect(DT1$id1,DT3$id3)
common_id <- as.vector(common_id)
sub_DT1 = DT1[DT1$id1 %in% common_id,]
sub_DT3 = DT3[DT3$id3 %in% common_id,]

sub_DT1 = sub_DT1[match(common_id, sub_DT1$id1),]
sub_DT3 = sub_DT3[match(common_id, sub_DT3$id3),]

dat <- data.frame(chr=sub_DT1$CHROM, start=sub_DT1$BIN_START, start=sub_DT1$BIN_END, mid_pos=sub_DT1$mid_pos, SNP=common_id, D1_pi=sub_DT1$PI, D3_pi=sub_DT3$PI)
dat$chr <- as.numeric(dat$chr)
dat$mid_pos <- as.numeric(dat$mid_pos)
dat$D1_pi <- as.numeric(dat$D1_pi)
dat$D3_pi <- as.numeric(dat$D3_pi)
#write.table(dat, file = "500K_chr5_pi.csv", sep = ",", quote = FALSE,
#            row.names = FALSE, col.names = FALSE)
# import the highlight snps list
file1 = 'snpsOfInterest_chr3.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
file2 = 'snpsOfInterest_chr5.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)
print('Prepare Data Done')
#jpeg("test.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(2,2))
manhattan(chr="chr",bp="mid_pos",p="D1_pi", subset(dat, chr == 3), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "CS pi (500K bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D3_pi", subset(dat, chr == 3), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "NEH pi (500K bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D1_pi", subset(dat, chr == 5), highlight1 = h2$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "CS pi (500K bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D3_pi", subset(dat, chr == 5), highlight1 = h2$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "NEH pi (500K bp/window)", cex.main=1.5)
graph2ppt(file="pi_CS_NEH_500K.pptx", width=9, height=6.5)
#dev.off()
#  print("Pi Plotting done")

###################### start SL_OBOYS2 plot ######################
common_id = intersect(DT5$id5,DT4$id4)
common_id <- as.vector(common_id)
sub_DT5 = DT5[DT5$id5 %in% common_id,]
sub_DT4 = DT4[DT4$id4 %in% common_id,]

sub_DT5 = sub_DT5[match(common_id, sub_DT5$id5),]
sub_DT4 = sub_DT4[match(common_id, sub_DT4$id4),]

dat <- data.frame(chr=sub_DT5$CHROM, start=sub_DT5$BIN_START, start=sub_DT5$BIN_END, mid_pos=sub_DT5$mid_pos, SNP=common_id, D5_pi=sub_DT5$PI, D4_pi=sub_DT4$PI)
dat$chr <- as.numeric(dat$chr)
dat$mid_pos <- as.numeric(dat$mid_pos)
dat$D5_pi <- as.numeric(dat$D5_pi)
dat$D4_pi <- as.numeric(dat$D4_pi)
#write.table(dat, file = "500K_chr5_pi.csv", sep = ",", quote = FALSE,
#            row.names = FALSE, col.names = FALSE)
# import the highlight snps list
file1 = 'snpsOfInterest_chr3.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
file2 = 'snpsOfInterest_chr5.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)
print('Prepare Data Done')
#jpeg("test.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(2,2))
manhattan(chr="chr",bp="mid_pos",p="D5_pi", subset(dat, chr == 3), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "SL pi (500K bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D4_pi", subset(dat, chr == 3), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "OBOYS2 pi (500K bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D5_pi", subset(dat, chr == 5), highlight1 = h2$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "SL pi (500K bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D4_pi", subset(dat, chr == 5), highlight1 = h2$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "OBOYS2 pi (500K bp/window)", cex.main=1.5)
graph2ppt(file="pi_SL_OBOYS2_500K.pptx", width=9, height=6.5)
#dev.off()
#  print("Pi Plotting done")
#}


setwd("~/Documents/Ryan_workplace/CVreseq_dxy")
source("manhattan.R")
library(export)

# plot_pi <- function(name1, name2, plot_name, pop){
# load the original windowed pi
name1 = "CS_DEBY_chr8.csv"
DT1 = read.delim(name1, header = TRUE, sep=',')
mid_pos <- DT1$mid
id1 = paste0(DT1$scaffold,'_',mid_pos)
DT1 <- as.data.frame(cbind(DT1,mid_pos, id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]
dat1 <- data.frame(chr=DT1$scaffold, mid_pos=DT1$mid, SNP=DT1$id1, p1=DT1$pi_CS, p2=DT1$pi_DEBY)

name2 = "CS_NEH_chr8.csv"
DT2 = read.delim(name2, header = TRUE, sep=',')
mid_pos <-  DT2$mid
id2 = paste0(DT2$scaffold,'_',mid_pos)
DT2 <- as.data.frame(cbind(DT2,mid_pos,id2))
# remove NA rows in the dxy column
DT2 <- DT2[complete.cases(DT2), ]
dat2 <- data.frame(chr=DT2$scaffold, mid_pos=DT2$mid, SNP=DT2$id2, p1=DT2$pi_CS, p2=DT2$pi_NEH)

name3 = "SL_OBOYS2_chr8.csv"
DT3 = read.delim(name3, header = TRUE, sep=',')
mid_pos <-  DT3$mid
id3 = paste0(DT3$scaffold,'_',mid_pos)
DT3 <- as.data.frame(cbind(DT3,mid_pos,id3))
# remove NA rows in the dxy column
DT3 <- DT3[complete.cases(DT3), ]
dat3 <- data.frame(chr=DT3$scaffold, mid_pos=DT3$mid, SNP=DT3$id3, p1=DT3$pi_SL, p2=DT3$pi_OBOYS2)

file1 = 'snpsOfInterest_chr8.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
###################### start CS_DEBY plot ######################

# violin plot
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
library(ggplot2)

plotdat = data.frame(Pi=c(dat1$p1,dat1$p2), x=c(rep('CS',length(dat1$p1)), rep('DEBY',length(dat1$p2))))
p <- ggplot(plotdat, aes(x=x, y=Pi, fill=x)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.2)
p + scale_color_manual(values=c("#999999", "#E69F00"))
graph2ppt(file="pi_CS_DEBY_1K.pptx", width=9, height=6.5)


plotdat = data.frame(Pi=c(dat2$p1,dat2$p2), x=c(rep('CS',length(dat2$p1)), rep('NEH',length(dat2$p2))))
p <- ggplot(plotdat, aes(x=x, y=Pi, fill=x)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.2)
p + scale_color_manual(values=c("#999999", "#E69F00"))
graph2ppt(file="pi_CS_NEH_1K.pptx", width=9, height=6.5)

plotdat = data.frame(Pi=c(dat3$p1,dat3$p2), x=c(rep('SL',length(dat3$p1)), rep('OBOYS2',length(dat3$p2))))
p <- ggplot(plotdat, aes(x=x, y=Pi, fill=x)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.2)
p + scale_color_manual(values=c("#999999", "#E69F00"))
graph2ppt(file="pi_SL_OBOYS2_1K.pptx", width=9, height=6.5)


# Continusous colors
dp + scale_fill_brewer(palette="Blues") + theme_classic()
# Discrete colors
dp + scale_fill_brewer(palette="Dark2") + theme_minimal()
# Gradient colors
dp + scale_fill_brewer(palette="RdBu") + theme_minimal()


jpeg("CS_DEBY_pi_chr8.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="mid_pos",p="p1", subset(dat1, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.04), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "CS pi (1000 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="p2",subset(dat1, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.04), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "DEBY pi (1000 bp/window)", cex.main=1.5)
dev.off()
#  print("Pi Plotting done")

###################### start CS_NEH plot ######################
jpeg("CS_NEH_pi_chr8.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="mid_pos",p="p1", subset(dat2, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.04), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "CS pi (1000 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="p2",subset(dat2, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.04), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "NEH pi (1000 bp/window)", cex.main=1.5)
dev.off()
#  print("Pi Plotting done")

###################### start SL_OBOYS2 plot ######################
jpeg("SL_OBOYS2_chr8.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="mid_pos",p="p1", subset(dat3, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.04), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "SL pi (1000 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="p2",subset(dat3, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.04), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "OBOYS2 pi (1000 bp/window)", cex.main=1.5)
dev.off()
#  print("Pi Plotting done")

###################################### 100bp pi from SM output #####################################
###################################### 100bp pi from SM output #####################################
###################################### 100bp pi from SM output #####################################
source("manhattan.R")
library(export)

# plot_pi <- function(name1, name2, plot_name, pop){
# load the original windowed pi
name1 = "CS_DEBY_chr8_100bp.csv"
DT1 = read.delim(name1, header = TRUE, sep=',')
mid_pos <- DT1$mid
id1 = paste0(DT1$scaffold,'_',mid_pos)
DT1 <- as.data.frame(cbind(DT1,mid_pos, id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]
dat1 <- data.frame(chr=DT1$scaffold, mid_pos=DT1$mid, SNP=DT1$id1, p1=DT1$pi_CS, p2=DT1$pi_DEBY)

name2 = "CS_NEH_chr8_100bp.csv"
DT2 = read.delim(name2, header = TRUE, sep=',')
mid_pos <-  DT2$mid
id2 = paste0(DT2$scaffold,'_',mid_pos)
DT2 <- as.data.frame(cbind(DT2,mid_pos,id2))
# remove NA rows in the dxy column
DT2 <- DT2[complete.cases(DT2), ]
dat2 <- data.frame(chr=DT2$scaffold, mid_pos=DT2$mid, SNP=DT2$id2, p1=DT2$pi_CS, p2=DT2$pi_NEH)

name3 = "SL_OBOYS2_chr8_100bp.csv"
DT3 = read.delim(name3, header = TRUE, sep=',')
mid_pos <-  DT3$mid
id3 = paste0(DT3$scaffold,'_',mid_pos)
DT3 <- as.data.frame(cbind(DT3,mid_pos,id3))
# remove NA rows in the dxy column
DT3 <- DT3[complete.cases(DT3), ]
dat3 <- data.frame(chr=DT3$scaffold, mid_pos=DT3$mid, SNP=DT3$id3, p1=DT3$pi_SL, p2=DT3$pi_OBOYS2)

write.table(DT1, file = "SM_chr8_pi_100bp.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

file1 = 'snpsOfInterest_chr8_100bp.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
###################### start CS_DEBY plot ######################

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
library(ggplot2)

jpeg("CS_DEBY_pi_chr8_100bp_violin.jpg", width = 6, height = 4, units = 'in', res = 300)
plotdat = data.frame(Pi=c(dat1$p1,dat1$p2), Population=c(rep('CS',length(dat1$p1)), rep('DEBY',length(dat1$p2))))
p <- ggplot(plotdat, aes(x=Population, y=Pi, fill=Population)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.2)
p + scale_color_manual(values=c("#999999", "#E69F00")) + theme(axis.text.x = element_text(size=12),axis.title.x = element_text(size=15), axis.text.y = element_text(size=12), axis.title.y = element_text(size=15))
#graph2ppt(file="pi_CS_DEBY_100bp.pptx", width=9, height=6.5)
dev.off()

jpeg("CS_NEH_pi_chr8_100bp_violin.jpg", width = 6, height = 4, units = 'in', res = 300)
plotdat = data.frame(Pi=c(dat2$p1,dat2$p2), Population=c(rep('CS',length(dat2$p1)), rep('NEH',length(dat2$p2))))
p <- ggplot(plotdat, aes(x=Population, y=Pi, fill=Population)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.2)
p + scale_color_manual(values=c("#999999", "#E69F00")) + theme(axis.text.x = element_text(size=12),axis.title.x = element_text(size=15), axis.text.y = element_text(size=12), axis.title.y = element_text(size=15))
#graph2ppt(file="pi_CS_NEH_100bp.pptx", width=9, height=6.5)
dev.off()

jpeg("SL_OBOYS2_pi_chr8_100bp_violin.jpg", width = 6, height = 4, units = 'in', res = 300)
plotdat = data.frame(Pi=c(dat3$p1,dat3$p2), Population=c(rep('SL',length(dat3$p1)), rep('OBOYS2',length(dat3$p2))))
p <- ggplot(plotdat, aes(x=Population, y=Pi, fill=Population)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.2)
p + scale_color_manual(values=c("#999999", "#E69F00")) + theme(axis.text.x = element_text(size=12),axis.title.x = element_text(size=15), axis.text.y = element_text(size=12), axis.title.y = element_text(size=15))
#graph2ppt(file="pi_SL_OBOYS2_100bp.pptx", width=9, height=6.5)
dev.off()

print("Start Pi Plotting")

jpeg("CS_DEBY_pi_chr8_100bp.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="mid_pos",p="p1", subset(dat1, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "CS pi (100 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="p2",subset(dat1, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "DEBY pi (100 bp/window)", cex.main=1.5)
dev.off()

jpeg("CS_NEH_pi_chr8_100bp.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="mid_pos",p="p1", subset(dat2, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "CS pi (100 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="p2",subset(dat2, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "NEH pi (100 bp/window)", cex.main=1.5)
dev.off()

jpeg("SL_OBOYS2_chr8_100bp.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="mid_pos",p="p1", subset(dat3, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "SL pi (100 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="p2",subset(dat3, chr == 8), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.01), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(paste(pi)), cex.lab=1.5, main = "OBOYS2 pi (100 bp/window)", cex.main=1.5)
dev.off()
print("Pi Plotting done")


