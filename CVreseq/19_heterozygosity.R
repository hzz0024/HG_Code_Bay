# load the R function
setwd("~/Documents/Ryan_workplace/CVreseq_He")
source("manhattan.R")
library(export)

#plot_pi <- function(name1, name2, plot_name, pop){
# load the original windowed pi
name1 = "CS"
DT1 = read.delim(paste0(name1,"_He.txt"), header = FALSE, sep='\t')
id1 = paste0(DT1$V1,'_',DT1$V2)
DT1 <- as.data.frame(cbind(DT1,id1))
name2 = "HC"
DT2 = read.delim(paste0(name2,"_He.txt"), header = FALSE, sep='\t')
id2 = paste0(DT2$V1,'_',DT2$V2)
DT2 <- as.data.frame(cbind(DT2,id2))
common_id = intersect(DT1$id1,DT2$id2)
common_id <- as.vector(common_id)
sub_DT1 = DT1[DT1$id1 %in% common_id,]
sub_DT2 = DT2[DT2$id2 %in% common_id,]

dat <- data.frame(chr=sub_DT1$V1, pos=sub_DT1$V2, SNP=common_id, D1_He=sub_DT1$V5, D2_He=sub_DT2$V5, diff_H=sub_DT1$V5-sub_DT2$V5)
dat$chr <- as.numeric(dat$chr)
dat$pos <- as.numeric(dat$pos)
dat$diff_H <- as.numeric(dat$diff_dxy)
#write.table(dat, file = "chr1_HC_CLP-CS_HCVA.csv", sep = ",", quote = FALSE,
#            row.names = FALSE, col.names = FALSE)
# import the highlight snps list
#file1 = 'snpsOfInterest.csv' 
#h1 = read.delim(file1, header = FALSE, sep=',')
#h1 = as.list(h1)
#print('Prepare Data Done')
###################### start He plot ######################
jpeg(paste0(name1,"-",name2,"_He.jpg"), width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(chr="chr",bp="pos",p="diff_H", subset(dat, chr == 2),  logp=FALSE, 
          cex.axis = 1, ylim = c(-1, 1), xlim = c(33400000, 33750000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F, #highlight1 = h1$V1,
          ylab="He difference", xlab="", cex.lab=1.5, main = "CS-HC Observed He", cex.main=1.5)
#graph2ppt(file="Dxy_1K.pptx", width=9, height=6.5)
dev.off()
  print("Pi Plotting done")
#}

  
###################### make igv file ######################
make_bedgraph <- function(name1, name2){
  DT1 = read.delim(paste0(name1,"_He.txt"), header = FALSE, sep='\t')
  id1 = paste0(DT1$V1,'_',DT1$V2)
  DT1 <- as.data.frame(cbind(DT1,id1))
  DT2 = read.delim(paste0(name2,"_He.txt"), header = FALSE, sep='\t')
  id2 = paste0(DT2$V1,'_',DT2$V2)
  DT2 <- as.data.frame(cbind(DT2,id2))
  common_id = intersect(DT1$id1,DT2$id2)
  common_id <- as.vector(common_id)
  sub_DT1 = DT1[DT1$id1 %in% common_id,]
  sub_DT2 = DT2[DT2$id2 %in% common_id,]
  dat <- data.frame(chr=sub_DT1$V1, pos=sub_DT1$V2, SNP=common_id, D1_He=sub_DT1$V5, D2_He=sub_DT2$V5, diff_H=sub_DT1$V5-sub_DT2$V5)
  dat$chr <- as.numeric(dat$chr)
  dat$pos <- as.numeric(dat$pos)
  dat$diff_H <- as.numeric(dat$diff_H)
  dat.bedgraph <- data.frame(Chromosome=as.numeric(dat$chr), Start=as.numeric(dat$pos)-1, End=as.numeric(dat$pos), Diff_H=as.numeric(dat$diff_H))
  dat.bedgraph <- dat.bedgraph[complete.cases(dat.bedgraph), ]
  dat.bedgraph$Chromosome <- as.numeric(as.character(dat.bedgraph$Chromosome))
  dat.bedgraph$Start <- as.numeric(as.character(dat.bedgraph$Start))
  dat.bedgraph$End <- as.numeric(as.character(dat.bedgraph$End))
  dat.bedgraph$Diff_H <- as.numeric(as.character(dat.bedgraph$Diff_H))
  dat.sort.bedgraph = dat.bedgraph[order(dat.bedgraph$Chromosome, dat.bedgraph$Start),]
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    dat.sort.bedgraph$Chromosome[dat.sort.bedgraph$Chromosome==i] = chr_str_list[i]
  write.table(dat.sort.bedgraph, file = paste0(name1,"-",name2,"_diff_H.sort.bedgraph"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

make_bedgraph("CS", "HC")
make_bedgraph("CS", "HCVA")
make_bedgraph("HC", "CLP")
make_bedgraph("HCVA", "CLP")
make_bedgraph("CS", "NEH")
make_bedgraph("CS", "DEBY")
make_bedgraph("CS", "UMFS")
  