# load the R function
setwd("~/Documents/Ryan_workplace/CVreseq_zoom_in")
source("manhattan.R")
library(KRIS)
options(scipen=999)
library(export)


format <- function(prefix){
  bed = paste0("zoom_in_",prefix,"_fst.bed")
  bim = paste0("zoom_in_",prefix,"_fst.bim")
  fam = paste0("zoom_in_",prefix,"_fst.fam")
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
  #chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  #for(i in seq(10)) 
  #  dat.sort.bedgraph$Chromosome[dat.sort.bedgraph$Chromosome==i] = chr_str_list[i]
  write.table(dat.sort.bedgraph, file = paste0(prefix, ".sort.bedgraph"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

format("CS_HC")
format("HCVA_CLP")


file1 = 'SNP_of_interest_CS_HC.txt' 
h1 = read.delim(file1, header = FALSE, sep='\t')
h1 = as.list(h1)

name1 = "CS_HC.sort.bedgraph"
DT1 = read.delim(name1, header = FALSE, sep='\t')
id1 = paste0(DT1$V1,'_',DT1$V3)
DT1 <- as.data.frame(cbind(DT1, id1))
DT1 <- DT1[complete.cases(DT1), ]
DT1[,4][DT1[,4]<0] = 0 #@chnage
dat1 <- data.frame(chr=DT1$V1, pos=DT1$V2, SNP=DT1$id, fst=DT1[,4])
dat1$chr <- as.numeric(dat1$chr)
dat1$pos <- as.numeric(dat1$pos)
dat1$fst <- as.numeric(dat1$fst)


file2 = 'SNP_of_interest_HCVA_CLP.txt' 
h2 = read.delim(file2, header = FALSE, sep='\t')
h2 = as.list(h2)

name2 = "HCVA_CLP.sort.bedgraph"
DT2 = read.delim(name2, header = FALSE, sep='\t')
id2 = paste0(DT2$V1,'_',DT2$V3)
DT2 <- as.data.frame(cbind(DT2, id2))
DT2 <- DT2[complete.cases(DT2), ]
DT2[,4][DT2[,4]<0] = 0 #@chnage
dat2 <- data.frame(chr=DT2$V1, pos=DT2$V2, SNP=DT2$id, fst=DT2[,4])
dat2$chr <- as.numeric(dat2$chr)
dat2$pos <- as.numeric(dat2$pos)
dat2$fst <- as.numeric(dat2$fst)

jpeg("zoom_in1.jpg", width = 18, height = 10, units = 'in', res = 300)
par(mar=c(5,6,4,4)+.1)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="pos",p="fst", subset(dat1, chr == 2), xlim =c(33580000, 33650000),  highlight1 = h1$V1, 
          logp=FALSE, cex.axis = 1.5, ylim = c(0, 1.02), #
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=NA, cex.lab=3, cex.main=3)
require(splines)
a = dat1$chr==2
b = dat1$pos>33580000
c = dat1$pos<33650000
dat11 = dat1[a&b&c,]
fit1<-lm(fst ~ bs(pos),data = dat11)
summary(fit1)
xlims<-range(dat11$pos)
pos.grid1<-seq(from=xlims[1], to = xlims[2])
points(pos.grid,predict(fit1,newdata = list(pos=pos.grid1)),col="black",lwd=1,type="l")
abline(h=0, col="red", lwd=0.6, lty=2)

manhattan(chr="chr",bp="pos",p="fst", subset(dat2, chr == 2), xlim =c(33580000, 33650000), highlight2 = h2$V1, 
          logp=FALSE, cex.axis = 1.5, ylim = c(0, 1.02), #
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=NA, cex.lab=3, cex.main=3)
require(splines)
a = dat2$chr==2
b = dat2$pos>33580000
c = dat2$pos<33650000
dat22 = dat2[a&b&c,]
fit2<-lm(fst ~ bs(pos),data = dat22 )
summary(fit2)
xlims<-range(dat22$pos)
pos.grid1<-seq(from=xlims[1], to = xlims[2])
points(pos.grid,predict(fit2,newdata = list(pos=pos.grid1)),col="black",lwd=1,type="l")
abline(h=0, col="red", lwd=0.6, lty=2)
dev.off()

##### plot with average H difference ##### 
jpeg("zoom_in4.jpg", width = 18, height = 10, units = 'in', res = 300)
par(mar=c(5,6,4,4)+.1)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="pos",p="fst", subset(dat1, chr == 2), xlim =c(33580000, 33650000),  highlight1 = h1$V1, 
          logp=FALSE, cex.axis = 1.5, ylim = c(0, 1.02), #
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=NA, cex.lab=3, cex.main=3)
require(splines)
a = dat1$chr==2
b = dat1$pos>33580000
c = dat1$pos<33650000
dat11 = dat1[a&b&c,]
intervals = seq(33580000, 33650000, 1000)
xs = c()
ys = c()
for(i in seq(length(intervals)-1)){
  idxs = dat11$pos>=intervals[i] & dat11$pos<=intervals[i+1]
  if(sum(idxs)>0){
    ys = c(ys, mean(dat11[idxs,]$fst))
    xs = c(xs, (intervals[i]+intervals[i+1])/2)
  }
}
lines(xs, ys)

manhattan(chr="chr",bp="pos",p="fst", subset(dat2, chr == 2), xlim =c(33580000, 33650000), highlight2 = h2$V1, 
          logp=FALSE, cex.axis = 1.5, ylim = c(0, 1.02), #
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=NA, cex.lab=3, cex.main=3)
require(splines)
a = dat2$chr==2
b = dat2$pos>33580000
c = dat2$pos<33650000
dat22 = dat2[a&b&c,]
xs = c()
ys = c()
for(i in seq(length(intervals)-1)){
  idxs = dat22$pos>=intervals[i] & dat22$pos<=intervals[i+1]
  if(sum(idxs)>0){
    ys = c(ys, mean(dat22[idxs,]$fst))
    xs = c(xs, (intervals[i]+intervals[i+1])/2)
  }
}
lines(xs, ys)

dev.off()


#################### He diff ####################

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
  write.table(dat.sort.bedgraph, file = paste0(name1,"-",name2,"_diff_H.sort.bedgraph"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

make_bedgraph("CS", "HC")
make_bedgraph("HCVA", "CLP")

###################### start He plot ######################
file1 = 'SNP_of_interest_CS_HC.txt' 
h1 = read.delim(file1, header = FALSE, sep='\t')
h1 = as.list(h1)

name = "CS-HC_diff_H.sort.bedgraph"
DT = read.delim(name, header = FALSE, sep='\t')
id = paste0(DT$V1,'_',DT$V3)
DT <- as.data.frame(cbind(DT, id))
dat3 <- data.frame(chr=DT$V1, pos=DT$V2, SNP=DT$id, Hd=DT[,4])

file2 = 'SNP_of_interest_HCVA_CLP.txt' 
h2 = read.delim(file2, header = FALSE, sep='\t')
h2 = as.list(h2)

name = "HCVA-CLP_diff_H.sort.bedgraph"
DT = read.delim(name, header = FALSE, sep='\t')
id = paste0(DT$V1,'_',DT$V3)
DT <- as.data.frame(cbind(DT, id))
dat4 <- data.frame(chr=DT$V1, pos=DT$V2, SNP=DT$id, Hd=DT[,4])


##### plot with spline ##### 
jpeg("zoom_in2.jpg", width = 18, height = 10, units = 'in', res = 300)
par(mar=c(5,6,4,4)+.1)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="pos",p="Hd", subset(dat3, chr == 2), xlim = c(33580000, 33650000),  highlight1 = h1$V1, 
          logp=FALSE, cex.axis = 1.5, ylim = c(-1.02, 1.02), #
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=NA, cex.lab=3, cex.main=3) #"Heterozygosity Difference"
# fitting the data with spline (linear regression with lm function)
require(splines)
a = dat3$chr==2
b = dat3$pos>33580000
c = dat3$pos<33650000
dat33 = dat3[a&b&c,]
fit1<-lm(Hd ~ bs(pos),data = dat33 )
summary(fit1)
xlims<-range(dat33$pos)
pos.grid1<-seq(from=xlims[1], to = xlims[2])
points(pos.grid,predict(fit1,newdata = list(pos=pos.grid1)),col="black",lwd=1,type="l")
abline(h=0, col="red", lwd=0.6, lty=2)

manhattan(chr="chr",bp="pos",p="Hd", subset(dat4, chr == 2), xlim = c(33580000, 33650000),  highlight2 = h2$V1, 
          logp=FALSE, cex.axis = 1.5, ylim = c(-1.02, 1.02), #
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=NA, cex.lab=3, cex.main=3)
#graph2ppt(file="fst_He_outlier",width=6,height=10)
require(splines)
a = dat4$chr==2
b = dat4$pos>33580000
c = dat4$pos<33650000
dat44 = dat4[a&b&c,]
fit2<-lm(Hd ~ bs(pos),data = dat44 )
summary(fit2)
xlims<-range(dat44$pos)
pos.grid2<-seq(from=xlims[1], to = xlims[2])
points(pos.grid,predict(fit2,newdata = list(pos=pos.grid2)),col="black",lwd=1,type="l")
abline(h=0, col="red", lwd=0.6, lty=2)

dev.off()

##### plot with average H difference ##### 
#calculate mean every 1000bp
jpeg("zoom_in3.jpg", width = 18, height = 10, units = 'in', res = 300)
par(mar=c(5,6,4,4)+.1)
par(mfrow=c(2,1))
manhattan(chr="chr",bp="pos",p="Hd", subset(dat3, chr == 2), xlim = c(33580000, 33650000),  highlight1 = h1$V1, 
          logp=FALSE, cex.axis = 1.5, ylim = c(-1.02, 1.02), #
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=NA, cex.lab=3, cex.main=3) #"Heterozygosity Difference"
a = dat3$chr==2
b = dat3$pos>33580000
c = dat3$pos<33650000
dat33 = dat3[a&b&c,]
intervals = seq(33580000, 33650000, 1000)
xs = c()
ys = c()
for(i in seq(length(intervals)-1)){
  idxs = dat33$pos>=intervals[i] & dat33$pos<=intervals[i+1]
  if(sum(idxs)>0){
    ys = c(ys, mean(dat33[idxs,]$Hd))
    xs = c(xs, (intervals[i]+intervals[i+1])/2)
  }
}
lines(xs, ys)
abline(h=0, col="red", lwd=0.6, lty=2)
manhattan(chr="chr",bp="pos",p="Hd", subset(dat4, chr == 2), xlim = c(33580000, 33650000),  highlight2 = h2$V1, 
          logp=FALSE, cex.axis = 1.5, ylim = c(-1.02, 1.02), #
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=NA, cex.lab=3, cex.main=3)
a = dat4$chr==2
b = dat4$pos>33580000
c = dat4$pos<33650000
dat44 = dat4[a&b&c,]
intervals = seq(33580000, 33650000, 1000)
xs = c()
ys = c()
for(i in seq(length(intervals)-1)){
  idxs = dat44$pos>=intervals[i] & dat44$pos<=intervals[i+1]
  if(sum(idxs)>0){
    ys = c(ys, mean(dat44[idxs,]$Hd))
    xs = c(xs, (intervals[i]+intervals[i+1])/2)
  }
}
lines(xs, ys)
abline(h=0, col="red", lwd=0.6, lty=2)
dev.off()



# jpeg("zoom_in.jpg", width = 18, height = 9, units = 'in', res = 300)
# par(mar=c(5,6,4,4)+.1)
# par(mfrow=c(4,1))
# manhattan(chr="chr",bp="pos",p="fst", subset(dat1, chr == 2), xlim =c(33530000, 33700000),  highlight1 = h1$V1, 
#           logp=FALSE, cex.axis = 1.5, ylim = c(0, 1.02), #
#           col=c("grey","black"),genomewideline=F, suggestiveline=F,
#           ylab="Fst", cex.lab=3, cex.main=3)
# manhattan(chr="chr",bp="pos",p="fst", subset(dat2, chr == 2), xlim =c(33530000, 33700000), highlight2 = h2$V1, 
#           logp=FALSE, cex.axis = 1.5, ylim = c(0, 1.02), #
#           col=c("grey","black"),genomewideline=F, suggestiveline=F,
#           ylab=NA, cex.lab=3, cex.main=3)
# #jpeg("zoom_in2.jpg", width = 18, height = 6, units = 'in', res = 300)
# manhattan(chr="chr",bp="pos",p="Hd", subset(dat3, chr == 2), xlim = c(33530000, 33700000),  highlight1 = h1$V1, 
#           logp=FALSE, cex.axis = 1.5, ylim = c(-1.02, 1.02), #
#           col=c("grey","black"),genomewideline=F, suggestiveline=F,
#           ylab=expression(Delta~Heterozygosity), cex.lab=3, cex.main=3)
# manhattan(chr="chr",bp="pos",p="Hd", subset(dat4, chr == 2), xlim = c(33530000, 33700000),  highlight2 = h2$V1, 
#           logp=FALSE, cex.axis = 1.5, ylim = c(-1.02, 1.02), #
#           col=c("grey","black"),genomewideline=F, suggestiveline=F,
#           ylab=NA, cex.lab=3, cex.main=3)
# #graph2ppt(file="fst_He_outlier",width=6,height=10)
# dev.off()

################################## genotype plot ####################################
#install.packages("remotes")
#remotes::install_github("JimWhiting91/genotype_plot")
library(devtools)
#install('/Users/ryan/Downloads/genotype_plot-master')
#install.packages("glue")
library(glue)
library(vcfR)
library(ggplot2)
library(GenotypePlot)
library(cowplot)
library(rlang)
library(export)

#my_vcf <- read.vcfR("CS_HC-HCVA_CLP.outlier.vcf.gz")
pop_name <- "CS_HC-HCVA_CLP.txt"
popmap = read.delim(pop_name, header = TRUE, sep='\t')
new_plot <- genotype_plot(vcf    =  "CS_HC-HCVA_CLP.outlier.vcf.gz",
                          chr    = 2,
                          start  = 33597600,  
                          end    = 33629600,
                          popmap = popmap,                              
                          cluster        = TRUE,                           
                          snp_label_size = 10000,                          
                          colour_scheme=c("#820041","#9782AD","#EFD3B5"))   

cowplot::plot_grid(new_plot$positions, new_plot$genotypes, axis="tblr",
                   align="v", nrow=2, ncol=1, rel_heights=c(1,9))

# Add dendrogram tips with points
Population = popmap$pop
group = strsplit(new_plot$dendro_labels, "_\\s*(?=[^_]+$)", perl=TRUE)
group = sapply(group, "[[", 1)
Origin = c(rep('High salinity', 24))
#Origin[group %in% Population[1:12]] = 'Del'
Origin[group %in% Population[7:12]] = 'Low salinity'
Origin[group %in% Population[19:24]] = 'Low salinity'
#tmp = data.frame(group, Origin)
Population = group
df = data.frame(x=1:length(new_plot$dendro_labels), y=rep(-1,24))
dendro_with_tips <- new_plot$dendrogram +
  geom_point(aes(x=df$x, y=df$y, color=Population, shape=Origin), size = 5) + #c("#1F78B4", "#33A02C", "#86a395","#8a584a")
  theme(legend.position="left") +
  theme(legend.key = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0.6, "cm"),
        legend.key.width = unit(0.8,"cm"))+
  scale_color_manual(values = c("#33A02C", "#FB9A99", "#E31A1c","#B2DF8A")) +
  scale_shape_manual(values=c(15, 17))+
  theme(legend.title = element_text(size = 12),
        legend.text=element_text(size=12))

dendro_with_tips
legend <- get_legend(dendro_with_tips)  
dendro_with_tips <- dendro_with_tips + theme(legend.position='none')
# For e.g. plotting dendrogram and clustered genotypes:
geno_and_dendro <- cowplot::plot_grid(#new_plot$dendrogram,
  dendro_with_tips, 
  new_plot$genotypes,
  axis="tblr",align="h",nrow=1,ncol=2,rel_widths=c(3,7)) 
plot_grid(geno_and_dendro, legend, rel_widths = c(5, 0.2))

graph2ppt(file="CS_HC-HCVA_CLP_outlier",width=12,height=10)






jpeg("gene.jpg", width = 18, height = 10, units = 'in', res = 300)

plot(1, type="n", xlab="", ylab="", yaxt='n', xlim=c(33597600, 33629600), ylim=c(0, 10))
#plot(1, type="n", xlab="", ylab="", yaxt='n', xlim=c(33580000, 33650000), ylim=c(0, 10))
fill <- function(x1, x2, color){
  delta = 10
  for(i in seq((x2-x1)/delta)){
    lines(x1+i*delta+c(0,0), 1+c(0+0.05,1-0.05), col=color)
  }
}

filename = '3gene_exon.txt'
fill_dat = read.delim(filename, sep='\t', header=F)
X = fill_dat$V4
Y = fill_dat$V5
for(i in seq(length(X))){
  fill(X[i], Y[i], 'red')
}

X = c(33603020, 33610182, 33582628)
Y = c(33609907, 33630413, 33601234)
for(i in seq(length(X))){
  lines(X[i]+c(0,Y[i]-X[i],Y[i]-X[i],0,0),1+c(0,0,1,1,0))
}
dev.off()
#graph2ppt(file="gene.pptx", width=9, height=6.5)

