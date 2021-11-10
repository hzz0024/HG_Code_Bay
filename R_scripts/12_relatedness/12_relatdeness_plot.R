####### 

# script used for relatedness estimation
# #!/bin/sh
# ngsRelate=/programs/ngsRelate-20190304/ngsRelate
# Del19_challenge_cnt=97
# Del20_challenge_cnt=101
# HC_cnt=47
# ARN_cnt=47
# COH_cnt=44
# SR_cnt=48
# NB_cnt=48
# 
# for pop in Del20_challenge HC ARN COH SR NB; do # Del19_challenge
# # extract the snp frequency data from mafs zip file
# zcat $pop'_minq20_minmq30_1x_CV30_masked.mafs.gz' | cut -f6 |sed 1d > $pop'_freq'
# cnt=$pop'_cnt'
# # calculate the stats using ngsrelate
# echo ${!cnt}
# $ngsRelate -f $pop'_freq' -g $pop'_minq20_minmq30_1x_CV30_masked.glf.gz' -n ${!cnt} -z $pop'.list' -O $pop'.res'
# done

############################### below are analysis across all sequenced bam files ##############################
####### format data and produce the csv file #######
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/12_relatedness")
library(data.table)
files <- list.files(pattern="*.res", full.names=TRUE, recursive=FALSE)
df_all = data.frame(ida = character(), idb = character(), rab = numeric())
for (x in files) {
  t <- read.delim(x, header = TRUE, sep='\t')
  df <- subset(t, select=c("ida", "idb" , "rab"))
  write.table(df, paste0(x, ".txt"), sep="\t", quote=F, row.names=F, col.names=F)
}


####### plot the het for different minq results #######
library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)
setwd("~/Documents/Ryan_workplace/DelBay_adult/12_relatedness")
file = "Summary_all.txt"
df <- read.delim(file, header = TRUE, sep='\t')
ggplot(df, aes(x=rab, y=Pop, shape=Pop)) +
  geom_point()

# reroder the population
df$Pop <- factor(df$Pop , levels=c("HC", "ARN", "COH", "SR", "NB", "Survivor", "Reference"))
# set up colors for plots
cbPalette <- c("#A71B4B", "#E97302", "#EAC728", "#0BC0B3", "#4461A8", "#F5191C", "#7A7A7A")

p1 <- ggplot(df, aes(x=Pop, y=rab, fill=as.factor(Pop))) +
  #ylim(0,0.0001)+
  # color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.5), outlier.alpha = 0.3)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="right")+
  labs(fill="Population")+
  ylab("Pairwise relatedness") +
  xlab("") +
  scale_fill_manual(values=cbPalette, breaks=c("HC", "ARN", "COH", "SR", "NB", "Survivor", "Reference")) ## @HG changed the name for your own plots
  
p1

graph2ppt(file="relatedness",width=18,height=9)









######### below are codes from original ngsrelate git (with minor edits) #############

args <- commandArgs(trailing=T)
df <- read.table(file = "Del20_challenge.res",h=T)

## df <- head(df[order(df$rab,decreasing=T),],n=10)
## df <- head(df[order(df$rab,decreasing=T),],n=2000)

# head(df,n=200)

## use ids if present
if ("ida" %in% colnames(df)){
  namea = 3
  nameb = 4
} else {
  namea = 1
  nameb = 2
}

df$names <- paste(df[,namea], df[,nameb],sep="-")

par.old <- par()
threshold <- 0.00
######################
## RELATEDNESS PLOT ##
######################
min.rab.threshold <- threshold
df.rab <- df[df$rab>=min.rab.threshold, ]

rab.height=nrow(df.rab)/30
if(rab.height<7)
  rab.height=7

rab.order <- order(df.rab$rab, decreasing=T)
#pdf(paste0(args[1], ".relatedness.pdf"), useDingbats=F, height=rab.height,width=2)
par(lty=0,xaxs="i",yaxs="i")
#names.arg=df.rab$names[rab.order], cex.names=0.2
xpos <- barplot(df.rab$rab[rab.order], horiz=T, cex.axis=0.3, cex.names=0.2, las=1, xlim=c(0,.2),space = 0.01, main="Relatedness")
#axis(side=1,cex=0.2)
axis(3, cex.axis=0.3, las=2)
mtext(text=df.rab$names[rab.order], line=0.1,side=2,at=xpos, cex=0.25,las=2)
abline(v = c(0.0, 0.05, 0.1, 0.25, 0.5, 1),lwd=0.1,lty='dotted')

library(export)
graph2ppt(file="chr_relatedness",width=4,height=6)

dev.off()


##############
## 2of3 IBD ##
##############

min.two.of.three.IBD.threshold <- threshold
df.two.of.three.IBD <- df[df[,"X2of3_IDB"]>=min.two.of.three.IBD.threshold, ]

two.of.three.IBD.height <- nrow(df.two.of.three.IBD)/30
if(two.of.three.IBD.height<7)
  two.of.three.IBD.height=7

two.of.three.IBD.order <- order(df.two.of.three.IBD[,"X2of3_IDB"], decreasing=T)
pdf(paste0(args[1], ".two_of_three_IBD.pdf"), useDingbats=F, height=two.of.three.IBD.height,width=2)
par(lty=0,xaxs="i",yaxs="i")
## names.arg=df.two.of.three.IBD$names[two.of.three.IBD.order],cex.names=0.2
xpos <- barplot(df.two.of.three.IBD[,"X2of3_IDB"][two.of.three.IBD.order], horiz=T, cex.axis=0.3, las=2, xlim=c(0,1),space = 0.01, main="2 of 3 alleles IBD")
axis(3, cex.axis=0.3, las=2)
mtext(text=df.two.of.three.IBD$names[two.of.three.IBD.order], line=0,side=2,at=xpos, cex=0.2,las=2)
abline(v = c(0.0, 0.05, 0.1, 0.25, 0.5, 1),lwd=0.1,lty='dotted')
dev.off()


################
## INBREEDING ##
################
df.inbreeding <- data.frame(name=as.character(c(as.character(df[,namea]), as.character(df[,nameb]))), F=c(df$Fa, df$Fb))
inbreeding.height <- length(unique(df.inbreeding$name))/10
if (inbreeding.height<7)
  inbreeding.height <- 7
l <- list()
for (name in unique(df.inbreeding$name)){
  l[[name]] <- df.inbreeding[df.inbreeding$name == name,"F"]
}
l <- l[order(sapply(l, function(x){median(x)}),decreasing=T)]
pdf(paste0(args[1], ".inbreeding.pdf"), useDingbats=F, height=inbreeding.height, width=2)
boxplot(l,horizontal=T, las=2, space=0.01, cex.axis=0.3, cex=0.2, ylim=c(0,1),space=0.01, boxlwd=1, main="Inbreeding\nCoefficient")
abline(v = c(0.0, 0.05, 0.1, 0.25, 0.5, 1),lwd=0.1,lty='dotted')
dev.off()


################
## R0 R1 KING ##
################
## df = df[df$Fa<0.1 & df$Fb < 0.1,]
pdf(paste0(args[1], ".R0_R1_KING.pdf"), useDingbats=F) ## , height=inbreeding.height, width=2)
par(mfrow=c(2,2))
plot(df$R1, df$R0, pch=19, cex=0.1, xlab="R1", ylab="R0")
plot(df$R1, df$KING, pch=19,cex=0.1, xlab="R1", ylab="KING")
plot(df$R1, df$R0, xlim=c(0,1), ylim=c(-0.2,0.5), pch=19, cex=0.1, xlab="R1", ylab="R0")
plot(df$R1, df$KING, xlim=c(0,1), ylim=c(-0.2, 0.3), pch=19,cex=0.1, xlab="R1", ylab="KING")
dev.off()
