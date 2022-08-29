library(gtools)
library(reshape2)
library(LDheatmap)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(Rmisc)
library(export)

####################################
##########  plot LD block   ########
####################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD")
TMP_FILE = "NB.ngsld.output"

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/15_LD_outlier_snps/")
TMP_FILE = "HC_19.ngsld.output"

df <- read.table(TMP_FILE, header=FALSE, stringsAsFactors=FALSE)
colnames(df)=c('snp1', 'snp2', 'dis', 'r1', 'D1', 'D2', 'r')
r <- df[complete.cases(df), ]
id <- unique(mixedsort(c(r[,"snp1"],r[,"snp2"])))
posStart <- head(id,1)
posEnd <- tail(id,1)
r <- rbind(r, c(posStart,posStart,0,NA,NA,NA,NA), c(posEnd,posEnd,0,NA,NA,NA,NA))

SNPs = read.delim("chr5_ADF1.txt", header = FALSE, sep='\t')

for (ld in c("r")) {
  m <- apply(acast(r, snp1 ~ snp2, value.var=ld, drop=FALSE),2,as.numeric)
  rownames(m) <- colnames(m)
  m <- m[mixedorder(rownames(m)),mixedorder(colnames(m))]
  id <- rownames(m)
  dist <- as.numeric(sub(".*:","",id))
  
  # Save plot
  tiff(paste("LD_blocks", ld,"jpg", sep="."), units="in", width=6, height=6, res=300)
  #pdf(paste("LD_blocks", ld,"pdf", sep="."), width=10, height=10)
  LDheatmap(m, genetic.distances=dist, geneMapLabelX=0.75, geneMapLabelY=0.25, title = NULL, color="blueToRed", LDmeasure=ld, SNP.name = SNPs$V1)
  require(grid)
  grid.edit("symbols", pch = 20, gp = gpar(cex = 1, col = "red")) # change the symbol size
  grid.edit(gPath("ldheatmap", "geneMap","SNPnames"), gp = gpar(cex=0.2, col = "red"))  # change the SNP label size
  dev.off()
}

#########################################
##########  plot LD distribution ########
#########################################
library(ggpubr)
library(cowplot)
library(Rmisc)
library(export)


# calculate means
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/15_LD_outlier_snps/")
cal_mean <- function(pname){
  #pname1 = "./Challenge_20.ngsld.output"
  dat1 = read.delim(pname, header = FALSE, sep='\t')
  dat1 = dat1[,c(3,7)]
  dat1 = na.omit(dat1)
  dat1 = dat1[which(dat1$V3 < 500),]
  message(paste0("The mean of r2 in ", pname, ":", mean(dat1$V7)))
}

for (pop in c("Challenge_19", "Challenge_20", 
              "Sur_19", "Ref_19", "Sur_20", "Ref_20",
              "HC_18", "HC_19", "HC_21", 
              "COH_18", "COH_19", "COH_21",
              "ARN_18", "ARN_19", "ARN_21", 
              "SR_18", "SR_19", "SR_21",
              "NB_18", "NB_19", "NB_21", 
              "Wild_18", "Wild_21", "Wild_19")){
  cal_mean(paste0(pop, ".ngsld.output"))
}


setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/15_LD_GEA/")
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/")

format <- function (pname1, pname2){
  #pname1 = "./15_LD_outlier_snps/Challenge_19.ngsld.output"
  dat1 = read.delim(pname1, header = FALSE, sep='\t')
  dat1 = dat1[,c(3,7)]
  dat1 = na.omit(dat1)
  dat1 = dat1[which(dat1$V3 < 500),]
  #pname2 = "./15_LD_random_snps/HC_19.97.ngsld.output"
  dat2 = read.delim(pname2, header = FALSE, sep='\t')
  dat2 = dat2[,c(3,7)]
  dat2 = na.omit(dat2)
  dat2 = dat2[which(dat2$V3 < 500),]
  if(dim(dat1)[1] < dim(dat2)[1]){
    idx = seq(1:dim(dat2)[1])
    random = sort(sample(idx, dim(dat1)[1]))
    dat2 = dat2[random,]
  } else {
    idx = seq(1:dim(dat1)[1])
    random = sort(sample(idx, dim(dat2)[1]))
    dat1 = dat1[random,]
  }
  df_total = rbind(dat1, dat2)
  df_total$group = rep(c("outlier", "random"), c(dim(dat1)[1], dim(dat2)[1]))
  return(df_total)
}

comp_mean <- function(pop_name){
  random_mean <- c()
  for(i in seq(1,70)){
    df <- format(paste0("./15_LD_outlier_snps/", pop_name, ".ngsld.output"), paste0("./15_LD_random_snps/", pop_name ,"_", i, ".ngsld.output"))
    df2 <- summarySE(df, measurevar="V7", groupvars=c("group"))                            
    message(i)
    print(df2$V7[2])
    random_mean = c(random_mean, df2$V7[2])
  }
  
  df<-format(paste0("./15_LD_outlier_snps/", pop_name, ".ngsld.output"), paste0("./15_LD_random_snps/", pop_name ,"_1.ngsld.output"))
  observed_mean  <- summarySE(df, measurevar="V7", groupvars=c("group"))$V7[1]
  random_mean <- as.data.frame(random_mean)
  random_mean_point <- mean(random_mean$random_mean)
  message("observed_mean: ", observed_mean)
  message("random_mean_point: ", random_mean_point)
  
  p <- sum(random_mean > observed_mean)/dim(random_mean)[1]
  message("p-value: ", p)
  hist_p <- ggplot(random_mean, aes(x=random_mean, fill=NULL)) +
    geom_histogram(fill="white", color="grey", bins = 15)+
    geom_vline(aes(xintercept=mean(random_mean)), color="grey",
               linetype="dashed")+
    labs(x=expression(paste(Mean~Linkage~disequilibrium~(r^2))), y = "Count")+
    theme_classic() +
    theme(text = element_text(size = 15)) +
    geom_text(x=observed_mean-0.0022, label=paste("SGS candidate mean \n", "=", round(observed_mean, digits = 3)), y=15, color="red", size = 4)+
    geom_text(x=random_mean_point+0.0027, label=paste("Mean of null distribution =\n", "=", round(random_mean_point, digits = 3)), y=18, color="Black", size = 4)+
    geom_vline(xintercept=observed_mean, linetype="dashed", color = "red") 
  
  print(hist_p)
  graph2ppt(file=paste0(pop_name, "_r2_comp"), width=4, height=3) 
  return(random_mean)
  
}

comp_mean("NB_18")
comp_mean("NB_19")
comp_mean("NB_21")
comp_mean("Ref_19")
comp_mean("HC_19")
comp_mean("HC_18")
comp_mean("HC_21")
comp_mean("Sur_19")
comp_mean("Sur_20")

#####################################
# potentail uility for zoom-in test #
#####################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/zoom_in")
random_mean <- c()
for(i in seq(1,100)){
  df <- format("./outlier_output/CHR19.chr5.ngsld.output", paste0("./random_output/CHR19.", i, ".ngsld.output"))
  df2 <- summarySE(df, measurevar="V7", groupvars=c("group"))
  message(i)
  print(df2$V7[2])
  random_mean = c(random_mean, df2$V7[2])
}

df<-format("./outlier_output/CHR19.chr5.ngsld.output", "./random_output/CHR19.1.ngsld.output")
observed_mean  <- summarySE(df, measurevar="V7", groupvars=c("group"))$V7[1]
random_mean <- as.data.frame(random_mean)
random_mean_point <- mean(random_mean$random_mean)

sum(random_mean$random_mean > observed_mean)

hist_p <- ggplot(random_mean, aes(x=random_mean, fill=NULL)) +
  geom_histogram(fill="white", color="grey", bins = 15)+
  geom_vline(aes(xintercept=mean(random_mean)), color="grey",
             linetype="dashed")+
  labs(x=expression(paste(Mean~Linkage~disequilibrium~(r^2))), y = "Count")+
  theme_classic() +
  theme(text = element_text(size = 15)) +
  geom_text(x=observed_mean-0.0022, label=paste("SGS candidate mean \n", "=", round(observed_mean, digits = 3)), y=15, color="red", size = 4)+
  geom_text(x=random_mean_point+0.0027, label=paste("Mean of null distribution =\n", "=", round(random_mean_point, digits = 3)), y=18, color="Black", size = 4)+
  geom_vline(xintercept=observed_mean, linetype="dashed", color = "red") 

hist_p

graph2ppt(file="CHR19_REF19_r2_comp", width=4, height=3) 
