###############################################################
############### Plot the allele frequency trend ###############
###############################################################
setwd("~/Documents/HG/DelBay19_adult/13_env_gen_association/plot")
# load the outliers from SGS with same directionality
fname1 = 'SGS_outlier_same_dir.txt'  
outlier = read.delim(fname1, header = FALSE, sep='\t')
#outlier$id = paste0(outlier$chromo,'_',outlier$position)

# specify the number of decimal point
options("scipen"=100, "digits"=6)

# load the snp info and allele frequency data from each population 
fname = "by_pop_0.05_pctind0.7_maxdepth3.mafs"  
df_maf = read.delim(fname, header = TRUE, sep=' ')
head(df_maf)
dim(df_maf)
df_maf$id = paste0(df_maf$chromo,'_',df_maf$position)
df_maf <- df_maf[with(df_maf, order(id)),]
head(df_maf)
# # load the delta_p information from DelBay19 challenge
fname = "ps_Del19_challenge.txt"
df_deltap = read.delim(fname, header = FALSE, sep='\t')
dim(df_deltap)
colnames(df_deltap) <- c("chromo", "position", "p1", "p0", "delta_p", "ps", "outlier")
df_deltap$id = paste0(df_deltap$chromo,'_',df_deltap$position)
df_deltap <- df_deltap[with(df_deltap, order(id)),]
head(df_deltap)
sum(df_maf$id == df_deltap$id)
# only extract GEA outliers

data1 = cbind(df_maf$chromo[df_maf$id %in% outlier$V1], 
              df_maf$position[df_maf$id %in% outlier$V1], 
              df_maf$major[df_maf$id %in% outlier$V1], 
              df_maf$minor[df_maf$id %in% outlier$V1],
              df_maf$id[df_maf$id %in% outlier$V1],
              df_maf$freqHC[df_maf$id %in% outlier$V1],
              df_maf$freqARN[df_maf$id %in% outlier$V1],
              df_maf$freqCOH[df_maf$id %in% outlier$V1],
              df_maf$freqSR[df_maf$id %in% outlier$V1],
              df_maf$freqNB[df_maf$id %in% outlier$V1],
              df_maf$freqHC[df_maf$id %in% outlier$V1]-df_maf$freqNB[df_maf$id %in% outlier$V1])

data2 =  cbind(df_deltap$id[df_deltap$id %in% outlier$V1],
               df_deltap$p1[df_deltap$id %in% outlier$V1],
               df_deltap$p0[df_deltap$id %in% outlier$V1],
               df_deltap$delta_p[df_deltap$id %in% outlier$V1],
               df_deltap$ps[df_deltap$id %in% outlier$V1])

outlier_list <- as.data.frame(cbind(data1, data2))

colnames(outlier_list) <- c("chromo", "position", "major", "minor", "id", 
                            "freqHC", "freqARN", "freqCOH", "freqSR", "freqNB", "Wdelta_p",
                            "id2", "freqSurv", "freqRef", "Cdelta_p", "ps" )


i <- c(seq(6,11), seq(13,16))
outlier_list[ , i] <- apply(outlier_list[ , i], 2,            # Specify own function within apply
                            function(x) as.numeric(as.character(x)))
sapply(outlier_list, class)
#write.table(outlier_list,"SGS_outlier.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)

# count how many SNPs with positive or negative delta_p change
pos_challenge <- length(outlier_list$id[which(outlier_list$Cdelta_p>0)])
neg_challenge <- length(outlier_list$id[which(outlier_list$Cdelta_p<0)])
message("outlier has ", pos_challenge, " SNP with positive delta_p (Surv - Ref)")
message("outlier has ", neg_challenge, " SNP with negative delta_p (Surv - Ref)")

pos_wild <- length(outlier_list$id[which(outlier_list$Wdelta_p>0)])
neg_wild <- length(outlier_list$id[which(outlier_list$Wdelta_p<0)])
message("outlier has ", pos_wild, " SNP with positive delta_p (HC - NB)")
message("outlier has ", neg_wild, " SNP with negative delta_p (HC - NB)")

pos_wild <- length(outlier_list$id[which((outlier_list$freqHC-outlier_list$freqSR)>0)])
neg_wild <- length(outlier_list$id[which((outlier_list$freqHC-outlier_list$freqSR)<0)])
message("outlier has ", pos_wild, " SNP with positive delta_p (HC - SR)")
message("outlier has ", neg_wild, " SNP with negative delta_p (HC - SR)")

# extract the snps with same directionality in allele frequency change 
df_pos_direction <- outlier_list[which(outlier_list$Wdelta_p>0 & outlier_list$Cdelta_p > 0), ]
df_neg_direction <- outlier_list[which(outlier_list$Wdelta_p<0 & outlier_list$Cdelta_p < 0), ]
message("SNPs share the same pos direction: ", length(df_pos_direction$chromo))
message("SNPs share the same neg direction: ", length(df_neg_direction$chromo))
df_same_direction <- rbind(df_pos_direction, df_neg_direction)
message("SNPs share the same direction: ", length(df_same_direction$chromo))

#write.table(df_same_direction,"SGS_outlier_same_direction.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
# show the distribution of delta_ps
hist(df_same_direction$Wdelta_p)
hist(df_same_direction$Cdelta_p)

# polarize the "minor" alleles
# for (i in 1:length(outlier_list$position)){
#   if (outlier_list$M_Pearson[i] <0){
#     outlier_list$tmp[i] = outlier_list$major[i]
#     outlier_list$major[i] = outlier_list$minor[i]
#     outlier_list$minor[i] = outlier_list$tmp[i]
#     outlier_list$freqHC[i] = 1-outlier_list$freqHC[i]
#     outlier_list$freqARN[i] = 1-outlier_list$freqARN[i]
#     outlier_list$freqCOH[i] = 1-outlier_list$freqCOH[i]
#     outlier_list$freqSR[i] = 1-outlier_list$freqSR[i]
#     outlier_list$freqNB[i] = 1-outlier_list$freqNB[i]
#   }
#   else {
#     outlier_list$tmp[i] = outlier_list$major[i]
#   }
# }
#write.table(outlier_list,"GEA_outlier_226_corrected_by_pearson.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)

# plot the allele frequency changes across the clines 
new_outlier = data.frame(id=rep(seq(length(outlier_list[,1])), each = 5))
all_dat_ = c()
for(i in seq(length(outlier_list[,1]))){
  dat_ = c(outlier_list[i,]$freqHC, outlier_list[i,]$freqARN, outlier_list[i,]$freqCOH, outlier_list[i,]$freqSR, outlier_list[i,]$freqNB)
  all_dat_ = c(all_dat_, dat_)
}
new_outlier$index = rep(c(10,20,30,40,50), length(outlier_list[,1]))
#new_outlier$index = rep(c("HC","ARN","COH","SR","NB"), length(outlier_list[,1]))
new_outlier$af = all_dat_

dat = new_outlier

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="plain", size=9),
          axis.text=element_text(face="plain",size=18),axis.title = element_text(face="plain",size=20),plot.title = element_text(face = "bold", hjust = 0.5,size=20))
)

#dat=read.table("outlier_for_lm_plot.txt", header=T)
ggplot(dat) + 
  mytheme+
  geom_smooth(method = lm, aes(x=index, y=af, group=id), colour="grey", size=0.01, se=FALSE) + 
  geom_smooth(method = lm, aes(x=index, y=af), colour="red", fill='red') +
  ylab("Allele frequency") +
  scale_x_discrete(name = 'Population', limits = c(10,20,30,40,50), labels=c("HC","ARN","COH","SR","NB"))

ggplot(dat) + 
  mytheme+
  #geom_point(aes(x=index, y=af, group=id), colour="grey", size=0.01) + 
  geom_line(aes(x=index, y=af, group=id), color="grey", size = 0.05) +
  stat_summary(aes(x=index, y = af,group = 1), fun = mean, geom="line", colour="red", size = 0.5)+
  ylab("Allele frequency") +
  scale_x_discrete(name = 'Population', limits = c(10,20,30,40,50), labels=c("HC","ARN","COH","SR","NB"))


###############################################################
###################### Bayes factor plot ######################
###############################################################

library(ggplot2)
#setwd("~/Documents/Ryan_workplace/DelBay_adult/13_env_gen_association/plot_outlier_trend/")
#setwd("/Volumes/cornell/Cohort_adaptation/DelBay_adult/GE_association")
#load bf values
#bf_allsnps<-read.table("allsnps.controlled.env.output_summary_betai_reg.out", header = T)
setwd("~/Documents/HG/DelBay19_adult/13_env_gen_association/salinity1")
setwd("~/Documents/HG/DelBay19_adult/13_env_gen_association/salinity2")
setwd("~/Documents/HG/DelBay19_adult/13_env_gen_association/salinity3")
bf_allsnps1<-read.table("allsnps.controlled.env.1.output_summary_betai_reg.out", header = T)
bf_allsnps2<-read.table("allsnps.controlled.env.2.output_summary_betai_reg.out", header = T)
bf_allsnps3<-read.table("allsnps.controlled.env.3.output_summary_betai_reg.out", header = T)
bf_allsnps4<-read.table("allsnps.controlled.env.4.output_summary_betai_reg.out", header = T)
bf_allsnps5<-read.table("allsnps.controlled.env.5.output_summary_betai_reg.out", header = T)
bf_allsnps6<-read.table("allsnps.controlled.env.6.output_summary_betai_reg.out", header = T)
bf_allsnps7<-read.table("allsnps.controlled.env.7.output_summary_betai_reg.out", header = T)
bf_allsnps8<-read.table("allsnps.controlled.env.8.output_summary_betai_reg.out", header = T)
bf_allsnps9<-read.table("allsnps.controlled.env.9.output_summary_betai_reg.out", header = T)
bf_allsnps10<-read.table("allsnps.controlled.env.10.output_summary_betai_reg.out", header = T)
head(bf_allsnps1)
#colnames(bf_allsnps) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
colnames(bf_allsnps1) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
colnames(bf_allsnps2) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
colnames(bf_allsnps3) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
colnames(bf_allsnps4) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
colnames(bf_allsnps5) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
colnames(bf_allsnps6) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
colnames(bf_allsnps7) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
colnames(bf_allsnps8) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
colnames(bf_allsnps9) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
colnames(bf_allsnps10) <- c("COVARIABLE","MRK","M_Pearson","SD_Pearson", "Bayes_Factor", "Beta_is","SD_Beta_is","eBPis")
#head(bf_allsnps)
head(bf_allsnps1)
cor(bf_allsnps1$Bayes_Factor,bf_allsnps3$Bayes_Factor) 
cor(bf_allsnps1$Bayes_Factor,bf_allsnps2$Bayes_Factor)
#load position info about the SNPs.
SNP_pos<-read.table("by_pop_0.05_pctind0.7_maxdepth3.snps", header=T)
SNP_pos$position = SNP_pos$position/1e+7
#should be same nb of rows
dim(bf_allsnps1)
dim(SNP_pos)

#tmpdata = cbind(bf_allsnps1$Bayes_Factor, bf_allsnps2$Bayes_Factor,bf_allsnps3$Bayes_Factor,bf_allsnps4$Bayes_Factor,bf_allsnps5$Bayes_Factor
#                ,bf_allsnps6$Bayes_Factor,bf_allsnps7$Bayes_Factor)

tmpdata = cbind(bf_allsnps1$Bayes_Factor, bf_allsnps2$Bayes_Factor,bf_allsnps3$Bayes_Factor,bf_allsnps4$Bayes_Factor,bf_allsnps5$Bayes_Factor
                ,bf_allsnps6$Bayes_Factor,bf_allsnps7$Bayes_Factor,bf_allsnps8$Bayes_Factor,bf_allsnps9$Bayes_Factor,bf_allsnps10$Bayes_Factor)
BF_median <- apply(tmpdata, 1, median, na.rm=T)

bf_pos<-cbind(SNP_pos, BF_median)
bf_pos$chromo <- factor(bf_pos$chromo, levels = c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10'))
bf_pos_sort <- bf_pos[with(bf_pos, order(chromo, position)),]

jpeg("BF_Manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
ggplot(bf_pos_sort, aes(x=position, y=BF_median))+ 
  geom_point(aes(colour = cut(BF_median, c(-Inf, 10, 20, Inf))),size = 0.8,show.legend = F)+
  scale_color_manual(name = "BF_median",
                     values = c("(-Inf,10]" = "grey",
                                "(10,20]" = "orange",
                                "(20, Inf]" = "red"))+
  theme_classic()+
  facet_grid(cols = vars(chromo), scales = "free_x", space="free_x") +
  theme(text = element_text(size=20)) +
  theme(axis.text.x = element_text(color = "grey20", size = 10)) +
  geom_hline(aes(yintercept =20), linetype="dotted", size=1, col="red", show.legend = FALSE)+
  geom_hline(aes(yintercept =10), linetype="dotted", size=1, show.legend = FALSE) +
  #scale_x_continuous(labels = scales::scientific) +
  xlab("Position (10M bp)") +
  ylab("Bayes Factor")
dev.off()


###############################################################
###################### simulation plot ########################
###############################################################

#load bf values from simulatd data
# bf_simu<-read.table("simul.controlled.env.output_summary_betai_reg.out", header=T)
# head(bf_simu)
# 
# #calculate the threshold
# threshold_fdr0.01 = quantile(bf_simu$BF.dB,probs=0.99)
# threshold_fdr0.05 = quantile(bf_simu$BF.dB,probs=0.95)
# 
# #add it on the plot
# ggplot(bf_pos, aes(x=position, y=BF.dB., colour=chromosome))+ 
#   geom_point()+
#   theme_classic()+
#   facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x")+
#   geom_hline(aes(yintercept =threshold_fdr0.05), linetype="dotted", size=1, col="red", show.legend = FALSE)+
#   geom_hline(aes(yintercept =threshold_fdr0.01), linetype="dotted", size=1, show.legend = FALSE)
# 
# #output outliers
# bf_pos[bf_pos$BF.dB.>= threshold_fdr0.05,]

###############################################################
###################### Extract outlier ########################
###############################################################

#load bf values
#bf_allsnps<-read.table("allsnps.controlled.env.output_summary_betai_reg.out", header = T)
SNP_pos<-read.table("by_pop_0.05_pctind0.7_maxdepth3.snps", header=T)
bf_pos<-cbind(SNP_pos, bf_pos_sort$BF_median)

#load bf values from simulatd data
# bf_simu<-read.table("simul.controlled.env.output_summary_betai_reg.out", header=T)

#calculate the threshold from simu (or you cna use BF = 10)
# threshold_fdr0.01 = quantile(bf_simu$BF.dB,probs=0.99)
# threshold_fdr0.05 = quantile(bf_simu$BF.dB,probs=0.95)
# 
# outlier<-bf_pos[bf_pos$BF.dB.>=threshold_fdr0.05,]
#outlier<-bf_pos[bf_pos$BF.dB.>=10,]
#outlier<-bf_pos[bf_pos$`bf_pos_sort$BF_median`>= 20,]
outlier<-bf_pos[bf_pos$`bf_pos_sort$BF_median`>= 10,]
length(outlier$chromo)
outlier$id = paste0(outlier$chromo,'_',outlier$position) 
write.table(outlier, "./salinity_2032113_outlier_BF10.txt", row.names=F, quote=F, sep="\t")

bf_allsnps0<-read.table("allsnps.controlled.env.3.output_summary_betai_reg.out", header = T)
bf_pos0<-cbind(SNP_pos, bf_allsnps0)
outlier0<-bf_pos0[bf_pos0$BF.dB.>=20,]
length(outlier0$chromo)
outlier0$id = paste0(outlier0$chromo,'_',outlier0$position) 
write.table(outlier0, "./salinity_2032113_outlier_BF20_0.txt", row.names=F, quote=F, sep="\t")

bf_allsnps1<-read.table("allsnps.controlled.env.1.output_summary_betai_reg.out", header = T)
bf_pos1<-cbind(SNP_pos, bf_allsnps1)
outlier1<-bf_pos1[bf_pos1$BF.dB.>=20,]
length(outlier1$chromo)
outlier1$id = paste0(outlier1$chromo,'_',outlier1$position) 
write.table(outlier1, "./salinity_2032113_outlier_BF20_1.txt", row.names=F, quote=F, sep="\t")

bf_allsnps2<-read.table("allsnps.controlled.env.2.output_summary_betai_reg.out", header = T)
bf_pos2<-cbind(SNP_pos, bf_allsnps2)
outlier2<-bf_pos2[bf_pos2$BF.dB.>=20,]
length(outlier2$chromo)
outlier2$id = paste0(outlier2$chromo,'_',outlier2$position) 
write.table(outlier2, "./salinity_2032113_outlier_BF20_2.txt", row.names=F, quote=F, sep="\t")

s2 <- intersect(outlier0$id, (intersect(outlier1$id, outlier2$id)))
length(s2)

s1 <- intersect(outlier1$id, outlier2$id)
length(s1)


###############################################################
###############         2d density plot          ###############
###############################################################

highlight1 = read.delim("GEA_BF_20_Del19_FDR_2K.intersect", header = FALSE, sep='\t')
highlight_plot1 = outlier_list[which(outlier_list$id %in% highlight1$V4),] # note the snp count is less than GEA_BF_20_Del19_FDR_10K.intersect because of duplicate intersects

SGS_Del19_FDR<-read.table("Del19_FDR_outlier.list", header = T)
GEA_BF20<-read.table("salinity_2032113_outlier_BF20.txt", header = T)
GEA_BF20$id = paste0(GEA_BF20$chromo,'_',GEA_BF20$position)
SGS_Del19_FDR$id[which(SGS_Del19_FDR$id %in% GEA_BF20$id)]

highlight2 = read.delim("GEA_BF_20_Del19_FDR_2K.intersect", header = FALSE, sep='\t')
highlight_plot2 = outlier_list[which(outlier_list$id %in% "NC_035786.1_8772371"),] # note the snp count is less than GEA_BF_20_Del19_FDR_10K.intersect because of duplicate intersects


library(tidyverse)
# https://www.r-graph-gallery.com/2d-density-plot-with-ggplot2.html#hex
# Area + contour
tiff("Environment_association_outlier_delta_p_plain.tiff", units="in", width=8, height=8, res=300)
ggplot(outlier_list, aes(x=Cdelta_p, y=Wdelta_p) ) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed", color = "grey", size=0.5) +
  geom_density_2d(colour="grey") +
  theme(legend.position='right')+
  geom_point(color="red")+
  #geom_point(data=highlight_plot1, aes(x=Cdelta_p, y=Wdelta_p) , color='red', size=2) +
  #geom_point(data=highlight_plot2, aes(x=Cdelta_p, y=Wdelta_p) , color='blue', size=2) +
  scale_x_continuous(expression(Delta~italic(p)~"in"~challenge~"("~Surv~-~Ref~")"), limits = c(-0.5, 0.5)) + 
  scale_y_continuous(expression(Delta~italic(p)~"in"~wild~transect~"("~HC~-~NB~")"), limits = c(-0.5, 0.5))+ 
  theme(axis.text=element_text(size=20),
      text = element_text(size=20,family="Times"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_rect(colour = "black", size=0.5))
dev.off()

df_same_direction<-read.table("GEA_outlier_same_direction.txt", header = T)
# polarize the "minor" alleles
for (i in 1:length(df_same_direction$position)){
  if (df_same_direction$Wdelta_p[i] < 0){
    tmp = df_same_direction$major[i]
    df_same_direction$major[i] = df_same_direction$minor[i]
    df_same_direction$minor[i] = tmp
    df_same_direction$freqHC[i] = 1-df_same_direction$freqHC[i]
    df_same_direction$freqARN[i] = 1-df_same_direction$freqARN[i]
    df_same_direction$freqCOH[i] = 1-df_same_direction$freqCOH[i]
    df_same_direction$freqSR[i] = 1-df_same_direction$freqSR[i]
    df_same_direction$freqNB[i] = 1-df_same_direction$freqNB[i]
    df_same_direction$Wdelta_p[i] = -df_same_direction$Wdelta_p[i]
    df_same_direction$freqSurv[i] = 1-df_same_direction$freqSurv[i]
    df_same_direction$freqRef[i] = 1-df_same_direction$freqRef[i]
    df_same_direction$Cdelta_p[i] = -df_same_direction$Cdelta_p[i]
  }
}

write.table(df_same_direction,"GEA_outlier_same_direction_corrected.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)

tiff("Environment_association_outlier_same_direction_delta_p.tiff", units="in", width=8, height=8, res=300)
ggplot(df_same_direction, aes(x=Cdelta_p, y=Wdelta_p) ) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed", color = "grey", size=0.5) +
  geom_density_2d(colour="grey") +
  theme(legend.position='right')+
  geom_point(color="red")+
  scale_x_continuous(expression(Delta~italic(p)~"in"~challenge~"("~Surv~-~Ref~")"), limits = c(0, 0.5)) + 
  scale_y_continuous(expression(Delta~italic(p)~"in"~wild~transect~"("~HC~-~NB~")"), limits = c(0, 0.5))+ 
  theme(axis.text=element_text(size=20),
        text = element_text(size=20,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))
dev.off()
#######################
#  plot delta_p      #
#######################
# plot distribution

jpeg("EA_32_outlier_19_direction_trend.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
# rank order based NB frequency
df_same_direction_rank <- df_same_direction[with(df_same_direction, order(freqRef)),]
mycol <- t_col("grey", perc = 50, name = "lt.grey")

hist(df_same_direction_rank$Wdelta_p, main=NULL, breaks = 50,
     xlab = expression(Delta~italic(p)), ylab = "Frequency")
# plot the trend
cnt = length(df_same_direction_rank$chromo)
cnt

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  ## Save the color
  invisible(t.col)
}
order = seq(1,length(df_same_direction_rank$chromo),1)
df_same_direction_rank$order = order
plot(df_same_direction_rank$order, df_same_direction_rank$freqNB,col="red",pch=20,cex=1.5,ylim=c(0,1),xlab="Outliers",ylab="Allele frequency",xlim=c(0,cnt+1)) # ylab=expression(italic("p"))
for (l in 1:nrow(df_same_direction_rank))
{
  segments(l,df_same_direction_rank$freqHC[l],l,df_same_direction_rank$freqNB[l],col=mycol,lwd=1.2, lty = "dotted")
  if (df_same_direction_rank$Wdelta_p[l]<=0)
  {
    points(l,df_same_direction_rank$freqHC[l],col="blue",pch=6,cex=1.5)
  }
  if (df_same_direction_rank$Wdelta_p[l]>0)
  {
    points(l,df_same_direction_rank$freqHC[l],col="blue",pch=17,cex=1.5)
  }
}

points(df_same_direction_rank$order+0.24, df_same_direction_rank$freqRef,col="red",pch=20,cex=1.5,ylim=c(0,1),ylab=expression(italic("p")),xlab="SNPs",xlim=c(0,cnt+1))
for (l in 1:nrow(df_same_direction_rank))
{
  segments(l+0.24,df_same_direction_rank$freqSurv[l],l+0.24,df_same_direction_rank$freqRef[l],col=mycol,lwd=1.2, lty = "dotted") #@ add by HG
  if (df_same_direction_rank$Wdelta_p[l]<=0)
  {
    points(l+0.24,df_same_direction_rank$freqSurv[l],col="green",pch=6,cex=1.5)
  }
  if (df_same_direction_rank$Wdelta_p[l]>0)
  {
    points(l+0.24,df_same_direction_rank$freqSurv[l],col="green",pch=17,cex=1.5) #@ add by HG
  }
}
dev.off()

###########################
## test on random shares ##
###########################
daf<-read.table("ps_Del19_challenge.txt", header = F)
head(daf)
colnames(daf) <- c("chr","pos","freqSurv","freqRef", "delta_p", "ps","outlier")
head(daf)
daf$id = paste0(daf$chr,'_',daf$pos)
snp_id = as.data.frame(daf$id)
library(dplyr)
n_bootstraps=1000
boot_cnt = rep(0, n_bootstraps)
a = seq(1,2032113,1)
for (i in 1:n_bootstraps){
  o1 = sample_n(snp_id, 2979)
  o2 = sample_n(snp_id, 33)
  cnt = length(intersect(o1$`daf$id`, o2$`daf$id`))
  boot_cnt[i] = cnt
}
boot_cnt
p = sum(boot_cnt >= 1)/length(boot_cnt)
p
hist(boot_cnt)

# test if overlaps exceed the expection

bed_format <- function(input, output){
  outlier_temp_SGS<- input
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  window<-2000
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$pos-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$pos+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,9,10,8)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

daf<-read.table("ps_Del19_challenge.txt", header = F)
colnames(daf) <- c("chr","pos","freqSurv","freqRef", "delta_p", "ps","outlier")
head(daf)
daf$id = paste0(daf$chr,'_',daf$pos)

n_bootstraps=1000
boot_cnt = rep(0, n_bootstraps)
#a = seq(1,2032113,1)
for (i in 1:n_bootstraps){
  idx1 = sort(sample.int(dim(daf)[1], 2979))
  idx2 = sort(sample.int(dim(daf)[1], 33))
  o1 = daf[idx1,]
  o2 = daf[idx2,]
  #o1 = sample_n(daf, 2979)
  #o2 = sample_n(daf, 226)
  bed_format(o1, "test1.bed")
  bed_format(o2, "test2.bed")
  #system("bedtools intersect -a test1.bed -b test2.bed > test.intersect", intern = TRUE)
  #cnt = system("cat test.intersect | wc -l", intern = TRUE)
  system("bedtools intersect -a test1.bed -b test2.bed -wo > test.intersect", intern = TRUE)
  if (file.size("test.intersect") == 0) next
  temp <- read.delim("test.intersect", header = FALSE, sep='\t')
  cnt <- length(unique(c(temp$V4,temp$V8)))
  print(cnt)
  boot_cnt[i] = cnt
}

boot_cnt <- as.numeric(as.character(boot_cnt))
hist(boot_cnt)
abline(v = mean(boot_cnt), col = "blue", lwd = 2)
abline(v = 13, col = "red", lwd = 2, lty=2)
# This test confirms whether the normal distribution of the data is violated.
shapiro.test(boot_cnt)
p = sum(boot_cnt >= 13)/length(boot_cnt)
p
t.test(boot_cnt,conf.level=0.95)


###########################
##### find intersects #####
###########################
setwd("~/Documents/Ryan_workplace/DelBay_adult/22_find_intersect")
setwd("/Volumes/cornell/Cohort_adaptation/DelBay_adult/GE_association")
bed_format <- function(input, output){
  outlier_temp_SGS<- read.delim(input, header = TRUE, sep='\t')
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  window<-10000
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$pos-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$pos+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,9,10,8)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format("Del19_FDR_outlier.list", "Del19_FDR_2K.bed")
bed_format("Del20_FDR_outlier.list", "Del20_FDR_2K.bed")
bed_format("NB_HC_FDR_outlier.list", "NB_HC_FDR_2K.bed")
bed_format("SR_HC_FDR_outlier.list", "SR_HC_FDR_2K.bed")

bed_format1 <- function(input, output){
  outlier_temp_SGS<- read.delim(input, header = TRUE, sep='\t')
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  window<-10000
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$position-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$position+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,17,18,5)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format1("GEA_BF20_outlier.txt", "GEA_BF20_10K.bed")

system("bedtools intersect -a GEA_BF20_10K.bed -b Del19_FDR_10K.bed -wo > GEA_BF_20_Del19_FDR_10K.intersect", intern = TRUE)
cnt = system("cat GEA_BF_20_Del19_FDR_10K.intersect | wc -l", intern = TRUE)
temp <- read.delim("GEA_BF_20_Del19_FDR_10K.intersect", header = FALSE, sep='\t')
cnt <- length(unique(c(temp$V4,temp$V8)))
cnt
