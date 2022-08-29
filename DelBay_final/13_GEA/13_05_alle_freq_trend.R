##################
## CD5 outliers ##
##################
library(wesanderson)
rm(list=ls())

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/07_trend_plot/")
snp_df <- read.delim("by_pop_0.05_pctind0.7_maxdepth3.mafs", header = T, sep=" ")
snp_df$snps <- paste0(snp_df$chromo, "_", snp_df$position)

outlier_df <- read.delim("CD5_shared_lfmm_rda.bed", header = F, sep="\t")
outliers <- paste0(outlier_df$V1, "_", outlier_df$V2)

snp_outliers_df <- snp_df[which(snp_df$snps %in% outliers),]

extract_df <- data.frame(snps = snp_outliers_df$snps, snp_outliers_df$freqHC_18, snp_outliers_df$freqARN_18, snp_outliers_df$freqCOH_18, snp_outliers_df$freqSR_18, snp_outliers_df$freqNB_18,
                         snp_outliers_df$freqHC_19, snp_outliers_df$freqARN_19, snp_outliers_df$freqCOH_19, snp_outliers_df$freqSR_19, snp_outliers_df$freqNB_19,
                         snp_outliers_df$freqHC_21, snp_outliers_df$freqARN_21, snp_outliers_df$freqCOH_21, snp_outliers_df$freqSR_21, snp_outliers_df$freqBS_21, snp_outliers_df$freqBEN_21, snp_outliers_df$freqNAN_21, snp_outliers_df$freqNB_21)

colnames(extract_df) = c("SNP", "HC_18", "ARN_18", "COH_18", "SR_18", "NB_18",
                         "HC_19", "ARN_19", "COH_19", "SR_19", "NB_19",
                         "HC_21", "ARN_21", "COH_21", "SR_21", "BS_21", "BEN_21", "NAN_21", "NB_21")
plot_df <- reshape2::melt(extract_df, id.vars=c("SNP"),
                variable.name = "Interval", 
                value.name = "Values")
plot_df$year <- c(rep("Cohort_2018", dim(snp_outliers_df)[1]*5), rep("Cohort_2019", dim(snp_outliers_df)[1]*5), rep("Cohort_2021", dim(snp_outliers_df)[1]*8))
colnames(plot_df) = c("SNP", "Pop", "AF", "Year")
order <- c("HC_18", "ARN_18", "COH_18", "SR_18", "NB_18",
           "HC_19", "ARN_19", "COH_19", "SR_19", "NB_19",
           "HC_21", "ARN_21", "COH_21", "SR_21", "BS_21", "BEN_21", "NAN_21", "NB_21")
plot_df$Pop <-factor(plot_df$Pop, levels=order)
plot_df['row']=rep(seq(1,119),18)

mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="", strip.text = element_text(face="plain", size=12),
          axis.text=element_text(face="plain",size=12),axis.title = element_text(face="plain",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=12))
)

cbPalette <- wes_palette("FantasticFox1", 5, type = "continuous")

# create mean by group
mean <- plot_df%>% group_by(Pop)%>%summarise(mean_val=mean(AF))
mean$Year = c(rep("Cohort_2018", 5), rep("Cohort_2019", 5), rep("Cohort_2021", 8))

str(plot_df) 

p1 <- ggplot(plot_df, aes(x=Pop, y=AF, shape=Year, color=Year)) +
  geom_line(aes(group = row),color="grey",alpha=0.5)+
  geom_point(aes(fill=Year),size=2, alpha=0.8)+
  #geom_jitter(aes(fill=Year),size=2)+
  scale_shape_manual(values=c(15,16,17), name="Year") + 
  scale_color_manual(values=cbPalette , name="Year") +
  geom_point(data= mean, aes(x= Pop, y=mean_val), col="red", size=4)+
  geom_line(data= mean, aes(x= Pop, y=mean_val), col="red")+
  mytheme+
  ylab("Allele frequency")+
  theme(axis.title.x=element_blank(),
        text=element_text(family="Times New Roman", face="bold", size=8, colour="black"),
        axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))+
  facet_grid(~Year,scales = 'free_x', space = 'free_x', switch = 'x')
  
p1
  
graph2ppt(file="CD5_outlier_trend",width=10,height=6)

####################
## MAX10 outliers ##
####################

rm(list=ls())

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/07_trend_plot/")
snp_df <- read.delim("by_pop_0.05_pctind0.7_maxdepth3.mafs", header = T, sep=" ")
snp_df$snps <- paste0(snp_df$chromo, "_", snp_df$position)

outlier_df <- read.delim("MAX10_shared_lfmm_rda.bed", header = F, sep="\t")
outliers <- paste0(outlier_df$V1, "_", outlier_df$V2)

snp_outliers_df <- snp_df[which(snp_df$snps %in% outliers),]

extract_df <- data.frame(snps = snp_outliers_df$snps, snp_outliers_df$freqHC_18, snp_outliers_df$freqARN_18, snp_outliers_df$freqCOH_18, snp_outliers_df$freqSR_18, snp_outliers_df$freqNB_18,
                         snp_outliers_df$freqHC_19, snp_outliers_df$freqARN_19, snp_outliers_df$freqCOH_19, snp_outliers_df$freqSR_19, snp_outliers_df$freqNB_19,
                         snp_outliers_df$freqHC_21, snp_outliers_df$freqARN_21, snp_outliers_df$freqCOH_21, snp_outliers_df$freqSR_21, snp_outliers_df$freqBS_21, snp_outliers_df$freqBEN_21, snp_outliers_df$freqNAN_21, snp_outliers_df$freqNB_21)

colnames(extract_df) = c("SNP", "HC_18", "ARN_18", "COH_18", "SR_18", "NB_18",
                         "HC_19", "ARN_19", "COH_19", "SR_19", "NB_19",
                         "HC_21", "ARN_21", "COH_21", "SR_21", "BS_21", "BEN_21", "NAN_21", "NB_21")
plot_df <- reshape2::melt(extract_df, id.vars=c("SNP"),
                          variable.name = "Interval", 
                          value.name = "Values")
plot_df$year <- c(rep("Cohort_2018", dim(snp_outliers_df)[1]*5), rep("Cohort_2019", dim(snp_outliers_df)[1]*5), rep("Cohort_2021", dim(snp_outliers_df)[1]*8))
colnames(plot_df) = c("SNP", "Pop", "AF", "Year")
order <- c("HC_18", "ARN_18", "COH_18", "SR_18", "NB_18",
           "HC_19", "ARN_19", "COH_19", "SR_19", "NB_19",
           "HC_21", "ARN_21", "COH_21", "SR_21", "BS_21", "BEN_21", "NAN_21", "NB_21")
plot_df$Pop <-factor(plot_df$Pop, levels=order)
plot_df['row']=rep(seq(1,45),18)

mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="", strip.text = element_text(face="plain", size=12),
          axis.text=element_text(face="plain",size=12),axis.title = element_text(face="plain",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=12))
)

cbPalette <- wes_palette("FantasticFox1", 5, type = "continuous")

# create mean by group
mean <- plot_df%>% group_by(Pop)%>%summarise(mean_val=mean(AF))
mean$Year = c(rep("Cohort_2018", 5), rep("Cohort_2019", 5), rep("Cohort_2021", 8))

p1 <- ggplot(plot_df, aes(x=Pop, y=AF, shape=Year, color=Year)) +
  geom_line(aes(group = row),color="grey",alpha=0.5)+
  geom_point(aes(fill=Year),size=2, alpha=0.8)+
  #geom_jitter(aes(fill=Year),size=3)+
  scale_shape_manual(values=c(15,16,17), name="Year") + 
  scale_color_manual(values=cbPalette , name="Year") +
  geom_point(data= mean, aes(x= Pop, y=mean_val), col="red", size=4)+
  geom_line(data= mean, aes(x= Pop, y=mean_val), col="red")+
  mytheme+
  ylab("Allele frequency")+
  theme(axis.title.x=element_blank(),
        text=element_text(family="Times New Roman", face="bold", size=8, colour="black"),
        axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))+
  facet_grid(~Year,scales = 'free_x', space = 'free_x', switch = 'x')
p1

graph2ppt(file="MAX10_outlier_trend",width=10,height=6)


