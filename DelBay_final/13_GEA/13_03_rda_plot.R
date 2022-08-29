####################################
##########  Delta_p plot ###########
####################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/03_rda")

library(ggplot2)
library(stringr)
###################### load maf values for delta_p ##########
# need to add the SNP_ID header to the loading file 
file1 = 'info_pop_env_rda.txt_rda_loading.txt' # loadings for SNPs
dat1 = read.delim(file1, header = T, sep=' ')
head(dat1)
new = str_split_fixed(dat1$SNP_ID, "_", 2)
dat1$chromo = new[,1]
dat1$position = new[,2]
dat1$chromo <- factor(dat1$chromo, levels = c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10'))
dat1_sort <- dat1[with(dat1, order(chromo, position)),]
# change the columns into numbers
#i <- c(2,6)
#dat1_sort[ , i] <- apply(dat1_sort[ , i], 2,function(x) as.numeric(as.character(x)))   # Specify own function within apply
# covert the position to 10M/unit
dat1_sort$RDA1 <- as.numeric(dat1_sort$RDA1)
dat1_sort$position <- as.numeric(dat1_sort$position)
dat1_sort$position = dat1_sort$position/1e+7
#dat1_sort = dat1_sort[1:10000,]
lims <- mean(dat1_sort$RDA1) + c(-1, 1) * 3 * sd(dat1_sort$RDA1)     # find loadings +/-z sd from mean loading
lims[1]
lims[2]

rda1 <- dat1_sort$RDA1
outlier <- cbind.data.frame(dat1_sort$chromo[rda1>lims[2]], dat1_sort$position[rda1>lims[2]], dat1_sort$RDA1[rda1>lims[2]], dat1_sort$RDA2[rda1>lims[2]])
colnames(outlier) <- c("chr","pos","RDA1","RDA2")
jpeg("RDA_Manhattan.jpg", width = 16, height = 6, units = 'in', res = 300)
ggplot(dat1_sort, aes(x=position, y=RDA1))+ 
  geom_point(aes(colour = cut(RDA1, c(-Inf, -0.0413278, 0.0398489, Inf))),size = 1, show.legend = F)+
  scale_color_manual(name = "RDA1",
                     values = c("(-Inf,-0.0413278]" = "grey",
                                "(-0.0413278,0.0398489]" = "orange",
                                "(0.0398489, Inf]" = "red"))+
  theme_classic()+
  facet_grid(cols = vars(chromo), scales = "free_x", space="free_x") +
  #theme(strip.text.x = element_blank())+
  theme(text = element_text(size=15)) +
  theme(axis.text.x = element_text(color = "grey20", size = 10)) +
  geom_hline(aes(yintercept =0.0398489), linetype="dotted", size=0.8, col="red", show.legend = FALSE)+
  geom_hline(aes(yintercept =-0.0413278), linetype="dotted", size=0.8, col="orange", show.legend = FALSE) +
  xlab("Position (10M bp)") +
  ylab("SNP loadings on RDA1")

dev.off()
