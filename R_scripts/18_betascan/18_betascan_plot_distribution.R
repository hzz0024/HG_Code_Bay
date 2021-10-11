###############################################
##########  delta_p distribution plot  ########
###############################################
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(grid)
library(gtable)
library(export)
setwd("~/Documents/HG/DelBay19_adult/18_betascan/betascore/extract_betascore/REF19")
pname = "./REF19_outlier.txt"
pname = "./NB_outlier.txt"
pname = "./HC_outlier.txt"
pname = "./CHR19_outlier.txt"
dat1 = read.delim(pname, header = TRUE, sep='\t')
dat1$group = "SGS outliers"

pname_all = "REF19_all.txt"
pname_all = "./NB_all.txt"
pname_all = "HC_all.txt"
pname_all = "CHR19_all.txt"
dat2 = read.delim(pname_all, header = TRUE, sep='\t')
dat2$group = "genome-wide"

dat_plot = rbind(dat1, dat2)

p2 <- ggplot(dat_plot, aes(x=Beta1, fill=group)) +  
  geom_density(alpha = 0.4) + theme(legend.position = "right") + scale_fill_manual(values=c("#F5F5F5", "#666666")) + 
  ggtitle(expression(HC~beta*~score~distribution~between~genome-wide~and~SGS~outliers)) +
  xlab(expression(beta*~score)) + ylab("Density")
p2

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
tiff("11_SGS_indep_p_zcomb_FDR005_1.tiff", units="in", width=8, height=8, res=300)
comb_plot<-grid.arrange(p2,empty,p1, p3, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()
gtable_show_layout(comb_plot)

graph2ppt(file="comb_plot.pptx", width=10, height=10)
annotate_figure(figure,
                top = text_grob("Visualizing mpg", color = "red", face = "bold", size = 14),
                bottom = text_grob("Data source: \n mtcars data set", color = "blue",
                                   hjust = 1, x = 1, face = "italic", size = 10),
                left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                right = "I'm done, thanks :-)!",
                fig.lab = "Figure 1", fig.lab.face = "bold"
)