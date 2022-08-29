#install.packages("remotes")
#remotes::install_github("JimWhiting91/genotype_plot")
library(devtools)
#install('/Users/ryan/Downloads/genotype_plot-master')
#install.packages("glue")
library(glue)
library(vcfR)
#library(GenotypePlot)
library(cowplot)
library(rlang)
library(export)

library(data.table)
library(ggdendro)
library(ggplot2)
library(vcfR)
library(tidyr)
library(dplyr)

setwd("~/Dropbox/Mac/Documents/Backup/Ryan_workplace/CVreseq_genotype_plot")
install('/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/genotype_plot-master')
library(GenotypePlot)
install.packages("extrafont")
library(extrafont)
font_import()
loadfonts(device="win")       #Register fonts for Windows bitmap output
fonts()   

setwd("~/Dropbox/Mac/Documents/Backup/Ryan_workplace/CVreseq_genotype_plot")
#my_vcf <- read.vcfR("CS_HC-HCVA_CLP.outlier.vcf.gz")
pop_name <- "HC_CS-CLP-HCVA.txt"
popmap = read.delim(pop_name, header = TRUE, sep='\t')
new_plot <- genotype_plot(vcf    =  "CS_HC-HCVA_CLP.outlier.vcf.gz",
                          chr    = 2,
                          start  = 33597600,  
                          end    = 33629600,
                          popmap = popmap,                              
                          cluster        = F,                           
                          snp_label_size = 10000,                          
                          colour_scheme=c("#FADCD9","#9986A5","#79402E"))   
cowplot::plot_grid(new_plot$positions, new_plot$genotypes, axis="tblr",
                   align="v", nrow=2, ncol=1, rel_heights=c(1,9))

genos = new_plot$genos
# Add dendrogram tips with points
pop_name <- "HC_CS-CLP-HCVA.txt"
popmap = read.delim(pop_name, header = TRUE, sep='\t')
Population = popmap$pop

# tmp = new_plot$dendro_labels
# group = strsplit(tmp, "_\\s*(?=[^_]+$)", perl=TRUE)
# #group = strsplit(new_plot$dendro_labels, "_\\s*(?=[^_]+$)", perl=TRUE)
# group = sapply(group, "[[", 1)
Origin = c(rep('Low salinity', 24))
#Origin[group %in% Population[1:12]] = 'Del'
Origin[7:12] = 'High salinity'
Origin[19:24] = 'High salinity'
#tmp = data.frame(group, Origin)
#Population = group
df = data.frame(y=24:1, x=rep(c(min(genos$snp)-800, min(genos$snp)-1400, min(genos$snp)-2000, min(genos$snp)-2000, min(genos$snp)-1400, min(genos$snp)-800),4))
Population <- factor(Population, levels = c("HC",  "CS", "CLP", "HC-VA"))
df$Origin = Origin
df$Population = Population
#df = data.frame(x=1:length(new_plot$dendro_labels), y=rep(-1,24))
dendro_with_tips <-  new_plot$genotypes +
  geom_point(data=df,aes(x=min(ggplot_build(new_plot$genotypes)$data[[1]]$x-1000),y=y,color=Population,shape=Origin),size=6, alpha=0.95)  +
  scale_shape_manual(values=c(18, 20)) + 
  scale_color_manual(values = c("#E31A1c", "#FB9A99", "#33A02C","#B2DF8A"),  name="Population/Line") +
  theme(legend.position="right") +
  theme(legend.key = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0.6, "cm"),
        legend.key.width = unit(0.8,"cm"))+
  guides(color = guide_legend(override.aes = list(size=6, shape = c(20,18,20,18),alpha = 0.95), order = 2), 
         alpha=FALSE,size=FALSE,shape= guide_legend(override.aes = list(size=4, shape = c(23,21)), order = 1), 
         fill = guide_legend(override.aes = list(size=6))) +
  # 1 CS #FB9A99 2 CLP #33A02C 3 HCVA #FB9A99 4 HC #E31A1c
  theme(legend.title = element_text(size = 12),
        legend.text=element_text(size=12)) + 
  theme(text=element_text(family="Times New Roman", face="bold", size=12, colour="black"))
  


dendro_with_tips
ggsave(filename="Final_plot.png", plot=dendro_with_tips, device="png",
       height=12, width=8, units="in", dpi=500)

graph2ppt(file="Final_genotype_plot",width=6,height=10)

