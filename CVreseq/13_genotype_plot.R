install.packages("remotes")
#remotes::install_github("JimWhiting91/genotype_plot")
library(devtools)
install('/Users/ryan/Downloads/genotype_plot-master')
install.packages("glue")
library(glue)
library(vcfR)
library(ggplot2)
library(GenotypePlot)
library(cowplot)
library(rlang)
library(export)

setwd("~/Documents/Ryan_workplace/CVreseq_genotype_plot")
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
                          colour_scheme=c("#d4b9da","#e7298a","#980043"))   

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
  geom_point(aes(x=df$x, y=df$y, color=Population, shape=Origin), size = 5) +
  theme(legend.position="left") +
  theme(legend.key = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0.6, "cm"),
        legend.key.width = unit(0.8,"cm"))+
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



######## code from Jon
#Libraries and graphing parameters
library(data.table)
library(ggdendro)
library(wesanderson)
library(directlabels)
library(patchwork)
# I made a custom version of the genotype plot because it wasn't properyl reading the vcf files.  I think the newest version is fixed for this:
source("/home/jpuritz/src/genotype_plot/R/genotype_plot.R")
p_noinbred <- c("#FF7F00","#FDBF6F","#fda56f","#871a1a","#E31A1c","#FB9A99","#8a584a" ,"#33A02C","#B2DF8A","#86a395","#1F78B4","#A6CEE3","#A6CEE3")
# Code below is way more complicated than it needs to be.  I'll share these files.
pops <- read.table("population_coordinates.full.tab", header = TRUE)
popind <- read.table("popmap", header=TRUE)
popmap <-join(pops,popind,type = "inner")
popmap <- popmap[order(popmap$IND),]
popmap$Type <- factor(popmap$Type, levels=c("Inbred","Domesticated","Wild"))
popmap$POP <- factor(popmap$POP, levels=c("HI","SM","UMFS","NEH","NG", "HG","HC","CS","DEBY" ,"CLP","HC-VA","LOLA","CL","SL","OBOYS","LM"))
popmap.ni <- droplevels(subset(popmap, POP != "NG" & POP !="HG" & POP != "LM"))
gp.popmap <-setNames(data.frame(matrix(ncol = 2, nrow = length(popmap.ni$POP))), c("ind", "pop"))
gp.popmap$ind <- popmap.ni$IND
gp.popmap$pop <- popmap.ni$POP
# Function for  GP plot code ( No legend )
Chrom_Outlier_INV_GP_Plot <- function(chrom,ncbi) {
  plot_pal <- c(col_pal[as.integer(chrom)],wes_palette("IsleofDogs1"))
  outlier.inv.gp <- genotype_plot(vcf="inversion.outliers.vcf.gz",chr= ncbi,popmap = gp.popmap, start=1, end=94778955,cluster=TRUE,colour_scheme=plot_pal) 
  outlier.inv.gp$genotypes <- outlier.inv.gp$genotypes + theme(legend.text = element_text(size=6),legend.position = "none", legend.title=element_text(size=8)) + guides(fill = guide_legend(override.aes = list(size=2))) 
  new_df <- data.frame(matrix(ncol=1,nrow=78, dimnames=list(NULL, c("IND"))))
  new_df$IND <- outlier.inv.gp$dendro_labels
  new_df <- join(new_df, popmap.ni)
  dendro_with_tips.outlier <- outlier.inv.gp$dendrogram + 
    geom_jitter(data=new_df,aes(x=1:length(IND),y=-2.5,color=POP,fill=POP,shape=Type),size=3, alpha=0.95, height=2.5, width = 0)+
    scale_shape_manual(values=c(23,21),name="Origin") +
    scale_fill_manual(values=p_noinbred, name="Population/Line") + scale_color_manual(values=p_noinbred , name="Population/Line") +
    guides(alpha="none",size="none" , shape="none",fill="none",color="none")+
    theme(legend.position = "none", plot.title = element_blank(), legend.title=element_text(size=12))
  fig <- (dendro_with_tips.outlier + outlier.inv.gp$genotypes)  + plot_layout(widths = c(4, 3)) +  theme(plot.margin = margin(0, 0, 0, 0, "cm"), plot.title = element_blank())
  assign(paste("gp.plot.inv",chrom,"out", sep="."),fig, envir = globalenv())
}
# With Legend
inv.gp <- genotype_plot(vcf="inversion.outliers.vcf.gz",chr= "NC_035780.1",popmap = gp.popmap, start=1, end=94778955,cluster=TRUE,colour_scheme=c(col_pal[1],wes_palette("IsleofDogs1"))) 
#dendro_with_tips.dsa <- inv.gp$dendrogram +
#geom_text(aes(x=1:length(inv.gp$dendro_labels),y=-2.5,label=inv.gp$dendro_labels))
inv.gp$genotypes <- inv.gp$genotypes + theme(legend.text = element_text(size=10),legend.position = "right", legend.title=element_text(size=12)) + guides(fill = guide_legend(override.aes = list(size=4))) 
new_df <- data.frame(matrix(ncol=1,nrow=78, dimnames=list(NULL, c("IND"))))
new_df$IND <- inv.gp$dendro_labels
new_df <- join(new_df, popmap.ni)
dendro_with_tips.dsa <- inv.gp$dendrogram + 
  geom_jitter(data=new_df,aes(x=1:length(IND),y=-2.5,color=POP,fill=POP,shape=Type),size=3, alpha=0.95, height=1.5, width = 0)+
  scale_shape_manual(values=c(23,21),name="Origin") +
  scale_fill_manual(values=p_noinbred, name="Population/Line") + scale_color_manual(values=p_noinbred , name="Population/Line") +
  guides(fill = guide_legend(override.aes = list(size =6, shape = c(21,21,23,23,21,21,23,21,21,23,21,21,23),alpha = 0.95)), alpha=FALSE,size=FALSE , shape= guide_legend(override.aes = list(size =6))) + 
  guides(alpha="none",size="none" , shape="none",fill="none",color="none")+
  theme(legend.position = "right", plot.title = element_blank(), legend.title=element_text(size=12))
gp.plot.inv.1.out <- dendro_with_tips.dsa + inv.gp$genotypes +  plot_layout(widths = c(4, 3)) + theme(plot.margin = margin(0, 0, 0, 0, "cm"), plot.title = element_blank())
