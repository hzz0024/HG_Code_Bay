###############################################################
###################### lfmm manhattan plot ####################
###############################################################
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)

rm(list=ls())
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/04_lfmm/all_env")

# original function with facet_grid and x axis features 
# ps_fdr <- function(fname) {
#   #fname = "CD11_k1" 
#   dat = read.delim(paste0(fname, ".pvalcalib"), header = T, sep=' ')
#   snps <- read.delim("by_pop_0.05_pctind0.7_maxdepth3.snps", header = T, sep=' ')
#   dat$id = paste0(snps$chromo,'_',snps$position)
#   dat$chr1 = snps$chromo
#   dat$pos = snps$position/1e+7
#   dat$fdr = p.adjust(dat$V1_1, method = 'BH')
#   head(dat)
#   fdr=dat$fdr
#   fdr_cut <- 0.05
#   outlier <- as.data.frame(cbind(dat$chr[fdr<fdr_cut], dat$pos[fdr<fdr_cut], dat$fdr[fdr<fdr_cut], dat$id[fdr<fdr_cut]))
#   colnames(outlier) <- c("chr","pos","fdr","id")
#   outlier <- outlier[with(outlier, order(chr, pos)),]
#   outlier$fdr <- as.numeric(outlier$fdr)
#   outlier$pos <- as.numeric(outlier$pos)
#   write.table(outlier, paste0(fname, ".FDR.outlier.txt"), col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)
#   # count how many SNPs with FDR < fdr_cut
#   message(paste0("FDR < ", fdr_cut, ": ", length(outlier$id)))
#   # ggplot for FDR values
#   dat$pos <- as.numeric(dat$pos)
#   dat$fdr <- as.numeric(dat$fdr)
#   dat$chr1 <- paste0("Chr",dat$chr1)
#   dat$chr <- factor(dat$chr1, levels = c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10'))
#   fig <- ggplot(outlier, aes(x=pos, y=-log10(fdr)))+ 
#     geom_point(aes(colour = cut(-log10(fdr), c(-Inf, 4, 5, Inf))),size = 1, show.legend = F)+
#     scale_color_manual(name = "fdr",
#                        values = c("(-Inf,4]" = "grey",
#                                   "(4,5]" = "orange",
#                                   "(5, Inf]" = "red"))+
#     theme_classic()+
#     facet_grid(cols = vars(chr), scales = "free_x", space="free_x") +
#     theme(text = element_text(size=20)) +
#     theme(axis.text.x = element_text(color = "grey20", size = 10)) +
#     geom_hline(aes(yintercept =5), linetype="dotted", size=0.8, col="red", show.legend = FALSE)+
#     geom_hline(aes(yintercept =4), linetype="dotted", size=0.8, col="orange", show.legend = FALSE) +
#     xlab("Position (10M bp)") +
#     ylab(paste0("-log10(p) ", fname))
#   # store the fig as variable
#   assign(paste0(fname,".out"), fig, envir = globalenv())
# }

ps_fdr_top <- function(fname) {
  #fname = "CD11_k1" 
  dat = read.delim(paste0(fname, ".pvalcalib"), header = T, sep=' ')
  snps <- read.delim("by_pop_0.05_pctind0.7_maxdepth3.snps", header = T, sep=' ')
  dat$id = paste0(snps$chromo,'_',snps$position)
  dat$chr1 = snps$chromo
  dat$pos = snps$position/1e+7
  dat$fdr = p.adjust(dat$V1_1, method = 'BH')
  head(dat)
  fdr=dat$fdr
  fdr_cut <- 0.05
  outlier <- as.data.frame(cbind(dat$chr[fdr<fdr_cut], dat$pos[fdr<fdr_cut], dat$fdr[fdr<fdr_cut], dat$id[fdr<fdr_cut]))
  colnames(outlier) <- c("chr","pos","fdr","id")
  outlier <- outlier[with(outlier, order(chr, pos)),]
  outlier$fdr <- as.numeric(outlier$fdr)
  outlier$pos <- as.numeric(outlier$pos)
  write.table(outlier, paste0(fname, ".FDR.outlier.txt"), col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)
  # count how many SNPs with FDR < fdr_cut
  message(paste0("FDR < ", fdr_cut, ": ", length(outlier$id)))
  # ggplot for FDR values
  dat$pos <- as.numeric(dat$pos)
  dat$fdr <- as.numeric(dat$fdr)
  dat$chr1 <- paste0("Chr",dat$chr1)
  dat$chr <- factor(dat$chr1, levels = c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10'))
  fig <-ggplot(dat, aes(x=pos, y=-log10(fdr)))+ 
    geom_point(aes(colour = cut(-log10(fdr), c(-Inf, 4, 5, Inf))),size = 1, show.legend = F)+
    scale_color_manual(name = "fdr",
                       values = c("(-Inf,4]" = "grey",
                                  "(4,5]" = "orange",
                                  "(5, Inf]" = "red"))+
    theme_classic()+
    facet_grid(cols = vars(chr), scales = "free_x", space="free_x") +
    theme(text = element_text(size=15)) +
    #theme(axis.text.x = element_text(color = "grey20", size = 10)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank()) + 
    geom_hline(aes(yintercept =5), linetype="dotted", size=0.8, col="red", show.legend = FALSE)+
    geom_hline(aes(yintercept =4), linetype="dotted", size=0.8, col="orange", show.legend = FALSE) +
    xlab("Position (10M bp)") +
    ylab(paste0("-log10(p) ", fname))
  
  # store the fig as variable
  assign(paste0(fname,".out"), fig, envir = globalenv())
}

ps_fdr_mid <- function(fname) {
  #fname = "CD11_k1" 
  dat = read.delim(paste0(fname, ".pvalcalib"), header = T, sep=' ')
  snps <- read.delim("by_pop_0.05_pctind0.7_maxdepth3.snps", header = T, sep=' ')
  dat$id = paste0(snps$chromo,'_',snps$position)
  dat$chr1 = snps$chromo
  dat$pos = snps$position/1e+7
  dat$fdr = p.adjust(dat$V1_1, method = 'BH')
  head(dat)
  fdr=dat$fdr
  fdr_cut <- 0.05
  outlier <- as.data.frame(cbind(dat$chr[fdr<fdr_cut], dat$pos[fdr<fdr_cut], dat$fdr[fdr<fdr_cut], dat$id[fdr<fdr_cut]))
  colnames(outlier) <- c("chr","pos","fdr","id")
  outlier <- outlier[with(outlier, order(chr, pos)),]
  outlier$fdr <- as.numeric(outlier$fdr)
  outlier$pos <- as.numeric(outlier$pos)
  write.table(outlier, paste0(fname, ".FDR.outlier.txt"), col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)
  # count how many SNPs with FDR < fdr_cut
  message(paste0("FDR < ", fdr_cut, ": ", length(outlier$id)))
  # ggplot for FDR values
  dat$pos <- as.numeric(dat$pos)
  dat$fdr <- as.numeric(dat$fdr)
  dat$chr1 <- paste0("Chr",dat$chr1)
  dat$chr <- factor(dat$chr1, levels = c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10'))
  fig <-ggplot(dat, aes(x=pos, y=-log10(fdr)))+ 
    geom_point(aes(colour = cut(-log10(fdr), c(-Inf, 4, 5, Inf))),size = 1, show.legend = F)+
    scale_color_manual(name = "fdr",
                       values = c("(-Inf,4]" = "grey",
                                  "(4,5]" = "orange",
                                  "(5, Inf]" = "red"))+
    theme_classic()+
    facet_grid(cols = vars(chr), scales = "free_x", space="free_x") +
    theme(strip.text.x = element_blank())+
    theme(text = element_text(size=15)) +
    #theme(axis.text.x = element_text(color = "grey20", size = 10)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank()) + 
    geom_hline(aes(yintercept =5), linetype="dotted", size=0.8, col="red", show.legend = FALSE)+
    geom_hline(aes(yintercept =4), linetype="dotted", size=0.8, col="orange", show.legend = FALSE) +
    xlab("Position (10M bp)") +
    ylab(paste0("-log10(p) ", fname))
  
  # store the fig as variable
  assign(paste0(fname,".out"), fig, envir = globalenv())
}

ps_fdr_bot <- function(fname) {
  #fname = "CD11_k1" 
  dat = read.delim(paste0(fname, ".pvalcalib"), header = T, sep=' ')
  snps <- read.delim("by_pop_0.05_pctind0.7_maxdepth3.snps", header = T, sep=' ')
  dat$id = paste0(snps$chromo,'_',snps$position)
  dat$chr1 = snps$chromo
  dat$pos = snps$position/1e+7
  dat$fdr = p.adjust(dat$V1_1, method = 'BH')
  head(dat)
  fdr=dat$fdr
  fdr_cut <- 0.05
  outlier <- as.data.frame(cbind(dat$chr[fdr<fdr_cut], dat$pos[fdr<fdr_cut], dat$fdr[fdr<fdr_cut], dat$id[fdr<fdr_cut]))
  colnames(outlier) <- c("chr","pos","fdr","id")
  outlier <- outlier[with(outlier, order(chr, pos)),]
  outlier$fdr <- as.numeric(outlier$fdr)
  outlier$pos <- as.numeric(outlier$pos)
  write.table(outlier, paste0(fname, ".FDR.outlier.txt"), col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)
  # count how many SNPs with FDR < fdr_cut
  message(paste0("FDR < ", fdr_cut, ": ", length(outlier$id)))
  # ggplot for FDR values
  dat$pos <- as.numeric(dat$pos)
  dat$fdr <- as.numeric(dat$fdr)
  dat$chr1 <- paste0("Chr",dat$chr1)
  dat$chr <- factor(dat$chr1, levels = c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10'))
  fig <-ggplot(dat, aes(x=pos, y=-log10(fdr)))+ 
    geom_point(aes(colour = cut(-log10(fdr), c(-Inf, 4, 5, Inf))),size = 1, show.legend = F)+
    scale_color_manual(name = "fdr",
                       values = c("(-Inf,4]" = "grey",
                                  "(4,5]" = "orange",
                                  "(5, Inf]" = "red"))+
    theme_classic()+
    facet_grid(cols = vars(chr), scales = "free_x", space="free_x") +
    theme(strip.text.x = element_blank())+
    theme(text = element_text(size=15)) +
    theme(axis.text.x = element_text(color = "grey20", size = 10)) +
    geom_hline(aes(yintercept =5), linetype="dotted", size=0.8, col="red", show.legend = FALSE)+
    geom_hline(aes(yintercept =4), linetype="dotted", size=0.8, col="orange", show.legend = FALSE) +
    xlab("Position (10M bp)") +
    ylab(paste0("-log10(p) ", fname))
  
  # store the fig as variable
  assign(paste0(fname,".out"), fig, envir = globalenv())
}

for(i in c("SA_k1", "SA_k2", "SA_k3")){
  ps_fdr_top(i)
}

for(i in c("CD5_k1", "CD5_k2", "CD5_k3", "CD7_k1", "CD7_k2", "CD7_k3", "CD9_k1", "CD9_k2", "CD9_k3", "CD11_k1", "CD11_k2", "CD11_k3")){
  ps_fdr_mid(i)
}

for(i in c("MAX10_k1", "MAX10_k2", "MAX10_k3")){
  ps_fdr_bot(i)
}

jpeg("./Figures/lfmm_K1_Manhattan.jpg", width = 18, height = 12, units = 'in', res = 300)
wrap_plots(SA_k1.out, CD5_k1.out,CD7_k1.out,CD9_k1.out,CD11_k1.out,MAX10_k1.out, nrow = 6) + plot_layout(guides = 'collect')
dev.off()

jpeg("./Figures/lfmm_k2_Manhattan.jpg", width = 18, height = 12, units = 'in', res = 300)
wrap_plots(SA_k2.out, CD5_k2.out,CD7_k2.out,CD9_k2.out,CD11_k2.out,MAX10_k2.out, nrow = 6) + plot_layout(guides = 'collect')
dev.off()

jpeg("./Figures/lfmm_K3_Manhattan.jpg", width = 18, height = 12, units = 'in', res = 300)
wrap_plots(SA_k3.out, CD5_k3.out,CD7_k3.out,CD9_k3.out,CD11_k3.out,MAX10_k3.out, nrow = 6) + plot_layout(guides = 'collect')
dev.off()


# SA_k1 <- ps_fdr("SA_k1")
# SA_k2 <- ps_fdr("SA_k2")
# SA_k3 <- ps_fdr()
# CD5_k1 <- ps_fdr("CD5_k1")
# CD5_k2 <- ps_fdr("CD5_k2")
# CD5_k3 <- ps_fdr("CD5_k3")
# CD7_k1 <- ps_fdr("CD7_k1")
# CD7_k2 <- ps_fdr("CD7_k2")
# CD7_k3 <- ps_fdr("CD7_k3")
# CD9_k1 <- ps_fdr("CD9_k1")
# CD9_k2 <- ps_fdr("CD9_k2")
# CD9_k3 <- ps_fdr("CD9_k3")
# CD11_k1 <- ps_fdr("CD11_k1")
# CD11_k2 <- ps_fdr("CD11_k2")
# CD11_k3 <- ps_fdr("CD11_k3")
# MAX10_k1 <- ps_fdr("MAX10_k1")
# MAX10_k2 <- ps_fdr("MAX10_k2")
# MAX10_k3 <- ps_fdr("MAX10_k3")


################################################
#######  process fisher outlier for ngsLD ######
################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/fisher_outlier")
########################################################
# Step 1: format the ps file and convert it to bed format
########################################################
# must sort the position for ngsLD running
format_bed <- function(pname, distance){
  #pname = "format/REF19_CHR19_NB_HC_out_0.05_fish.txt"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat = dat[with(dat, order(V1, V2)),]
  bed_list <- paste0(dat$V1, "\t" , dat$V2, "\t", dat$V2+distance)
  print(head(bed_list))
  write.table(bed_list, paste0(pname, ".bed"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_bed("./format/REF19_CHR19_NB_HC_out_0.05_fish.txt", 250)

########################################################
# Step 2: using bedtools to merge intervals
########################################################
# path in /workdir/hz269/DelBay_all_angsd_final/15_LD_prunning/process_bed_file

# for f in *.bed; do
# echo $f
# cat $f | wc -l
# bedtools merge -i $f > $f'.merged.txt'
# cat $f'.merged.txt' | wc -l
# done

# ps_Del19_ARN_COH.txt.bed
# 1113
# 922
# ps_Del19_challenge.txt.bed
# 2185
# 1854
# ps_Del19_HC_NB.txt.bed
# 2073
# 1696
# ps_Del20_challenge.txt.bed
# 3117
# 2525

########################################################
# Step 3: convert the bed format to rf input for Angsd
########################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/03_global/")
format_rf <- function(pname){
  #pname = "WILD.unlinked.id"
  dat = read.delim(pname, header = FALSE, sep=':')
  dat = dat[with(dat, order(V1, V2)),]
  angsd_list <- paste0(dat$V1, ":" , dat$V2, "-", dat$V2)
  write.table(angsd_list, paste0(strsplit(pname, split = ".txt")[[1]][1], ".rf.txt"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_rf("format/WILD.unlinked.id")

########################################################

       