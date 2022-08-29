library(dplyr)
library(stringr)

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/17_annotation_enrichment/06_final_format/")
##############
# bed format #
##############
bed_format <- function(snp_list){
  snp_list_df = read.delim(snp_list, header = FALSE, sep='\t')
  head(snp_list_df)
  df = data.frame(snp_list_df$V1, snp_list_df$V2, snp_list_df$V2)
  return(df)
}

SNP_all <- bed_format("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.txt")
write.table(SNP_all, file = "All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(SNP_all, file = "Challenge_20_all_minq20_minmq30_CV30_masked.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

##############
# annotation #
##############
# it is better to run command below directly in the terminal
system("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.bed -o All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.bed.gene.txt")

###################
# Extract GO data #
###################
# for whole SNPs

LOC_all <- read.table("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.bed.gene.txt", header = FALSE, sep = "\t", quote="", col.names = paste0("V",seq_len(6)), fill = TRUE)
#LOC_all <- read.delim("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.bed.gene.txt", header = FALSE, sep='\t')
dim(LOC_all)
LOC_all <- LOC_all[which(LOC_all$V6 != ""),]
dim(LOC_all)
# from NCBI protein database
#emapper = read.delim("out.emapper.annotation.withLOC.uniq.txt", header = TRUE, sep='\t')
emapper = read.table("out.emapper.annotation.withLOC.uniq.txt", header = TRUE, sep='\t', quote="", comment.char="#", col.names = paste0("V",seq_len(21)), fill = TRUE)
dim(emapper)
#emapper$X.Locus
# extract the genes mapped by whole SNPs
extract_df <-emapper[which(emapper$V1 %in% paste0("LOC",LOC_all$V6)),]
#extract_df <-emapper[which(emapper$X.Locus %in% paste0("LOC",LOC_all$V6)),]
dim(extract_df)
extract_df_filter <- extract_df[which(extract_df$V10 != "-"),]
dim(extract_df_filter)
write.table(extract_df_filter, file = "DB.wholeSNPs.emapper.annotation.withLOC.uniq.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#######################
# Format outliers     #
#######################
library(stringr)
library(dplyr)
format_FDR <- function(fname){
  output_list <- read.delim(fname, header = TRUE, sep='\t')
  df <- data.frame(output_list$chr, output_list$pos, output_list$pos)
  df <- as.data.frame(df)
  write.table(df, file = paste0(fname, ".bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#fname = 'permutation_outliers_n_45968.txt.bed.gene.txt'
format_LOC <- function(fname){
  gene_list <- read.delim(fname, header = FALSE, sep='\t')
  dim(gene_list)
  gene_list <- gene_list[which(gene_list$V6 != ""),]
  dim(gene_list)
  all_list <- str_split(gene_list$V6, ";")
  unique_list <- sort(paste0("LOC", unique(unlist(all_list))))
  write.table(unique_list, file = paste0(fname, ".unique.LOC.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

format_FDR("18_SGS_HC_NB_FDR_outlier.list")
system(paste("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t 18_SGS_HC_NB_FDR_outlier.list.bed -o 18_SGS_HC_NB_FDR_outlier.list.bed.gene.txt"))

format_FDR("19_SGS_HC_NB_FDR_outlier.list")
system(paste("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t 19_SGS_HC_NB_FDR_outlier.list.bed -o 19_SGS_HC_NB_FDR_outlier.list.bed.gene.txt"))

format_FDR("21_SGS_HC_NB_FDR_outlier.list")
system(paste("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t 21_SGS_HC_NB_FDR_outlier.list.bed -o 21_SGS_HC_NB_FDR_outlier.list.bed.gene.txt"))

format_FDR("19_SGS_Sur_Ref_FDR_outlier.list")
system(paste("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t 19_SGS_Sur_Ref_FDR_outlier.list.bed -o 19_SGS_Sur_Ref_FDR_outlier.list.bed.gene.txt"))

format_FDR("20_SGS_Sur_Ref_FDR_outlier.list")
system(paste("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t 20_SGS_Sur_Ref_FDR_outlier.list.bed -o 20_SGS_Sur_Ref_FDR_outlier.list.bed.gene.txt"))

format_FDR("permutation_outliers_n_45968.txt")
system(paste("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t permutation_outliers_n_45968.txt.bed -o permutation_outliers_n_45968.txt.bed.gene.txt"))

format_LOC("18_SGS_HC_NB_FDR_outlier.list.bed.gene.txt")
format_LOC("19_SGS_HC_NB_FDR_outlier.list.bed.gene.txt")
format_LOC("21_SGS_HC_NB_FDR_outlier.list.bed.gene.txt")
format_LOC("19_SGS_Sur_Ref_FDR_outlier.list.bed.gene.txt")
format_LOC("20_SGS_Sur_Ref_FDR_outlier.list.bed.gene.txt")
format_LOC("permutation_outliers_n_45968.txt.bed.gene.txt")
format_LOC("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.bed.gene.txt")

#######################
# enrichment analysis #
#######################

## load packages
library(tidyverse)
library(clusterProfiler)
library(stringr)
library(export)
# set up working location
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/17_annotation_enrichment")

dir.create('R_Library', recursive = T)

## prepare GO and KEGG lib
install.packages('./org.My.eg.db_1.0.tar.gz', 
                 repos = NULL,
                 lib = 'R_Library') 

# load OrgDB
library(org.My.eg.db, lib = 'R_Library')

# load gene ID
SGS_single_SNP <- read.delim("permutation_outliers_n_45968.txt.bed.gene.txt.unique.LOC.txt", header = F, sep='\t')$V1

# GO enrichment
SGS_singleSNP_go <- enrichGO(gene = SGS_single_SNP,
                             OrgDb = org.My.eg.db,
                             keyType = 'GID',
                             ont = 'ALL',
                             qvalueCutoff = 0.05,
                             pvalueCutoff = 0.05)


SGS_singleSNP_go_df <- as.data.frame(SGS_singleSNP_go)
write.table(SGS_singleSNP_go_df, file = "permutation_outliers_n_45968_go_df.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

head(SGS_singleSNP_go_df)

require(DOSE)
require(enrichplot)
require(viridis)
## plot

barplot(SGS_singleSNP_go, showCategory = 10, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale="free")

# key results for GO
dotplot(SGS_singleSNP_go, showCategory = 10, split="ONTOLOGY") + 
  scale_color_viridis(option="plasma") + 
  facet_grid(ONTOLOGY~., scale="free")+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=60)) + 
  labs(size="Counts",col="FDR")
graph2ppt(file="permutation_outliers_n_45968_go.pptx", width=10, height=10)

cnetplot(SGS_singleSNP_go, 
         showCategory = 5,
         node_label = "all", # category | gene | all | none
         circular = TRUE, 
         colorEdge = TRUE)

###################################
#########  KEGG enrichment ########
###################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/17_annotation_enrichment")
emapper <- read_delim('./DB.wholeSNPs.emapper.annotation.withLOC.uniq.txt', 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                KO = X12, 
                Pathway = X13)

# load gene ID
SGS_single_SNP <- read.delim("20_SGS_Sur_Ref_FDR_outlier.list.bed.gene.txt.unique.LOC.txt", header = F, sep='\t')$V1

pathway2gene <- dplyr::select(emapper, Pathway, GID) %>%
  separate_rows(Pathway, sep = ',', convert = F) %>%
  filter(str_detect(Pathway, 'ko')) %>%
  mutate(Pathway = str_remove(Pathway, 'ko'))

library(magrittr)
get_path2name <- function(){
  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
  keggpathid2name.df[,1] %<>% gsub("path:map", "", .)
  colnames(keggpathid2name.df) <- c("path_id","path_name")
  return(keggpathid2name.df)
}

pathway2name <- get_path2name()

SGS_singleSNP_ekp <- enricher(SGS_single_SNP,
                              TERM2GENE = pathway2gene,
                              TERM2NAME = pathway2name,
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)

SGS_singleSNP_ekp_df <- as.data.frame(SGS_singleSNP_ekp)

head(SGS_singleSNP_ekp_df)

## plot
barplot(SGS_singleSNP_ekp, showCategory = 10)

dotplot(SGS_singleSNP_ekp, showCategory = 10) + 
  scale_color_viridis(option="plasma") + 
  labs(size="Counts",col="FDR") +
  scale_size_area(max_size = 14)

graph2ppt(file="SGS_singleSNP_KEGG.pptx", width=8, height=8)

cnetplot(SGS_singleSNP_ekp, 
         #foldChange = geneList, 
         showCategory = 3,
         node_label = "all", # category | gene | all | none
         circular = FALSE, 
         colorEdge = TRUE)
## -----------------------------------------------------------
save(SGS_singleSNP_go, SGS_singleSNP_ekp, file = './SGS_singleSNP_enrich.rdata')


###################################
#########  Shared GO term  ########
###################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/17_annotation_enrichment")

format_go <- function(go_df_file) {
  go_df <- read.delim(go_df_file, header = T, sep='\t')
  go_id <- go_df$ID
  return(go_id)
}

HC_NB_18 <- format_go("18_SGS_HC_NB_go_df.txt")
HC_NB_19 <- format_go("19_SGS_HC_NB_go_df.txt")
HC_NB_21 <- format_go("21_SGS_HC_NB_go_df.txt")
Sur_Ref_19 <- format_go("19_SGS_Sur_Ref_go_df.txt")
Sur_Ref_20 <- format_go("20_SGS_Sur_Ref_go_df.txt")
Permutation <- format_go("permutation_outliers_n_45968_go_df.txt")
library("ggvenn")
x <- list(
  A = HC_NB_18, 
  B = HC_NB_19, 
  C = HC_NB_21,
  D = Sur_Ref_19,
  E = Sur_Ref_20,
  F = Permutation
)
names(x) <- c("HC_NB_18","HC_NB_19","HC_NB_21", "Sur_Ref_19", "2020_SGS_outliers", "2020_permutation_outliers")
ggvenn(
  x, columns = c("HC_NB_18", "HC_NB_19", "HC_NB_21"),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5,set_name_size = 6
)

ggvenn(
  x, columns = c("Sur_Ref_19", "Sur_Ref_20"),
  fill_color = c("#CD534CFF", "#bc5090"),
  stroke_size = 0.5,set_name_size = 6
)

ggvenn(
  x, columns = c("2020_SGS_outliers", "2020_permutation_outliers"),
  fill_color = c("#CD534CFF", "#EFC000FF"),
  stroke_size = 0.5,set_name_size = 6
)

ggvenn(
  x, columns = c("Sur_Ref_19", "HC_NB_19"),
  fill_color = c("#CD534CFF", "#EFC000FF"),
  stroke_size = 0.5,set_name_size = 6
)

length(intersect(Sur_Ref_20, intersect(Sur_Ref_19, intersect(HC_NB_21, intersect(HC_NB_18,  HC_NB_19)))))
shared <- intersect(Sur_Ref_20, intersect(Sur_Ref_19, intersect(HC_NB_21, intersect(HC_NB_18,  HC_NB_19))))
go_df <- read.delim("18_SGS_HC_NB_go_df.txt", header = T, sep='\t')

all_shared <- go_df[which(go_df$ID %in% shared),]

shared <- intersect(Sur_Ref_20, Permutation)
go_df <- read.delim("permutation_outliers_n_45968_go_df.txt", header = T, sep='\t')
all_shared <- go_df[which(go_df$ID %in% shared),]

library(gridExtra)
a <- tableGrob(all_shared[1:3])
grid.table(all_shared[1:3])


