library(dplyr)
library(stringr)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/22_enrichment/enrichment")
##############
# bed format #
##############
bed_format <- function(bim_file){
  bim = read.delim(bim_file, header = FALSE, sep='\t')
  head(bim)
  bim %>% 
    mutate(V1 = str_replace(V1, "1\\b", "NC_035780.1")) %>% 
    mutate(V1 = str_replace(V1, "2\\b", "NC_035781.1")) %>% 
    mutate(V1 = str_replace(V1, "3\\b", "NC_035782.1")) %>% 
    mutate(V1 = str_replace(V1, "4\\b", "NC_035783.1")) %>% 
    mutate(V1 = str_replace(V1, "5\\b", "NC_035784.1")) %>% 
    mutate(V1 = str_replace(V1, "6\\b", "NC_035785.1")) %>% 
    mutate(V1 = str_replace(V1, "7\\b", "NC_035786.1")) %>% 
    mutate(V1 = str_replace(V1, "8\\b", "NC_035787.1")) %>%
    mutate(V1 = str_replace(V1, "9\\b", "NC_035788.1")) %>% 
    mutate(V1 = str_replace(V1, "10\\b", "NC_035789.1"))  -> bim
  df = data.frame(bim$V1, bim$V4, bim$V4)
  return(df)
}

array_all <- bed_format("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.bim")
write.table(array_all, file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

bed_format <- function(bim_file){
  bim = read.delim(bim_file, header = FALSE, sep='\t')
  head(bim)
  bim %>% 
    mutate(V1 = str_replace(V1, "1\\b", "NC_035780.1")) %>% 
    mutate(V1 = str_replace(V1, "2\\b", "NC_035781.1")) %>% 
    mutate(V1 = str_replace(V1, "3\\b", "NC_035782.1")) %>% 
    mutate(V1 = str_replace(V1, "4\\b", "NC_035783.1")) %>% 
    mutate(V1 = str_replace(V1, "5\\b", "NC_035784.1")) %>% 
    mutate(V1 = str_replace(V1, "6\\b", "NC_035785.1")) %>% 
    mutate(V1 = str_replace(V1, "7\\b", "NC_035786.1")) %>% 
    mutate(V1 = str_replace(V1, "8\\b", "NC_035787.1")) %>%
    mutate(V1 = str_replace(V1, "9\\b", "NC_035788.1")) %>% 
    mutate(V1 = str_replace(V1, "10\\b", "NC_035789.1"))  -> bim
  df = data.frame(bim$V1, bim$V2, bim$V2)
  return(df)
}

PCAdapt <- bed_format("no_DBX1_UNC_MEH_MEW_MEH_n_340_PCAdapt_BP_q05_n_1103.bed")
write.table(PCAdapt, file = "no_DBX1_UNC_MEH_MEW_MEH_n_340_PCAdapt_BP_q05_n_1103.modified.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


##############
# annotation #
##############
# it is better to run  command below directly in the terminal
#system("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.bed -o genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.bed.gene.txt")

###################
# Extract GO data #
###################
# for whole SNPs
LOC_all <- read.delim("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.bed.gene.txt", header = FALSE, sep='\t')
dim(LOC_all)
LOC_all <- LOC_all[which(LOC_all$V6 != ""),]
dim(LOC_all)
# from NCBI protein database
emapper = read.delim("out.emapper.annotation.withLOC.uniq.txt", header = TRUE, sep='\t')
dim(emapper)
emapper$X.Locus
# extract the genes mapped by whole SNPs
extract_df <-emapper[which(emapper$X.Locus %in% paste0("LOC",LOC_all$V6)),]
dim(extract_df)
final_df<-extract_df[!(extract_df$GOs=="-"),]
dim(final_df)
write.table(final_df, file = "domestication.wholeSNPs.emapper.annotation.withLOC.uniq.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#######################
# enrichment analysis #
#######################

## load packages
library(tidyverse)
library(clusterProfiler)
library(stringr)
# set up working location
setwd("~/Dropbox/Mac/Documents/HG/Domestication/22_enrichment/enrichment")

dir.create('R_Library', recursive = T)

## prepare GO and KEGG lib
install.packages('./org.My.eg.db_1.0.tar.gz', 
                 repos = NULL,
                 lib = 'R_Library') 

# load OrgDB
library(org.My.eg.db, lib = 'R_Library')

# load gene ID
SGS_single_SNP <- read.delim("no_DBX1_UNC_MEH_MEW_MEH_n_340_PCAdapt_BP_q05_n_1103.txt", header = F, sep='\t')$V1

# GO enrichment
SGS_singleSNP_go <- enrichGO(gene = SGS_single_SNP,
                             OrgDb = org.My.eg.db,
                             keyType = 'GID',
                             ont = 'ALL',
                             qvalueCutoff = 0.05,
                             pvalueCutoff = 0.05)

SGS_singleSNP_go_df <- as.data.frame(SGS_singleSNP_go)
write.table(SGS_singleSNP_go_df, file = "SGS_wild_go_df.txt", sep = "\t", quote = FALSE,
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
graph2ppt(file="SGS_wild_go.pptx", width=10, height=10)

cnetplot(SGS_singleSNP_go, 
         showCategory = 5,
         node_label = "all", # category | gene | all | none
         circular = TRUE, 
         colorEdge = TRUE)

## -----------------------------------------------------------

###################################
#########  KEGG enrichment ########
###################################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/22_enrichment/enrichment")
emapper <- read_delim('./domestication.wholeSNPs.emapper.annotation.withLOC.uniq.txt', 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                KO = X12, 
                Pathway = X13)

# load gene ID
SGS_single_SNP <- read.delim("no_DBX1_UNC_MEH_MEW_MEH_n_340_PCAdapt_BP_q05_n_1103.txt", header = F, sep='\t')$V1

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
                              pvalueCutoff = 0.2,
                              qvalueCutoff = 0.2)

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

