## load packages
library(tidyverse)
library(clusterProfiler)
library(stringr)
library(viridis)
library(export)
# set up working location
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/17_annotation_enrichment")

dir.create('R_Library', recursive = T)

## prepare GO and KEGG lib
install.packages("./org.Cv.eg.db_1.0.tar.gz", 
                 repos = NULL,
                 lib = 'R_Library') 

# load OrgDB
library(org.Cv.eg.db, lib = 'R_Library')

# load gene ID
SGS_single_SNP <- read.delim("HC_NB_FDR_outlier.gene.txt", header = F, sep='\t')$V1

# GO enrichment
SGS_singleSNP_go <- enrichGO(gene = SGS_single_SNP,
         OrgDb = org.Cv.eg.db,
                   keyType = 'GID',
                   ont = 'ALL',
                   qvalueCutoff = 0.05,
                   pvalueCutoff = 0.05)


SGS_singleSNP_go_df <- as.data.frame(SGS_singleSNP_go)
write.table(SGS_singleSNP_go_df, file = "HC_NB_FDR_outlier_GO_output.txt", sep = "\t", quote = FALSE,
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
graph2ppt(file="SGS_singleSNP_challenge.pptx", width=10, height=10)

cnetplot(SGS_singleSNP_go, 
         showCategory = 5,
         node_label = "all", # category | gene | all | none
         circular = TRUE, 
         colorEdge = TRUE)

## -----------------------------------------------------------

###################################
#########  KEGG enrichment ########
###################################
setwd("~/Documents/HG/DelBay19_adult/17_annotation/Clusterprofiler_steps/05_formal_run")
emapper <- read_delim('./out.emapper.annotation.withLOC.uniq.txt', 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                KO = X12, 
                Pathway = X13)

# load gene ID
SGS_single_SNP <- read.delim("SGS_outliers_single_snp.txt", header = F, sep='\t')$V1

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
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1)

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
