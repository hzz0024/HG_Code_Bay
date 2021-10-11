library(plyr)
library(dplyr)
library(tidyverse)
library(readtext)

#set working directory
setwd("~/Documents/HG/DelBay19_adult/17_annotation/Clusterprofiler_steps/03_format")

#import NCBI protein table
#this table includes "YP" mitochondrial transcripts
LOC <- read.table("proteins_398_327504.txt", header=TRUE, sep="\t", dec=".", quote="") %>% 
  dplyr::select("GeneID", "Locus", "Protein_product")
names(LOC)[names(LOC)=="Protein_product"] <- "query" # rename the protein ID column to match the B2G annotation

#import original B2G annotation file created by Proestou and Sullivan using the standard blast2go workflow (GO weight parameter increased from 5 to 15) for all 60201 transcripts.
#This annotation does not include "YP" mitochondrial transcripts. 
#ann <- read.csv("all_60201_gcf_002022765_2_c_virginica_3_0_protein_WP_15_kevinInterpro_customExport.csv")
ann <- read.delim("out.emapper.annotations.txt", header = TRUE, sep='\t')
#merge the annotation and LOC dataframes so that we only keep the 60,201 non-mitochondrial transcripts. 
ann_with_LOC <- merge(ann, LOC, by="query",all.x=TRUE)
#rearrange the column order so that the geneIDs and locus IDs are first.
ann_with_LOC <- ann_with_LOC[, c(1, 23, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 )]

# Dimenesions
dim(ann_with_LOC)
# remove duplicate rows with dplyr
ann_with_LOC_uniq <- ann_with_LOC %>% 
  # Base the removal on the "Locus" column
  distinct(Locus, .keep_all = TRUE)
# Dimenesions
dim(ann_with_LOC_uniq)
# only keep the LOC for downstream formatting
ann_with_LOC_uniq <- ann_with_LOC_uniq[, c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 )]
#export
write.table(ann_with_LOC_uniq, file = "out.emapper.annotation.withLOC.uniq.txt", sep = "\t", row.names = FALSE, quote = F)
