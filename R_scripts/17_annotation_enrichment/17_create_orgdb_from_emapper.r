#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
library(tidyverse)

# Create a parser
p <- arg_parser("make OrgDB from emapper")

# Add command line arguments
p <- add_argument(p, "annotation", help="emapper annotation result", type="character")

# Parse the command line arguments
argv <- parse_args(p)

emapper <- read_delim(argv$annotation, 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                Gene_Symbol = X9, 
                GO = X10, 
                KO = X12, 
                Pathway = X13, 
                OG = X7, 
                Gene_Name = X8)

gene_info <- dplyr::select(emapper,  GID, Gene_Name) %>%
  dplyr::filter(!is.na(Gene_Name))

gene2go <- dplyr::select(emapper, GID, GO) %>%
  separate_rows(GO, sep = ',', convert = F) %>%
  filter(!is.na(GO)) %>%
  mutate(EVIDENCE = 'IEA') 

AnnotationForge::makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               maintainer='zhangsan <zhangsan@genek.tv>',
               author='zhangsan',
               outputDir="./",
               tax_id=0000,
               genus='M',
               species='y',
               goTable="go",
               version="1.0")

pkgbuild::build('.//org.My.eg.db', dest_path = ".")
