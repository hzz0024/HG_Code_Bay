# HG Code Bay

This repository is a collection of the the workflow and scripts during my postdoc at Hare lab

## WGS

scripts for whole-genome sequencing data trimming and processing (see Other_scripts->WGS->DelBay19_WGS_data_processing_HG.txt)

## Angsd 

scripts for Angsd analysis

## Detect the genetic signature of within-generation selection

### SGS

Scripts for SGS test, modified from Rey et al. 2020 (see R_scripts -> SGS_HG_local.R).  

A bash script allows to set the three parameters for the SGS test. They are 

1) num_split=100 # how many parts we want to divide the whole snp list  

2) num_sample=1000000 # sampling times

3) total_n=1732036 # total number of snps   

Scripts for p-value adjust (see R_scripts -> SGS_pvalue_adjust.R).

### Survival_vector

scripts for survival vector based outlier detection method 

### Fisher_exact

Scripts for Fisher's exact tests (see R_scripts -> Fisher_exact -> Fish_exact.R)

### Outlier evaluation

- plot deltap against start allele frequency

  1) Create the outlier list (e.g. fisher_exact_outlier.txt in the Data folder)    
  2) Using python script (extract.py) to extract the allele frequency values for the potential outliers    
  
  ```sh
  python3 extract.py
  ```

  3) Reveal the relationship between deltap and start p using R script Deltap_p_plot.R
