library(tidyverse)

basedir <- "~/Documents/GitHub/HG_Code_Bay/DelBay19/"
refname <- "CV30_masked"

fastq_list <- read_lines(paste0(basedir, "sample_lists/fastq_list.txt"))
fastq_table <- read_tsv(paste0(basedir, "sample_lists/fastq_table.tsv"))

#Subset the table to the samples we're currently interested in analyzing
target_samples <- filter(fastq_table, prefix %in% fastq_list)

## Find all duplicated samples
duplicated_samples <- (target_samples$sample_id)[duplicated(target_samples$sample_id)] %>% unique()

# Write cd
write_lines(c("BASEDIR=$1", "cd $BASEDIR'/bam'\n"), paste0(basedir, "scripts/merge_bam.sh"))


## Loop through all duplicated samples 
for (i in 1:length(duplicated_samples)){
  dup_sample_dat <- filter(fastq_table, sample_id==duplicated_samples[i])
  
  ## Extract the bam file names from the unmerged sample table
  input <- dup_sample_dat %>%
    mutate(unmerged_bam = paste(sample_id, population, seq_id, lane_number, data_type, "bt2", refname, "minq20_sorted.bam", sep = "_")) %>% 
    # We are reconstructing the $SAMPLE_SEQ_ID that is the first part of the separate bam file names
    .$unmerged_bam
  
  output <- paste(dup_sample_dat[1,"sample_id"], dup_sample_dat[1, "population"], "merged_bt2", refname, "minq20_sorted.bam", sep = "_")
  
  ## Paste together the command line
  #merging_script[i+1] <- paste("samtools merge", as.character(output), input[1], input[2], sep = " ")
  write_lines(paste("samtools merge", as.character(output), input[1], input[2], sep = " "), paste0(basedir, "scripts/merge_bam.sh"), append = TRUE)
  
  # Also write a target bam list that we'll use for downstream looping over merged bam files
  if (i == 1){
    write_lines(paste(dup_sample_dat[1,"sample_id"], dup_sample_dat[1, "population"], "merged", sep = "_"), paste0(basedir, "sample_lists/merged_bam_list.txt"))
  }
  
  if (i > 1){
    write_lines(paste(dup_sample_dat[1,"sample_id"], dup_sample_dat[1, "population"], "merged", sep = "_"), paste0(basedir, "sample_lists/merged_bam_list.txt"), append = TRUE)
  }
  
}