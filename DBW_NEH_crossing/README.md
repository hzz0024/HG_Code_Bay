# Data processing for population genomics with low-coverage whole-genome sequencing

This pipeline was built by Honggang Zhao with help and scripts from Nina's lab and many inputs from Claire Mérot github site.

# overview of the pipeline 
![overview_angsd_pipeline]()

## 00_DEPENDANCIES

A folder named DBW_NEH_crossing_data_processing, and run all commands from the DBW_NEH_crossing_data_processing folder

add angsd to the path in .bashrc

add the misc folder (containing RealSFS, theta stat etc to the path in .bashrc

export PATH="/home/camer78/Softwares/angsd2/angsd:$PATH"
export PATH="/home/camer78/Softwares/angsd2/angsd/misc:$PATH"

## 01_PREPARE_DATA

Obtain the raw fastq from sequencing output

```sh
for d in */ ; do
    echo "$d"
    cd $d
    mv * ..
    cd ..
done
```

create a fastq list

```sh
ls *_1.fq.gz > fastq_list.txt
```

create folders

```sh
mkdir sample_lists
mkdir adapter_clipped
mkdir bam_mtDNA
mkdir fastqc
mkdir bam
mkdir log
mkdir script
mkdir qual_filtered
mkdir reference
mkdir mt_mapped
```

## 02_raw_reads_filtering

```sh
#!/bin/bash
for i in {1..7}
do
    echo -e '#!/bin/sh\nstart=`date +%s`  ## date at start\nBASEDIR=/workdir/hz269/DBW_NEH_crossing_data_processing\nTRIMMOMATIC=/programs/trimmomatic/trimmomatic-0.39.jar\nSAMPLELIST=$BASEDIR/sample_lists/fastq_list_'$i'.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table.\nSAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.\nRAWFASTQDIR=$BASEDIR/raw_fastq/ # Path to raw fastq files.\nRAWFASTQSUFFIX1=_1.fq.gz # Suffix to raw fastq files. Use forward reads with paired-end data.\nRAWFASTQSUFFIX2=_2.fq.gz # Suffix to raw fastq files. Use reverse reads with paired-end data.\nADAPTERS=$BASEDIR/reference/NexteraPE-PE.fa # Path to a list of adapter/index sequences, copied from /programs/bbmap-38.86/resources/adapters.fa or /programs/trimmomatic/adapters/TruSeq3-PE-2.fa\n\n## Loop over each sample\nfor SAMPLEFILE in `cat $SAMPLELIST`; do\n\t## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library\n\tSAMPLE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 4`\n\tPOP_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 5`\n\tSEQ_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 3`\n\tLANE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 2`\n\tSAMPLE_UNIQ_ID=$SAMPLE_ID'\''_'\''$POP_ID'\''_'\''$SEQ_ID'\''_'\''$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely\n\techo "Sample: $SAMPLE_UNIQ_ID"\n\t## Extract data type from the sample table\n\tDATATYPE=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 6`\n\n\t## The input and output path and file prefix\n\tRAWFASTQ_ID=$RAWFASTQDIR$SAMPLEFILE\n\tSAMPLEADAPT=$BASEDIR'\''/adapter_clipped/'\''$SAMPLE_UNIQ_ID\n\n\t## Adapter clip the reads with Trimmomatic\n\t# The options for ILLUMINACLIP are: ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>\n\t# The MINLENGTH drops the read if it is below the specified length in bp\n\t# For definitions of these options, see http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf\n\t\tjava -jar $TRIMMOMATIC PE -threads 1 -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $RAWFASTQ_ID$RAWFASTQSUFFIX2 $SAMPLEADAPT'\''_adapter_clipped_f_paired.fastq.gz'\'' $SAMPLEADAPT'\''_adapter_clipped_f_unpaired.fastq.gz'\'' $SAMPLEADAPT'\''_adapter_clipped_r_paired.fastq.gz'\'' $SAMPLEADAPT'\''_adapter_clipped_r_unpaired.fastq.gz'\'' '\''ILLUMINACLIP:'\''$ADAPTERS'\'':2:30:10:1:true MINLENGTH:80'\''\n\ndone\n\nend=`date +%s` ## date at end\nruntime=$((end-start))\nhours=$((runtime / 3600))\nminutes=$(( (runtime % 3600) / 60 ))\nseconds=$(( (runtime % 3600) % 60 ))\necho "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"' >> formal/'2_trim_'$i'.sh'
done
```

## 04_remove_poly_g_tails

```sh
#!/bin/bash
for i in {1..7}
do
    echo -e '#!/bin/bash\nstart=`date +%s`  ## date at start\n## This script is used to quality filter and trim poly g tails. It can process both paired end and single end data.\nBASEDIR=/workdir/hz269/DBW_NEH_crossing_data_processing\nSAMPLELIST=$BASEDIR/sample_lists/fastq_list_'$i'.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table. \nSAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.\nFILTER=polyg # Type of filtering. Values can be: polyg (forced PolyG trimming only), quality (quality trimming, PolyG will be trimmed as well if processing NextSeq/NovaSeq data), or length (trim all reads to a maximum length)\nTHREAD=1 # Number of thread to use. Default is 10\nFASTP=/programs/fastp-0.20.0/bin/fastp ## Path to the fastp program. The default path is /workdir/programs/fastp_0.19.7/fastp\nMAXLENGTH=100 # Maximum length. This input is only relevant when FILTER=length, and its default value is 100.\n\n## Loop over each sample\nfor SAMPLEFILE in `cat $SAMPLELIST`; do\n\n\t## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library\n\tSAMPLE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 4`\n\tPOP_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 5`\n\tSEQ_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 3`\n\tLANE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 2`\n\tSAMPLE_UNIQ_ID=$SAMPLE_ID'\''_'\''$POP_ID'\''_'\''$SEQ_ID'\''_'\''$LANE_ID\n\techo "Sample: $SAMPLE_UNIQ_ID"\n\t## Extract data type from the sample table\n\tDATATYPE=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 6`\n\t## The input and output path and file prefix\n\tSAMPLEADAPT=$BASEDIR'\''/adapter_clipped/'\''$SAMPLE_UNIQ_ID\n\tSAMPLEQUAL=$BASEDIR'\''/qual_filtered/'\''$SAMPLE_UNIQ_ID\n\n\t## Trim polyg tail or low quality tail with fastp.\n\t# --trim_poly_g forces polyg trimming, --cut_right enables cut_right quality trimming\n\t# -Q disables quality filter, -L disables length filter, -A disables adapter trimming\n\t# Go to https://github.com/OpenGene/fastp for more information\n\t$FASTP --trim_poly_g --cut_right -L -A --thread $THREAD -i $SAMPLEADAPT'\''_adapter_clipped_f_paired.fastq.gz'\'' -I $SAMPLEADAPT'\''_adapter_clipped_r_paired.fastq.gz'\'' -o $SAMPLEQUAL'\''_adapter_clipped_qual_filtered_f_paired.fastq.gz'\'' -O $SAMPLEQUAL'\''_adapter_clipped_qual_filtered_r_paired.fastq.gz'\'' -h $SAMPLEQUAL'\''_adapter_clipped_fastp.html'\'' -j $SAMPLEQUAL'\''_adapter_clipped_fastp.json'\''\ndone\n\nend=`date +%s` ## date at end\nruntime=$((end-start))\nhours=$((runtime / 3600))\nminutes=$(( (runtime % 3600) / 60 ))\nseconds=$(( (runtime % 3600) % 60 ))\necho "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"' >> formal/'4_polyg_'$i'.sh'
done
```

## 05_map_to_mtDNA

```sh
for i in {1..7}
do
    echo -e '#!/bin/bash\nstart=`date +%s` \nBASEDIR=/workdir/hz269/DBW_NEH_crossing_data_processing\nBOWTIE=/programs/bowtie2-2.3.5.1-linux-x86_64/bowtie2\nSAMTOOLS=/programs/samtools-1.11/bin/samtools\nSAMPLELIST=$BASEDIR/sample_lists/fastq_list_'$i'.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table.\nSAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.\nFASTQDIR=$BASEDIR/qual_filtered/ # Path to the directory where fastq file are stored.\nFASTQSUFFIX1=_adapter_clipped_qual_filtered_f_paired.fastq.gz # Suffix to fastq files. Use forward reads with paired-end data.\nFASTQSUFFIX2=_adapter_clipped_qual_filtered_r_paired.fastq.gz # Suffix to fastq files. Use reverse reads with paired-end data.\nMAPPINGPRESET=very-sensitive # The pre-set option to use for mapping in bowtie2 (very-sensitive for end-to-end (global) mapping [typically used when we have a full genome reference], very-sensitive-local for partial read mapping that allows soft-clipping [typically used when mapping genomic reads to a transcriptome]\nREFERENCE=$BASEDIR/reference/CV30_mtDNA.fasta # Path to reference fasta file and file name\nREFNAME=CV30_mtDNA # Reference name to add to output files, e.g. gadMor2\n\n## Loop over each sample\nfor SAMPLEFILE in `cat $SAMPLELIST`; do\n\n\t## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library\n\tSAMPLE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 4`\n\tPOP_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 5`\n\tSEQ_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 3`\n\tLANE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 2`\n\tSAMPLE_UNIQ_ID=$SAMPLE_ID'\''_'\''$POP_ID'\''_'\''$SEQ_ID'\''_'\''$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely\n\n\t## Extract data type from the sample table\n\tDATATYPE=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 6`\n\n\t## The input and output path and file prefix\n\tSAMPLETOMAP=$FASTQDIR$SAMPLE_UNIQ_ID\n\tSAMPLEBAM=$BASEDIR'\''/bam_mtDNA/'\''$SAMPLE_UNIQ_ID\n\tSAMPLETOFASTQ=$BASEDIR'\''/mt_mapped/'\''$SAMPLE_UNIQ_ID\n\n\t## Define platform unit (PU), which is the lane number\n\tPU=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 2`\n\n\t## Define reference base name\n\tREFBASENAME="${REFERENCE%.*}"\n\n\t## Map reads to the reference\n\techo $SAMPLE_UNIQ_ID\n\n\t# Map the mtDNA\n\t$BOWTIE -q --phred33 --$MAPPINGPRESET --un-conc-gz $SAMPLETOFASTQ -p 1 -I 0 -X 1500 --fr --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 -S $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.sam'\''\n\n\t# Convert to bam file for storage (including all the mapped reads)\n\t$SAMTOOLS view -bS -F 4 $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.sam'\'' > $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.bam'\''\n\n\t# Remove the sam files\n\trm -f $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.sam'\''\n\n\t# Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20)\n\t$SAMTOOLS view -h -q 20 $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.bam'\'' | $SAMTOOLS view -buS - | $SAMTOOLS sort -o $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''_minq20_sorted.bam'\''\ndone\n\nend=`date +%s` ## date at end\nruntime=$((end-start))\nhours=$((runtime / 3600))\nminutes=$(( (runtime % 3600) / 60 ))\nseconds=$(( (runtime % 3600) % 60 ))\necho "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"' >> formal/'5_map_mtDNA_'$i'.sh'
done
```

## 06_map_to_genome

```sh
for i in {1..7}; do
echo -e '#!/bin/sh\nstart=`date +%s` \nBASEDIR=/workdir/hz269/DBW_NEH_crossing_data_processing\nBOWTIE=/programs/bowtie2-2.3.5.1-linux-x86_64/bowtie2\nSAMTOOLS=/programs/samtools-1.11/bin/samtools\nSAMPLELIST=$BASEDIR/sample_lists/fastq_list_'$i'.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table.\nSAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.\nFASTQDIR=$BASEDIR/mt_mapped/ # Path to the directory where fastq file are stored.\nFASTQSUFFIX1=.1 # Suffix to fastq files. Use forward reads with paired-end data.\nFASTQSUFFIX2=.2 # Suffix to fastq files. Use reverse reads with paired-end data.\nMAPPINGPRESET=very-sensitive # The pre-set option to use for mapping in bowtie2 (very-sensitive for end-to-end (global) mapping [typically used when we have a full genome reference], very-sensitive-local for partial read mapping that allows soft-clipping [typically used when mapping genomic reads to a transcriptome]\nREFERENCE=$BASEDIR/reference/CV30_masked.fasta # Path to reference fasta file and file name\nREFNAME=CV30_masked # Reference name to add to output files, e.g. gadMor2\n\n## Loop over each sample\nfor SAMPLEFILE in `cat $SAMPLELIST`; do\n\n\t## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library\n\tSAMPLE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 4`\n\tPOP_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 5`\n\tSEQ_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 3`\n\tLANE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 2`\n\tSAMPLE_UNIQ_ID=$SAMPLE_ID'\''_'\''$POP_ID'\''_'\''$SEQ_ID'\''_'\''$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely\n\n\t## Extract data type from the sample table\n\tDATATYPE=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | tr -d '\''\\r'\'' | cut -f 6`\n\n\t## The input and output path and file prefix\n\tSAMPLETOMAP=$FASTQDIR$SAMPLE_UNIQ_ID\n\tSAMPLEBAM=$BASEDIR'\''/bam/'\''$SAMPLE_UNIQ_ID\n\n\t## Define platform unit (PU), which is the lane number\n\tPU=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 2`\n\n\t## Define reference base name\n\tREFBASENAME="${REFERENCE%.*}"\n\n\t## Map reads to the reference\n\techo $SAMPLE_UNIQ_ID\n\n\t# Map the paired-end reads\n\t# We ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (loose their mate so they become single-end), these can be mapped in a separate step and the resulting bam files merged with the paired-end mapped reads\n\t$BOWTIE -q --phred33 --$MAPPINGPRESET -p 4 -I 0 -X 1500 --fr --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 -S $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.sam'\''\n\n\t## Convert to bam file for storage (including all the mapped reads)\n\t$SAMTOOLS view -bS -F 4 $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.sam'\'' > $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.bam'\''\n\trm -f $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.sam'\''\n\n\t## Filter the mapped reads (to onky retain reads with high mapping quality)\n\t# Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20)\n\t$SAMTOOLS view -h -q 20 $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.bam'\'' | $SAMTOOLS view -buS - | $SAMTOOLS sort -o $SAMPLEBAM'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''_minq20_sorted.bam'\''\n\ndone\n\nend=`date +%s` ## date at end\nruntime=$((end-start))\nhours=$((runtime / 3600))\nminutes=$(( (runtime % 3600) / 60 ))\nseconds=$(( (runtime % 3600) % 60 ))\necho "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"' >> formal/'6_map_'$i'.sh'
done
```

## 08_remove_duplicates

```sh
for i in {1..7}
do
    echo -e '#!/bin/bash\nstart=`date +%s` \n## This script is used to deduplicate bam files and clipped overlapping read pairs for paired end data. It can process both paired end and single end data.\nBAMLIST=/workdir/hz269/DBW_NEH_crossing_data_processing/sample_lists/bam_list_'$i'.txt # Path to a list of merged, deduplicated, and overlap clipped bam files. Full paths should be included.\nSAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.\nPICARD=/programs/picard-tools-2.19.2/picard.jar # Path to picard tools\nBAMUTIL=/programs/bamUtil/bam # Path to bamUtil\ncd bam\n\n## Loop over each sample\nfor SAMPLEBAM in `cat $BAMLIST`; do\n\n\t## Extract the file name prefix for this sample\n\tSAMPLESEQID=`echo $SAMPLEBAM | sed '\''s/_bt2_.*//'\'' | sed -e '\''s#.*/bam/\(\)#\1#'\''`\n\tSAMPLEPREFIX=`echo ${SAMPLEBAM%.bam}`\n\n\t## Remove duplicates and print dupstat file\n\t# We used to be able to just specify picard.jar on the CBSU server, but now we need to specify the path and version\n\tjava -Xmx60g -jar $PICARD MarkDuplicates I=$SAMPLEBAM O=$SAMPLEPREFIX'\''_dedup.bam'\'' M=$SAMPLEPREFIX'\''_dupstat.txt'\'' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true\n\n\t## Clip overlapping paired end reads (only necessary for paired end data)\n\t$BAMUTIL clipOverlap --in $SAMPLEPREFIX'\''_dedup.bam'\'' --out $SAMPLEPREFIX'\''_dedup_overlapclipped.bam'\'' --stats\n\ndone\n\nend=`date +%s` ## date at end\nruntime=$((end-start))\nhours=$((runtime / 3600))\nminutes=$(( (runtime % 3600) / 60 ))\nseconds=$(( (runtime % 3600) % 60 ))\necho "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"' >> formal/'8_dup_clip_'$i'.sh'
done
```

## 08_merge_bam_from_different_batch

```R
library(tidyverse)
setwd("~/Dropbox/Mac/Documents/HG/Github/HG_Code_Bay/DBW_NEH_crossing")

basedir <- "./"
refname <- "CV30_masked"

fastq_list <- read_lines(paste0(basedir, "sample_lists/fastq_list.txt"))
fastq_table <- read_tsv(paste0(basedir, "sample_lists/fastq_table.txt"))

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
```

## 09_realignment_around_indels

```sh
for i in {1..5}
do
    echo -e '#!/bin/bash\nstart=`date +%s`\nBAMLIST=/workdir/hz269/DelBay18_data_processing/sample_lists/DelBay18_dedup_overlapclipped_'$i'.list # Path to a list of merged, deduplicated, and overlap clipped bam files. Full paths should be included.\nBASEDIR=/workdir/hz269/DelBay18_data_processing\nREFERENCE=$BASEDIR/reference/CV30_masked.fasta # Path to reference fasta file and file name\nSAMTOOLS=/programs/samtools-1.11/bin/samtools # Path to samtools\nGATK=/programs/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar # Path to GATK\n\ncd bam\n\n## Loop over each sample\nfor SAMPLEBAM in `cat $BAMLIST`; do\n\nif [ -e $SAMPLEBAM'\''.bai'\'' ]; then\n\techo "the file already exists"\nelse\n\t## Index bam files\n\t$SAMTOOLS index $SAMPLEBAM\nfi\ndone\n\n## Realign around in-dels\n# This is done across all samples at once\n\n## Use an older version of Java\nexport JAVA_HOME=/usr/local/jdk1.8.0_121\nexport PATH=$JAVA_HOME/bin:$PATH\n\n## Create list of potential in-dels\n#if [ ! -f $BASEDIR'\''/bam/DelBay18_for_indel_realigner.intervals'\'' ]; then\n\t\t#java -Xmx120g -jar $GATK \\\n\t\t\t#-T RealignerTargetCreator \\\n\t\t\t#-R $REFERENCE \\\n\t\t\t#-I $BAMLIST \\\n\t\t\t#-o $BASEDIR'\''/bam/DelBay18_for_indel_realigner.intervals'\'' \\\n\t\t\t#-drf BadMate\n#fi\n\n## Run the indel realigner tool\njava -Xmx120g -jar $GATK \\\n\t-T IndelRealigner \\\n\t-R $REFERENCE \\\n\t-I $BAMLIST \\\n\t-targetIntervals $BASEDIR'\''/bam/DelBay18_for_indel_realigner.intervals'\'' \\\n\t--consensusDeterminationModel USE_READS  \\\n\t--nWayOut _realigned.bam\n\nend=`date +%s` ## date at end\nruntime=$((end-start))\nhours=$((runtime / 3600))\nminutes=$(( (runtime % 3600) / 60 ))\nseconds=$(( (runtime % 3600) % 60 ))\necho "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"' >> formal/'9_realign_'$i'.sh'
done
```

## script for summary statistics

```sh
wt_count_fastq.sh
for i in {1..8}
do
    echo -e '#!/bin/bash\nstart=`date +%s`\n#!/bin/bash\n\n## This script is used to count number of bases in raw, adapter clipped, and quality filtered fastq files. The result of this script will be stored in a nohup file.\nBASEDIR=/workdir/hz269/DelBay18_data_processing/\nSAMPLELIST=$BASEDIR/sample_lists/fastq_list_'$i'.txt\nSAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt\nRAWFASTQDIR=/workdir/hz269/DelBay18_data_processing/raw_fastq/\nQUALFILTERED=true # Whether the sample has gone through quality filtering. true or false\n\n# Create headers for the output\nif $QUALFILTERED; then\n\tprintf '\''sample_seq_id\\traw_reads\\traw_bases\\tadapter_clipped_bases\\tqual_filtered_bases\\n'\''\nelse\n\tprintf '\''sample_seq_id\\traw_reads\\traw_bases\\tadapter_clipped_bases\\n'\''\nfi\n\n# Loop over each sample in the sample table\nfor SAMPLEFILE in `cat $SAMPLELIST`; do\n\tRAWFASTQFILES=$RAWFASTQDIR$SAMPLEFILE'\''*.gz'\''  # The input path and file prefix\n\n\t# Count the number of reads in raw fastq files. We only need to count the forward reads, since the reverse will contain exactly the same number of reads. fastq files contain 4 lines per read, so the number of total reads will be half of this line number.\n\tRAWREADS=`zcat $RAWFASTQFILES | wc -l`\n\n\t# Count the number of bases in raw fastq files. We only need to count the forward reads, since the reverse will contain exactly the same number of bases. The total number of reads will be twice this count.\n\tRAWBASES=`zcat $RAWFASTQFILES | awk '\''NR%4==2'\'' | tr -d "\\n" | wc -m`\n\n\t# Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library\n\tPOP_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 5`\n\tSAMPLE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 4`\n\tSEQ_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 3`\n\tLANE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 2`\n\tSAMPLE_SEQ_ID=$SAMPLE_ID'\''_'\''$POP_ID'\''_'\''$SEQ_ID'\''_'\''$LANE_ID\n\n\t# Find all adapter clipped fastq files corresponding to this sample and store them in the object ADAPTERFILES.\n\tADAPTERFILES=$BASEDIR'\''adapter_clipped/'\''$SAMPLE_SEQ_ID'\''*.gz'\''\n\n\t# Count all bases in adapter clipped files.\n\tADPTERCLIPBASES=`zcat $ADAPTERFILES | awk '\''NR%4==2'\'' | tr -d "\\n" | wc -m`\n\n\t# If reads are quality filtered, count quality filtered files.\n\tif $QUALFILTERED; then\n\n\t\t# Find all quality trimmed fastq files corresponding to this sample and store them in the object QUALFILES.\n\t\tQUALFILES=$BASEDIR'\''qual_filtered/'\''$SAMPLE_SEQ_ID'\''*.gz'\''\n\n\t\t# Count bases in quality trimmed files.\n\t\tQUALFILTPBASES=`zcat $QUALFILES | awk '\''NR%4==2'\'' | tr -d "\\n" | wc -m`\n\n\t\t# Write the counts in appropriate order.\n\t\tprintf "%s\\t%s\\t%s\\t%s\\t%s\\n" $SAMPLE_SEQ_ID $((RAWREADS/4)) $RAWBASES $ADPTERCLIPBASES $QUALFILTPBASES\n\n\t\t# When reads are not quality filtered, directly write the output\n\telse\n\n\t\t# Write the counts in appropriate order.\n\t\tprintf "%s\\t%s\\t%s\\t%s\s\\n" $SAMPLE_SEQ_ID $((RAWREADS/4)) $RAWBASES $ADPTERCLIPBASES\n\n\tfi\n\ndone\n\nend=`date +%s` ## date at end\nruntime=$((end-start))\nhours=$((runtime / 3600))\nminutes=$(( (runtime % 3600) / 60 ))\nseconds=$(( (runtime % 3600) % 60 ))\necho "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"' >> formal/'count_fastq_'$i'.sh'
done

wt_count_unmerged.sh
for i in {1..8}
do
    echo -e '#!/bin/bash\nstart=`date +%s`\n#!/bin/bash\nBASEDIR=/workdir/hz269/DelBay18_data_processing/\nSAMPLELIST=$BASEDIR/sample_lists/fastq_list_'$i'.txt\nSAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt\nREFNAME=CV30_masked\nSAMTOOLS=/programs/samtools-1.11/bin/samtools\nMINQ=20 # Minimum mapping quality filter used in the sorting step. When left undefined, the sorted bam file will not be counted (e.g. when reads are not filtered in the sorting step).\n\nif [ -z "$MINQ" ]; then\n\tprintf '\''sample_id\\tmapped_bases\\n'\''\nelse\n\tprintf '\''sample_id\\tmapped_bases\\tqual_filtered_mapped_bases\\n'\''\nfi\n\nfor SAMPLEFILE in `cat $SAMPLELIST`; do\n\n\t# Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library\n\tPOP_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 5`\n\tSAMPLE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 4`\n\tSEQ_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 3`\n\tLANE_ID=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 2`\n\tSAMPLE_SEQ_ID=$SAMPLE_ID'\''_'\''$POP_ID'\''_'\''$SEQ_ID'\''_'\''$LANE_ID\n\n\t## Extract data type from the sample table\n\tDATATYPE=`grep -P "${SAMPLEFILE}\\t" $SAMPLETABLE | cut -f 6`\n\n\t## Count raw mapped bases\n\tRAWBAMFILE=$BASEDIR'\''bam/'\''$SAMPLE_SEQ_ID'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''.bam'\''\n\tMAPPEDBASES=`$SAMTOOLS stats $RAWBAMFILE | grep ^SN | cut -f 2- | grep "^bases mapped (cigar)" | cut -f 2`\n\n\tif [ -z "$MINQ" ]; then\n\t\tprintf "%s\\t%s\\n" $SAMPLE_SEQ_ID $MAPPEDBASES\n\telse\n\t\t## Count quality filtered mapped bases\n\t\tQUALFILTBAMFILE=$BASEDIR'\''bam/'\''$SAMPLE_SEQ_ID'\''_'\''$DATATYPE'\''_bt2_'\''$REFNAME'\''_minq'\''$MINQ'\''_sorted.bam'\''\n\t\tQUAFILTBASES=`$SAMTOOLS stats $QUALFILTBAMFILE | grep ^SN | cut -f 2- | grep "^bases mapped (cigar)" | cut -f 2`\n\t\tprintf "%s\\t%s\\t%s\\n" $SAMPLE_SEQ_ID $MAPPEDBASES $QUAFILTBASES\n\tfi\ndone\n\nend=`date +%s` ## date at end\nruntime=$((end-start))\nhours=$((runtime / 3600))\nminutes=$(( (runtime % 3600) / 60 ))\nseconds=$(( (runtime % 3600) % 60 ))\necho "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"' >> formal/'count_bam_unmerged_'$i'.sh'
done
```
