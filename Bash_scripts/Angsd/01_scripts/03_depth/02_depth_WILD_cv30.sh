#!/bin/bash
module load angsd/0.931

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
#maybe edit
# need to change $CHR at first echo part and angsd command
# also change the target below
target="WILD"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh
N_IND=$(wc -l $WILD | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "Ouput of -doQsDist can be used for depth evaluation with  all individuals listed in $WILD"
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is 70% of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

####Calculate the SAF, MAF and GL
angsd -P $NB_CPU \
-doMaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
-doDepth 1 \
-anc $ANC -remove_bads 1 -doQsDist 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 \
-minInd $MIN_IND -minMaf $MIN_MAF -maxDepth 2000 \
-b $WILD -out "/scratch/hzz0024/DelBay19_HG/03_depth/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"

#main features
# -P nb of threads
# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2) 
# -doMajorMinor 1 use the most frequent allele as major
# -anc provide a ancestral sequence = reference in our case
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

#main filters
#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ (minimum quality of reads?)
#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%
#filter on allele frequency -minMaf, set to 0.05 

