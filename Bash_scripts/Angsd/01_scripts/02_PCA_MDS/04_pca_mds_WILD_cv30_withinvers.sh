###this script will work on bamfiles by population and calculate saf  & maf
# maybe edit
target="WILD"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /workdir/hz269/01_scripts/01_config.sh
N_IND=$(wc -l $WILD | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "Ouput can be used for depth evaluation with all individuals listed in "$WILD
echo "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"
echo "filter on allele frequency = "$MIN_MAF

/programs/angsd_20180926/angsd/angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -doIBS 1 -makeMatrix 1 -doCov 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth 962 -SNP_pval 1e-6 -b $WILD -out "/workdir/hz269/02_PCA_MDS/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30_withinvers"

#main features
# -P nb of threads
# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2)
# -doMajorMinor 1 use the most frequent allele as major
# -anc provide a ancestral sequence = reference in our case
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

#main filters
#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ minimum quality of reads
#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%
#filter on allele frequency -minMaf, set to 0.05