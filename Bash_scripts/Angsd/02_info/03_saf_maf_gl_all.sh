#!/bin/bash
module load angsd/0.931

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
#maybe edit
target=CHR
NB_CPU=20 #change accordingly
REGIONS="-rf 02_info/regions_25kb_100snp.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh
N_IND=$(wc -l $CHR | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * 0.7)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*}

echo "Calculate the SAF, MAF and GL for all individuals listed in $CHR"
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is 70% of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

####Calculate the SAF, MAF and GL
angsd -P $NB_CPU \
-doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
-doDepth 1 -dumpCounts 1 \
-anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 \
-minInd $MIN_IND -minMaf $MIN_MAF -setMinDepth 49 -setMaxDepth 347 \
-b $CHR -rf chr.list -out /scratch/hzz0024/DelBay19_HG/03_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv

#main features
# -P nb of threads -nQueueSize maximum waiting in memory (necesary to optimize CPU usage
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

#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles
#order the sites file by chromosome names 
#makes a region file matching the sites files and with same order
#index sites file
echo "from the maf file, extract a list of SNP chr, positoin, major all, minor all"
cd /scratch/hzz0024/DelBay19_HG/03_saf_maf_gl_all
#gunzip "$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv.mafs.gz 
cat "$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv.mafs | tail -n +2 > FILE.tmp && mv FILE.tmp "$target"_snplist
awk '{print $1,$2,$3,$4}' "$target"_snplist > "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col
angsd sites index "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col

#Eric propose a much shorter version using bash to cut the 4 columns of the mafs.gz. but then angsd is unable to index it
#I can't find the problem, so I came back to older solution with R which works
#gunzip -c 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind_"$PERCENT_IND".mafs.gz | cut -f -4 > 02_info/sites_all_maf"$MIN_MAF"_pctind_"$PERCENT_IND"

