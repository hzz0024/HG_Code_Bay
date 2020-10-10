###this script will produce files for PCA and MDS plot using the identified SNP list
# maybe edit
target="WILD"
NB_CPU=32 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /workdir/hz269/01_scripts/01_config.sh
N_IND=$(wc -l $WILD | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "Ouput can be used for depth evaluation with all individuals listed in "$WILD
echo "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"
echo "filter on allele frequency = "$MIN_MAF

/programs/angsd_20180926/angsd/angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 3 -anc $ANC -b $WILD -sites ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 -out "/workdir/hz269/05_saf_maf_by_pop/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"
