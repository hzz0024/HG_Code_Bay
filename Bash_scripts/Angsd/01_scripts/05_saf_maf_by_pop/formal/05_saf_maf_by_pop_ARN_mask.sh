module load angsd/0.931
###this script will work on bamfiles by population and calculate saf  & maf
# maybe edit
target="ARN"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_July/01_scripts/01_config.sh
N_IND=$(wc -l $ARN | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "Ouput can be used for depth evaluation with all individuals listed in "$ARN
echo "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"
echo "filter on allele frequency = "$MIN_MAF

angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doMajorMinor 3 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -b $ARN -sites WILD_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_mask -out "/scratch/hzz0024/DelBay19_July/05_saf_maf_by_pop/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask"
