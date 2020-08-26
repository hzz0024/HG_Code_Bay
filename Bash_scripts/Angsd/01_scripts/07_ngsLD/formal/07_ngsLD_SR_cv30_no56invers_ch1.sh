module load angsd/0.931
#this script is used to create beagle files for ngsLD calculation
# maybe edit
target="SR"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_1.list" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh
N_IND=$(wc -l $SR | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 3 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -b $SR -sites WILD_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_no56invers_chr1.list -out "/scratch/hzz0024/DelBay19_HG/07_ngsLD/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30_nochr56invers_ch1"
