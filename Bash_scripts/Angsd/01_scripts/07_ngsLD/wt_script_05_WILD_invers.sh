# this script is used to create beagle files for ngsLD calculation
# Change the -site #### for CHR and WILD group
for i in {1..10}
do
    for pop in HC ARN COH SR NB
    do
    echo -e 'module load angsd/0.931\n#this script is used to create beagle files for ngsLD calculation\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_'$i'.list" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 3 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -b $'$pop' -sites WILD_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_no56invers_chr'$i'.list -out "/scratch/hzz0024/DelBay19_HG/07_ngsLD/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30_nochr56invers_ch'$i'"' >> formal/'07_ngsLD_'$pop'_cv30_no56invers_ch'$i'.sh'
    done
done
