module load angsd/0.31
for i in {0..9}; do
    j=`expr $i + 1`
    echo "Number of SNPs in chr$j is:"
    grep "NC_03578$i.1" ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 > "ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_chr"$j
    cat "ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_chr"$j | wc -l
    angsd sites index "ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_chr"$j
done
