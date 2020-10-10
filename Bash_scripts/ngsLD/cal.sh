# extrac the start and end information for each chromosome
for i in {0..9}; do
echo "chr $i"
cat ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 |grep "NC_03578$i.1" | awk -F ' ' '{print $2}' | head -n 1
cat ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 |grep "NC_03578$i.1" | awk -F ' ' '{print $2}' | tail -n 1
done
