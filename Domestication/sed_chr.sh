#for i in by_pop_0.05_pctind0.7_maxdepth3.snps; do
for i in salinity_2032113_outlier_BF10.txt; do
#for i in salinity_2032113_outlier_BF20.txt; do
sed -i.bak 's/Chr2/NC_035781.1/g;s/Chr3/NC_035782.1/g;s/Chr4/NC_035783.1/g;s/Chr5/NC_035784.1/g;s/Chr6/NC_035785.1/g;s/Chr7/NC_035786.1/g;s/Chr8/NC_035787.1/g;s/Chr9/NC_035788.1/g;s/Chr10/NC_035789.1/g' $i
#sed -i.bak 's/NC_035780.1/Chr1/g;s/NC_035781.1/Chr2/g;s/NC_035782.1/Chr3/g;s/NC_035783.1/Chr4/g;s/NC_035784.1/Chr5/g;s/NC_035785.1/Chr6/g;s/NC_035786.1/Chr7/g;s/NC_035787.1/Chr8/g;s/NC_035788.1/Chr9/g;s/NC_035789.1/Chr10/g;s/NC_007175.2/11/g' $i
#sed -i .bak 's/NC_035780.1/1/g;s/NC_035781.1/2/g;s/NC_035782.1/3/g;s/NC_035783.1/4/g;s/NC_035784.1/5/g;s/NC_035785.1/6/g;s/NC_035786.1/7/g;s/NC_035787.1/8/g;s/NC_035788.1/9/g;s/NC_035789.1/10/g;s/NC_007175.2/11/g' $i
done

for i in salinity_2032113_outlier_BF10.txt; do
#for i in salinity_2032113_outlier_BF20.txt; do
sed -i.bak 's/Chr1/NC_035780.1/g' $i
#sed -i.bak 's/NC_035780.1/Chr1/g;s/NC_035781.1/Chr2/g;s/NC_035782.1/Chr3/g;s/NC_035783.1/Chr4/g;s/NC_035784.1/Chr5/g;s/NC_035785.1/Chr6/g;s/NC_035786.1/Chr7/g;s/NC_035787.1/Chr8/g;s/NC_035788.1/Chr9/g;s/NC_035789.1/Chr10/g;s/NC_007175.2/11/g' $i
#sed -i .bak 's/NC_035780.1/1/g;s/NC_035781.1/2/g;s/NC_035782.1/3/g;s/NC_035783.1/4/g;s/NC_035784.1/5/g;s/NC_035785.1/6/g;s/NC_035786.1/7/g;s/NC_035787.1/8/g;s/NC_035788.1/9/g;s/NC_035789.1/10/g;s/NC_007175.2/11/g' $i
done
