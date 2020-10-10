# CH n_ind = 50
# REF n_ind = 48
# HC n_ind = 48
# ARN n_ind = 47
# COH n_ind = 44
# SR n_ind = 48
# NB n_ind = 48

for i in {5..5}
do
    for pop in SR # REF HC ARN COH SR NB
    do
	for k in 25 50
	do
	    zcat $pop"_maf0.05_pctind0.7_cv30_ch"$i".mafs.gz" | cut -f 1,2 | tail -n +2 > $pop"_maf0.05_pctind0.7_cv30_ch"$i"_pos.txt"
	    cnt=$(cat $pop"_maf0.05_pctind0.7_cv30_ch"$i"_pos.txt" | wc -l)
    	    /programs/ngsLD/ngsLD \
	    --geno $pop"_maf0.05_pctind0.7_cv30_ch"$i".beagle.gz" \
	    --pos $pop"_maf0.05_pctind0.7_cv30_ch"$i"_pos.txt" \
	    --n_ind 48 \
	    --n_sites $cnt \
	    --out $pop"_ngsLD_ch"$i"_k"$k"_output" \
	    --probs \
	    --max_kb_dist 0 \
            --max_snp_dist $k \
	    --min_maf 0.05 \
	    --n_threads 8
	done
    done
done
