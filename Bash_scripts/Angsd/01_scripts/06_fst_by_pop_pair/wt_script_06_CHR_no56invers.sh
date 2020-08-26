#!/bin/bash
for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_July/06_fst_by_pop_pair pop_ch.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_CHR_'$genome'_no56invers.sh'
done

