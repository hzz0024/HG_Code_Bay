#!/bin/bash

for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_HC_ARN.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_HC_ARN_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_HC_COH.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_HC_COH_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_HC_SR.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_HC_SR_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_HC_NB.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_HC_NB_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_ARN_COH.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_ARN_COH_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_ARN_SR.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_ARN_SR_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_ARN_NB.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_ARN_NB_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_COH_SR.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_COH_SR_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_COH_NB.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_COH_NB_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_SR_NB.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_SR_NB_'$genome'.sh'
done
