#!/bin/bash
# Read a string with spaces using for loop
for POP in CH REF HC ARN COH SR NB
do
    echo -e 'module load angsd/0.931\n\ntarget=\x27'$POP'\x27\nNB_CPU=20 #change accordingly\nsource /scratch/hzz0024/DelBay19_Sep/01_scripts/01_config.sh\nBASE_DIR=\x27/scratch/hzz0024/DelBay19_Sep/09_theta/\x27\n\nangsd -P $NB_CPU -doSaf 1 -doThetas 1 -GL 1 -fold 1 -pest $BASE_DIR$target\x27.sfs\x27 -anc $ANC -b $'$POP' -sites ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 -out $BASE_DIR$target\n\n## Print per-SNP theta\n/tools/angsd-0.931/misc/thetaStat print $BASE_DIR$target\x27.thetas.idx\x27\n\n## Do fixed window theta\n/tools/angsd-0.931/misc/thetaStat do_stat $BASE_DIR$target\x27.thetas.idx\x27 -win 10000 -step 10000 -outnames $BASE_DIR$target\x27.thetas.window.idx\x27\n\n## Do per-chromosome average theta\n/tools/angsd-0.931/misc/thetaStat do_stat $BASE_DIR$target\x27.thetas.idx\x27 -outnames $BASE_DIR$target\x27.thetas.average.idx\x27' >> formal/'09_theta_by_pop_'$POP'.sh'
done

####### an example

module load angsd/0.931

target='CH'
NB_CPU=20 #change accordingly
source /scratch/hzz0024/DelBay19_Sep/01_scripts/01_config.sh
BASE_DIR='/scratch/hzz0024/DelBay19_Sep/09_theta/'

angsd -P $NB_CPU -doSaf 1 -doThetas 1 -GL 1 -fold 1 -pest $BASE_DIR$target'.sfs' -anc $ANC -b $CH -sites ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 -out $BASE_DIR$target

## Print per-SNP theta
/tools/angsd-0.931/misc/thetaStat print $BASE_DIR$target'.thetas.idx'

## Do fixed window theta
/tools/angsd-0.931/misc/thetaStat do_stat $BASE_DIR$target'.thetas.idx' -win 10000 -step 10000 -outnames $BASE_DIR$target'.thetas.window.idx'

## Do per-chromosome average theta
/tools/angsd-0.931/misc/thetaStat do_stat $BASE_DIR$target'.thetas.idx' -outnames $BASE_DIR$target'.thetas.average.idx'
