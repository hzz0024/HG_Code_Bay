#!/bin/bash
# Read a string with spaces using for loop
for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_CHR_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_CHR_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_CHR_'$genome'.sh'
done
