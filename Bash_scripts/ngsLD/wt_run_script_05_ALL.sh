#!/bin/bash
# Read a string ith spaces using for loop
for i in {1..10}
do
    for pop in CH REF HC ARN COH SR NB
    do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 07_ngsLD_'$pop'_cv30_ch'$i' -j oe -e '$pop'.error -l nodes=1:ppn=16,mem=64GB,walltime=12:00:00 -m be -M $EMAIL -d $CWD -V 07_ngsLD_'$pop'_cv30_ch'$i'.sh\nexit 0;' >> formal/'run_07_ngsLD_'$pop'_cv30_ch'$i'.sh'
    done
done
