#!/bin/bash
# Read a string ith spaces using for loop
for i in {1..10}
do 
    for pop in CH REF
    do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 07_ngsLD_'$pop'_cv30_invers_ch'$i' -j oe -e '$pop'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 07_ngsLD_'$pop'_cv30_no56invers_ch'$i'.sh\nexit 0;' >> formal/'run_07_ngsLD_'$pop'_cv30_no56invers_ch'$i'.sh'
    done
done
