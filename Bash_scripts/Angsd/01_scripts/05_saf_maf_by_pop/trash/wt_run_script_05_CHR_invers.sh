#!/bin/bash
# Read a string with spaces using for loop
for pop in CH REF
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 05_saf_'$pop'_mask_invers -j oe -e '$pop'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 05_saf_maf_by_pop_'$pop'_mask_no56invers.sh\nexit 0;' >> formal/'run_05_saf_maf_by_pop_'$pop'_mask_no56invers.sh'
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 05_saf_'$pop'_cv30_invers -j oe -e '$pop'_cv30.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 05_saf_maf_by_pop_'$pop'_cv30_no56invers.sh\nexit 0;' >> formal/'run_05_saf_maf_by_pop_'$pop'_cv30_no56invers.sh'
done
