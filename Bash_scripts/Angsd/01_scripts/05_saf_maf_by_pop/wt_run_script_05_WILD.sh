#!/bin/bash
# Read a string with spaces using for loop
for pop in HC ARN COH SR NB
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 05_saf_'$pop'_mask -j oe -e '$pop'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 05_saf_maf_by_pop_'$pop'_mask.sh\nexit 0;' >> formal/'run_05_saf_maf_by_pop_'$pop'_mask.sh'
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 05_saf_'$pop'_cv30 -j oe -e '$pop'_cv30.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 05_saf_maf_by_pop_'$pop'_cv30.sh\nexit 0;' >> formal/'run_05_saf_maf_by_pop_'$pop'_cv30.sh'
done
