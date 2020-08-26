#!/bin/bash
# Read a string with spaces using for loop
for pop in CHR
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q fastfat -N 04_saf_'$pop'_mask -j oe -e '$pop'_mask.error -l nodes=1:ppn=16,mem=200GB,walltime=240:00:00,flags=ADVRES:liuzhan_ff -m be -M $EMAIL -d $CWD -V 04_saf_maf_al_all_'$pop'_mask_noinvers.sh\nexit 0;' >> formal/'run_04_saf_maf_al_all_'$pop'_mask_noinvers.sh'
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q fastfat -N 04_saf_'$pop'_cv30 -j oe -e '$pop'_cv30.error -l nodes=1:ppn=16,mem=200GB,walltime=240:00:00,flags=ADVRES:liuzhan_ff -m be -M $EMAIL -d $CWD -V 04_saf_maf_al_all_'$pop'_cv30_noinvers.sh\nexit 0;' >> formal/'run_04_saf_maf_al_all_'$pop'_cv30_noinvers.sh'
done
