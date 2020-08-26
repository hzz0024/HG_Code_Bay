#!/bin/bash
for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_HC_ARN_'$genome'_no56invers -j oe -e '$genome'_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_HC_ARN_'$genome'_no56invers.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_HC_ARN_'$genome'_no56invers.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_HC_COH_'$genome'_no56invers -j oe -e '$genome'_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_HC_COH_'$genome'_no56invers.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_HC_COH_'$genome'_no56invers.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_HC_SR_'$genome'_no56invers -j oe -e '$genome'_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_HC_SR_'$genome'_no56invers.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_HC_SR_'$genome'_no56invers.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_HC_NB_'$genome'_no56invers -j oe -e '$genome'_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_HC_NB_'$genome'_no56invers.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_HC_NB_'$genome'_no56invers.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_ARN_COH_'$genome'_no56invers -j oe -e '$genome'_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_ARN_COH_'$genome'_no56invers.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_ARN_COH_'$genome'_no56invers.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_ARN_SR_'$genome'_no56invers -j oe -e '$genome'_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_ARN_SR_'$genome'_no56invers.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_ARN_SR_'$genome'_no56invers.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_ARN_NB_'$genome'_no56invers -j oe -e '$genome'_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_ARN_NB_'$genome'_no56invers.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_ARN_NB_'$genome'_no56invers.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_COH_SR_'$genome'_no56invers -j oe -e '$genome'_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_COH_SR_'$genome'_no56invers.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_COH_SR_'$genome'_no56invers.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_COH_NB_'$genome'_no56invers -j oe -e '$genome'_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_COH_NB_'$genome'_no56invers.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_COH_NB_'$genome'_no56invers.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_SR_NB_'$genome'_no56invers -j oe -e '$genome'_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_SR_NB_'$genome'_no56invers.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_SR_NB_'$genome'_no56invers.sh'
done
