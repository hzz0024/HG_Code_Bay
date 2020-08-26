#!/bin/bash
for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_HC_ARN_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_HC_ARN_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_HC_ARN_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_HC_COH_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_HC_COH_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_HC_COH_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_HC_SR_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_HC_SR_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_HC_SR_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_HC_NB_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_HC_NB_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_HC_NB_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_ARN_COH_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_ARN_COH_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_ARN_COH_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_ARN_SR_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_ARN_SR_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_ARN_SR_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_ARN_NB_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_ARN_NB_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_ARN_NB_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_COH_SR_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_COH_SR_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_COH_SR_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_COH_NB_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_COH_NB_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_COH_NB_'$genome'.sh'
done

for genome in cv30 mask
do
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N 06_fst_SR_NB_'$genome' -j oe -e '$genome'_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_SR_NB_'$genome'.sh\nexit 0;' >> formal/'run_06_fst_by_pop_pair_SR_NB_'$genome'.sh'
done
