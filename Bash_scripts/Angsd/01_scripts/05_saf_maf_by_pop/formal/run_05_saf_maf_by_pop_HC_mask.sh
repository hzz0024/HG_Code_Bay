EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q general -N 05_saf_HC_mask -j oe -e HC_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 05_saf_maf_by_pop_HC_mask.sh
exit 0;