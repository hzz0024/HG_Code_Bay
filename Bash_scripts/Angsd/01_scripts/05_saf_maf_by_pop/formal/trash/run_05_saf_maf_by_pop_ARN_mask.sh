EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q general -N 05_saf_ARN_mask -j oe -e ARN_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 05_saf_maf_by_pop_ARN_mask.sh
exit 0;
