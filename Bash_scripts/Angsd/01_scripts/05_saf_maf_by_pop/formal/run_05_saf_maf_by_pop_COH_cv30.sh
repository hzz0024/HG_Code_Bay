EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q general -N 05_saf_COH_cv30 -j oe -e COH_cv30.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 05_saf_maf_by_pop_COH_cv30.sh
exit 0;