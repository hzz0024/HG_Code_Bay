EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q general -N 04_saf_CHR_mask -j oe -e CHR_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 04_saf_maf_al_all_CHR_mask.sh
exit 0;
