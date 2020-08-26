EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q general -N 06_fst_CHR_mask_no56invers -j oe -e mask_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_CHR_mask_nochr56invers.sh
exit 0;
