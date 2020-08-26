EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q general -N 06_fst_SR_NB_cv30_no56invers -j oe -e cv30_no56invers_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 06_fst_by_pop_pair_SR_NB_cv30_no56invers.sh
exit 0;
