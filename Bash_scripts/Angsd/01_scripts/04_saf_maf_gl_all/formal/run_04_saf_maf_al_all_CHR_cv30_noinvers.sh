EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q fastfat -N 04_saf_CHR_cv30 -j oe -e CHR_cv30.error -l nodes=1:ppn=16,mem=200GB,walltime=240:00:00,flags=ADVRES:liuzhan_ff -m be -M $EMAIL -d $CWD -V 04_saf_maf_al_all_CHR_cv30_noinvers.sh
exit 0;
