EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q fastfat -N 04_saf_WILD_mask -j oe -e WILD_mask.error -l nodes=1:ppn=16,mem=200GB,walltime=240:00:00,flags=ADVRES:liuzhan_ff -m be -M $EMAIL -d $CWD -V 04_saf_maf_al_all_WILD_mask.sh
exit 0;