EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q fastfat -N 04_theta_WILD_cv30 -j oe -e WILD_cv30.error -l nodes=1:ppn=16,mem=200GB,walltime=240:00:00,flags=ADVRES:liuzhan_ff -m be -M $EMAIL -d $CWD -V 04_saf_maf_theta_WILD_cv30.sh
exit 0;
