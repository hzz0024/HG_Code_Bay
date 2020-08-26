EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q general -N 04_theta_CHR_cv30 -j oe -e CHR_cv30.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 04_saf_maf_theta_CHR_cv30.sh
exit 0;
