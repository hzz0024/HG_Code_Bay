EMAIL=`whoami`"@auburn.edu";
CWD=`pwd`;
qsub -q general -N 07_ngsLD_ARN_cv30_invers_ch1 -j oe -e ARN_mask.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V 07_ngsLD_ARN_cv30_no56invers_ch1.sh
exit 0;
