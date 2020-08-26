# your email address from Hopper username

EMAIL=`whoami`"@auburn.edu";

# set the output directory to the current directory

CWD=`pwd`;

echo "Submitting a job to $NODES nodes with $CORES cores each..."

# submit a job to the system with the torque qsub command using

# all of the items we have just calculated...

qsub -q general -N test -j oe -e test.error -l nodes=1:ppn=20,mem=120GB,walltime=600:00:00 -m abe -M $EMAIL -d $CWD -V 03_saf_maf_gl_all.sh

exit 0;
