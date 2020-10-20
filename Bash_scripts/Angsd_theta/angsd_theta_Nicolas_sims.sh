## script from Nicolas simulation analysis scripts (03/25/20)
## calculating theta with angsd

## Get saf file
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80,160}; do
    nohup /workdir/programs/angsd0.931/angsd/angsd \
    -bam $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
    -out $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
    -doSaf 1 \
    -anc $BASE_DIR'slim/ancestral.fasta' \
    -GL 1 \
    -P 1 \
    > '/workdir/lcwgs-simulation/nohups/get_saf_neutral_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.nohup' &
  done
done
## Get SFS from saf (this need to be broken up into two separate runs due to memory limiations)
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80,160}; do
    nohup /workdir/programs/angsd0.931/angsd/misc/realSFS \
      $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
      -P 4 \
      >& $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.sfs' \
      2> '/workdir/lcwgs-simulation/nohups/get_sfs_neutral_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.nohup' < /dev/null &
  done
done
## Estimate theta
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80,160}; do
    nohup /workdir/programs/angsd0.931/angsd/angsd \
    -bam $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
    -out $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
    -doThetas 1 \
    -doSaf 1 \
    -pest $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.sfs' \
    -anc $BASE_DIR'slim/ancestral.fasta' \
    -GL 1 \
    -P 1 \
    > '/workdir/lcwgs-simulation/nohups/estimate_theta_neutral_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.nohup' &
  done
done
## Print per-SNP theta
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80,160}; do
    nohup /workdir/programs/angsd0.931/angsd/misc/thetaStat print \
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.tsv' &
  done
done
## Do fixed window theta
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80}; do
    nohup /workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat \
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    -win 10000 -step 10000 \
    -outnames $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.windowed_thetas.log' &
  done
done
## Do per-chromosome average theta
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80}; do
    nohup /workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat \
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    -outnames $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.average_thetas.idx' \
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.average_thetas.log' &
  done
done