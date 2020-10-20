module load angsd/0.931

## Step 1 Get SFS from saf
BASE_DIR='/scratch/hzz0024/DelBay19_Sep/09_theta/'
for POP in {ALL}; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 1 '$POP' sfs starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/realSFS \
    $BASE_DIR$POP'_maf0.05_pctind0.7_cv30.saf.idx' \
    -P 16 \
    -fold 1 \
    > $BASE_DIR$POP'.sfs'
done

## Step 2a Estimate theta
for POP in {ALL}; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 2a '$POP' theta estimate starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/realSFS saf2theta \
    $BASE_DIR$POP'_maf0.05_pctind0.7_cv30.saf.idx' \
    -sfs $BASE_DIR$POP'.sfs' \
    -fold 1 \
    -outname $BASE_DIR$POP
done

## Step 2b Print per-SNP theta
for POP in {ALL}; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 2a '$POP' per-SNP theta print starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/thetaStat print \
    $BASE_DIR$POP'.thetas.idx' \
    > $BASE_DIR$POP'.thetas.tsv'
done

## Step 3a do fixed window theta
for POP in {ALL}; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 3a '$POP' fixed window theta starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/thetaStat do_stat \
    $BASE_DIR$POP'.thetas.idx' \
    -win 5000 -step 5000 \
    -outnames $BASE_DIR$POP'.thetas.window.idx'
done

## Step 3b do per-chromosome average theta
for POP in {ALL}; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 3b '$POP' per-chromosome average theta starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/thetaStat do_stat \
    $BASE_DIR$POP'.thetas.idx' \
    -outnames $BASE_DIR$POP'.thetas.average.idx'
done