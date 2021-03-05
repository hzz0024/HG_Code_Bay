#### Theta estimates

Two methods have been proposed to calculate the theta values, see this post [Thetas: new saf2theta gives different results than old angsd -doThetas](https://github.com/ANGSD/angsd/issues/336).

For initial saf file, some parameter setting examples can be found [https://github.com/atigano/Peromyscus_eremicus_genome/tree/master/variant_calling_ANGSD](https://github.com/atigano/Peromyscus_eremicus_genome/tree/master/variant_calling_ANGSD)

- Old method

1) produce the saf file

```sh
#!/bin/bash
# Read a string with spaces using for loop

for pop in CH REF HC COH ARN SR NB
do
    echo -e 'module load angsd/0.931\n###this script will work on bamfiles by population and calculate saf  & maf\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_Sep/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -doCounts 1 -GL 1 -doMajorMinor 1 -anc $ANC -remove_bads 1 -minQ 20 -minMapQ 30 -b $'$pop' $REGIONS -out "/scratch/hzz0024/DelBay19_Sep/11_SFS/"$target"_SFS_cv30"' >> formal/'09_theta_'$pop'_cv30.sh'
done

# an example 
angsd -P $NB_CPU -doMaf 1 -dosaf 1 -doCounts 1 -GL 1 -doMajorMinor 1 -anc $ANC -setMaxDepth 200 -minQ 20 -minMapQ 30 -b $CH $REGIONS -out "/scratch/hzz0024/DelBay19_Sep/11_SFS/"$target"_SFS_cv30"
```
2) produce the sfs file

```sh
## Get SFS from saf
module load angsd/0.931

BASE_DIR='/scratch/hzz0024/DelBay19_Sep/11_SFS/'
for POP in {CH,REF,HC,ARN,COH,SR,NB}; do
    echo $POP' sfs starts'
    /tools/angsd-0.931/misc/realSFS \
      $BASE_DIR$POP'_SFS_cv30.saf.idx' \
      -P 16 -fold 1 \
      > $BASE_DIR$POP'.sfs'
done
```

Delete extra 0 values in sfs file. For example, CH pop has a total of 50 samples, and 101 if unfold values. We need to delete the extra 50 0 values and only retain the 51 values.

3) Estimate theta

```sh
#!/bin/bash
# Read a string with spaces using for loop
for POP in CH REF HC ARN COH SR NB
do
    echo -e 'module load angsd/0.931\n\ntarget=\x27'$POP'\x27\nNB_CPU=20 #change accordingly\nsource /scratch/hzz0024/DelBay19_Sep/01_scripts/01_config.sh\nBASE_DIR=\x27/scratch/hzz0024/DelBay19_Sep/09_theta/\x27\n\nangsd -P $NB_CPU -doSaf 1 -doThetas 1 -GL 1 -fold 1 -pest $BASE_DIR$target\x27.sfs\x27 -anc $ANC -b $'$POP' -sites ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 -out $BASE_DIR$target\n\n## Print per-SNP theta\n/tools/angsd-0.931/misc/thetaStat print $BASE_DIR$target\x27.thetas.idx\x27\n\n## Do fixed window theta\n/tools/angsd-0.931/misc/thetaStat do_stat $BASE_DIR$target\x27.thetas.idx\x27 -win 10000 -step 10000 -outnames $BASE_DIR$target\x27.thetas.window.idx\x27\n\n## Do per-chromosome average theta\n/tools/angsd-0.931/misc/thetaStat do_stat $BASE_DIR$target\x27.thetas.idx\x27 -outnames $BASE_DIR$target\x27.thetas.average.idx\x27' >> formal/'09_theta_by_pop_'$POP'.sh'
done

# example
module load angsd/0.931

target='CH'
NB_CPU=20 #change accordingly
source /scratch/hzz0024/DelBay19_Sep/01_scripts/01_config.sh
BASE_DIR='/scratch/hzz0024/DelBay19_Sep/09_theta/'

angsd -P $NB_CPU -doSaf 1 -doThetas 1 -GL 1 -fold 1 -pest $BASE_DIR$target'.sfs' -anc $ANC -b $CH -sites ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 -out $BASE_DIR$target

## Print per-SNP theta
/tools/angsd-0.931/misc/thetaStat print $BASE_DIR$target'.thetas.idx'

## Do fixed window theta
/tools/angsd-0.931/misc/thetaStat do_stat $BASE_DIR$target'.thetas.idx' -win 10000 -step 10000 -outnames $BASE_DIR$target'.thetas.window.idx'

## Do per-chromosome average theta
/tools/angsd-0.931/misc/thetaStat do_stat $BASE_DIR$target'.thetas.idx' -outnames $BASE_DIR$target'.thetas.average.idx'
```

- New method

Using the same sfs files in old method, run the script below (this is the most recent method in Angsd)

```sh
module load angsd/0.931

BASE_DIR='/scratch/hzz0024/DelBay19_Sep/09_theta_new_method/'

for POP in CH REF HC ARN COH SR NB; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 1 '$POP' sfs starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/realSFS \
    $BASE_DIR$POP'_SFS_cv30.saf.idx' \
    -P 16 \
    -fold 1 \
    > $BASE_DIR$POP'.sfs'
done

## Step 2a Estimate theta
for POP in CH REF HC ARN COH SR NB; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 2a '$POP' theta estimate starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/realSFS saf2theta \
    $BASE_DIR$POP'_maf0.05_pctind0.7_cv30.saf.idx' \
    -sfs $BASE_DIR$POP'.sfs' \
    -fold 1 \
    -outname $BASE_DIR$POP
done

## Step 2b Print per-SNP theta
for POP in CH REF HC ARN COH SR NB; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 2a '$POP' per-SNP theta print starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/thetaStat print \
    $BASE_DIR$POP'.thetas.idx' \
    > $BASE_DIR$POP'.thetas.tsv'
done

## Step 3a do fixed window theta
for POP in CH REF HC ARN COH SR NB; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 3a '$POP' fixed window theta starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/thetaStat do_stat \
    $BASE_DIR$POP'.thetas.idx' \
    -win 10000 -step 10000 \
    -outnames $BASE_DIR$POP'.thetas.window.idx'
done

## Step 3b do per-chromosome average theta
for POP in CH REF HC ARN COH SR NB; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 3b '$POP' per-chromosome average theta starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/thetaStat do_stat \
    $BASE_DIR$POP'.thetas.idx' \
    -outnames $BASE_DIR$POP'.thetas.average.idx'
done
```

#### Theta estimates correction

```sh
module load bedtools/2.29.0
WIN=5000 
###in the loop
for CHR in `cat chromosomes.txt`; do ###list of chromosomes in chromosomes.txt file; do
###make bed file for all variant and invariant sites for each chromosome
grep "$CHR" ALL_maf0.05_pctind0.7_cv30_allvar.mafs > ALL_maf0.05_pctind0.7_cv30_allvar_${CHR}.mafs
cut -f 1,2 ALL_maf0.05_pctind0.7_cv30_allvar_${CHR}.mafs > ALL_maf0.05_pctind0.7_cv30_allvar_${CHR}.txt
cut -f2 ALL_maf0.05_pctind0.7_cv30_allvar_${CHR}.txt | awk '{$1 = $1 + 1; print}' | paste ALL_maf0.05_pctind0.7_cv30_allvar_${CHR}.txt - | sed 's/ //g'> ALL_maf0.05_pctind0.7_cv30_allvar_${CHR}.bed

###split the genome window file in chr
grep "$CHR" genome_windows_${WIN}.bed > genome_windows_${WIN}_${CHR}.bed


###calculate the number of sites in each window for each chromosome
bedtools coverage -a genome_windows_${WIN}_${CHR}.bed -b ALL_maf0.05_pctind0.7_cv30_allvar_${CHR}.bed -counts > allvar_${WIN}bwin_${CHR}.txt
## not sure why replace the \t0 here
cut -f4 allvar_${WIN}bwin_${CHR}.txt | sed 's/\t0/NA/g' > allvar_${WIN}win_NA_${CHR}.txt

grep "$CHR" ALL_pi_global.bed > ALL_pi_global_${CHR}.bed
## pate - will add the new column to the end of data
awk '{print exp($4)}' ALL_pi_global_${CHR}.bed | paste ALL_pi_global_${CHR}.bed - > ALL_pi_global_log_${CHR}.bed
# bedtools is used to sum up the theta value in each window
bedtools map -a genome_windows_${WIN}_${CHR}.bed -b ALL_pi_global_log_${CHR}.bed -c 5 -o sum | sed 's/\t[.]/\tNA/g' - > ALL_pi_global_log_${WIN}bwin_${CHR}.txt
###pi_peer_global_noout_log_50kbwin_chr2_pilon.txt

paste ALL_pi_global_log_${WIN}bwin_${CHR}.txt allvar_${WIN}win_NA_${CHR}.txt | sed 's/[.]\t/NA\t/g' - > pi_peer_global_noout_log_${WIN}bwin_sites_${CHR}.txt
## divide the theta by number of SNPs in a window
awk '{if(/NA/)var="NA";else var=$4/$5;print var}' pi_peer_global_noout_log_${WIN}bwin_sites_${CHR}.txt | paste pi_peer_global_noout_log_${WIN}bwin_sites_${CHR}.txt - > pi_peer_global_noout_log_${WIN}bwin_sites_corrected_${CHR}.txt

done
```