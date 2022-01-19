# exclude inversions and mtDNA
vcftools --vcf genetyped_data_all_samples.vcf --exclude-bed invers.bed --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --recode --recode-INFO-all --out genetyped_data_all_samples_nochr156invers

# invers.bed list the three big inversions that used during filtering
#chrom chromStart  chromEnd
#1 40630000  42400000
#5 62190000  79000000
#6 30500000  43900000

# extract the sample id from vcf file
bcftools query -l genetyped_data_all_samples_nochr156invers.recode.vcf > sample_id_all_842.txt

# exclude LGF, PCs, VC familes, NYH, and CBW populations, and exclude one individual (DBW1-30) that potentailly mislabeled
vcftools --vcf genetyped_data_all_samples_nochr156invers.recode.vcf --keep sample_id_509.txt --recode --recode-INFO-all --out genetyped_data_n_509

VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_all_samples_nochr156invers.recode.vcf
  --keep sample_id_509.txt
  --recode-INFO-all
  --out genetyped_data_n_509
  --recode

Keeping individuals in 'keep' list
After filtering, kept 509 out of 842 Individuals
Outputting VCF file...
After filtering, kept 276327 out of a possible 276327 Sites
Run Time = 77.00 seconds

# check the missing ind
source /programs/miniconda3/bin/activate dDocent-2.8.13
./filter_missing_ind.sh genetyped_data_n_509.recode.vcf genetyped_data_n_509_indmiss

VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_all_samples_nochr156invers.recode.vcf
  --keep sample_id_509.txt
  --recode-INFO-all
  --out genetyped_data_n_509
  --recode

Keeping individuals in 'keep' list
After filtering, kept 509 out of 842 Individuals
Outputting VCF file...
After filtering, kept 276327 out of a possible 276327 Sites
Run Time = 76.00 seconds
[hz269@cbsuhare maf005]$ source /programs/miniconda3/bin/activate dDocent-2.8.13
(dDocent-2.8.13) [hz269@cbsuhare maf005]$ ./filter_missing_ind.sh genetyped_data_n_509.recode.vcf genetyped_data_n_509_indmiss

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_n_509.recode.vcf
  --missing-indv
  --out genetyped_data_n_509_indmiss

After filtering, kept 509 out of 509 Individuals
Outputting Individual Missingness
After filtering, kept 276327 out of a possible 276327 Sites
Run Time = 13.00 seconds



                                          Histogram of % missing data per individual
      400 +---------------------------------------------------------------------------------------------------------+
          |                 *                +                 *                 +                +                 |
          |                 *                                 '*otalmissing' using (bin($1,binwidth)):(1.0) ******* |
      350 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      300 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      250 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
      200 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      150 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      100 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  ************************************                 |
       50 |-+               *                                  *                                  *               +-|
          |******************                                  *                                  *                 |
          |                 *                +                 *                 +                ******************|
        0 +---------------------------------------------------------------------------------------------------------+
        0.025              0.03            0.035              0.04             0.045             0.05             0.055
                                                       % of missing data

The 85% cutoff would be 0.0411867
Would you like to set a different cutoff, yes or no
yes
Please enter new cutoff
0.1
All individuals with more than 10.0% missing data will be removed.

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_n_509.recode.vcf
  --remove lowDP.indv
  --recode-INFO-all
  --out genetyped_data_n_509_indmiss
  --recode

Excluding individuals in 'exclude' list
After filtering, kept 509 out of 509 Individuals
Outputting VCF file...
After filtering, kept 276327 out of a possible 276327 Sites

# maf and genotype rate filtering
vcftools --vcf genetyped_data_n_509_indmiss.recode.vcf --maf 0.05 --max-missing 0.95 --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_n_509_indmiss.recode.vcf
  --recode-INFO-all
  --maf 0.05
  --max-missing 0.95
  --out genetyped_data_n_509_maf05_maxmiss095
  --recode

After filtering, kept 509 out of 509 Individuals
Outputting VCF file...
After filtering, kept 148063 out of a possible 276327 Sites
Run Time = 40.00 seconds

# start population level filtering
# replace the origianl sample id with more clear sample strata (prepared by Excel)
# extract the original vcf sample id
bcftools query -l genetyped_data_n_509_maf05_maxmiss095.recode.vcf > sample_original_name
bcftools reheader -s sample_rename genetyped_data_n_509_maf05_maxmiss095.recode.vcf > genetyped_data_n_509_maf05_maxmiss095.rename.vcf

# call rate missing filtering in each population
./pop_missing_filter.sh genetyped_data_n_509_maf05_maxmiss095.rename.vcf popmap.txt 0.95 17 genetyped_data_n_509_maf05_maxmiss095_popmiss095

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_n_509_maf05_maxmiss095.rename.vcf
  --exclude-positions loci.to.remove
  --recode-INFO-all
  --out genetyped_data_n_509_maf05_maxmiss095_popmiss095
  --recode

After filtering, kept 509 out of 509 Individuals
Outputting VCF file...
After filtering, kept 148063 out of a possible 148063 Sites
Run Time = 30.00 seconds

# filter for HWE
./filter_hwe_by_pop.pl -v genetyped_data_n_509_maf05_maxmiss095_popmiss095.recode.vcf -p popmap.txt -h 0.01 -c 0.5 -o genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe
 
#filter_hwe_by_pop.pl -v <vcffile> -p <popmap> [options]
#
#    Options: -v <vcffile> input vcf file -p <popmap> tab-separated file of
#    samples and population designations -h [hwe] minimum Hardy-Weinberg
#    p-value cutoff for SNPs -c [cutoff] proportion of all populations that a
#    locus can be below HWE cutoff without being filtered -o [out] name of outfile

Processing population: DBW1 (31 inds)
Processing population: DBW2 (32 inds)
Processing population: DBX1 (32 inds)
Processing population: DBX2 (31 inds)
Processing population: DBX3 (31 inds)
Processing population: LIW1 (31 inds)
Processing population: LIW2 (30 inds)
Processing population: MEH2 (32 inds)
Processing population: MEW1 (30 inds)
Processing population: MEW2 (31 inds)
Processing population: NCW1 (32 inds)
Processing population: NCW2 (30 inds)
Processing population: NEH1 (32 inds)
Processing population: NEH2 (32 inds)
Processing population: UMFS (30 inds)
Processing population: UNC1 (20 inds)
Processing population: UNC2 (22 inds)
Outputting results of HWE test for filtered loci to 'filtered.hwe'
Kept 141960 of a possible 148063 loci (filtered 6103 loci)

# access the missing rate, call rate and allele frequency distribution
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --missing-indv --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --missing-site --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --freq2 --max-alleles 2 --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe

# code for plink formating
/programs/plink-a2.3LM/plink2 --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --allow-extra-chr --make-bed --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe

# start LD clumping process, see PCAdapt code
# extract pruned snp
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --snps all_pruned_SNP_list.txt --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned
/programs/plink-a2.3LM/plink2 --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned.recode.vcf --allow-extra-chr --make-bed --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned


# extract neutral snp. used for Ne estimate
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --exclude-bed Outflank_outliers_q05_n_1366.bed --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outFLANK
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outFLANK.recode.vcf --exclude-bed PCAdapt_outliers_q05_n_428.bed --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outlier

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf
  --recode-INFO-all
  --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outFLANK
  --recode
  --exclude-bed Outflank_outliers_q05_n_1366.bed

After filtering, kept 509 out of 509 Individuals
Outputting VCF file...
  Read 1366 BED file entries.
After filtering, kept 140595 out of a possible 141960 Sites
Run Time = 30.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outFLANK.recode.vcf
  --recode-INFO-all
  --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outlier
  --recode
  --exclude-bed PCAdapt_outliers_q05_n_428.bed

After filtering, kept 509 out of 509 Individuals
Outputting VCF file...
  Read 428 BED file entries.
After filtering, kept 140354 out of a possible 140595 Sites
Run Time = 29.00 seconds

# extract neutral SNPs for Ne
bcftools query -f '%CHROM\t%POS\t%ID\n' genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outlier.recode.vcf > neutral_SNPs_n_140354.list
shuf -n 500 neutral_SNPs_n_140376.list | cut -f3 -d$'\t' > neutral_SNPs_n_500_list.txt
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outlier.recode.vcf --snps neutral_SNPs_n_500_list.txt --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_500

bcftools query -f '%CHROM\t%POS\t%ID\n' genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outlier.recode.vcf > neutral_SNPs_n_140376.list
shuf -n 1000 neutral_SNPs_n_140376.list | cut -f3 -d$'\t' > neutral_SNPs_n_1K_list.txt
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outlier.recode.vcf --snps neutral_SNPs_n_1K_list.txt --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_1K

bcftools query -f '%CHROM\t%POS\t%ID\n' genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outlier.recode.vcf > neutral_SNPs_n_140376.list
shuf -n 5000 neutral_SNPs_n_140376.list | cut -f3 -d$'\t' > neutral_SNPs_n_5K_list.txt
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_no_outlier.recode.vcf --snps neutral_SNPs_n_5K_list.txt --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_5K

# format with perl script vcf2genepop.pl
./vcf2genepop_hg.pl vcf=genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_500.recode.vcf pops=MEW1,MEW2,UMFS,MEH2,LIW1,LIW2,NEH1,NEH2,DBW1,DBW2,DBX1,DBX2,DBX3,NCW1,NCW2,UNC1,UNC2 > genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_500.gen
./vcf2genepop_hg.pl vcf=genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_1K.recode.vcf pops=MEW1,MEW2,UMFS,MEH2,LIW1,LIW2,NEH1,NEH2,DBW1,DBW2,DBX1,DBX2,DBX3,NCW1,NCW2,UNC1,UNC2 > genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_1K.gen
./vcf2genepop_hg.pl vcf=genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_5K.recode.vcf pops=MEW1,MEW2,UMFS,MEH2,LIW1,LIW2,NEH1,NEH2,DBW1,DBW2,DBX1,DBX2,DBX3,NCW1,NCW2,UNC1,UNC2 > genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_5K.gen

# extract thinned neutral SNPs. Used for STRUCTURE analysis
./vcf.sh
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned.recode.vcf --exclude-bed Outflank_outliers_q05_n_1366.bed --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outFLANK
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outFLANK.recode.vcf --exclude-bed PCAdapt_outliers_q05_n_428.bed --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outlier

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned.recode.vcf
  --recode-INFO-all
  --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outFLANK
  --recode
  --exclude-bed Outflank_outliers_q05_n_1366.bed

After filtering, kept 509 out of 509 Individuals
Outputting VCF file...
  Read 1366 BED file entries.
After filtering, kept 106044 out of a possible 106456 Sites

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outFLANK.recode.vcf
  --recode-INFO-all
  --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outlier
  --recode
  --exclude-bed PCAdapt_outliers_q05_n_428.bed

After filtering, kept 509 out of 509 Individuals
Outputting VCF file...
  Read 428 BED file entries.
After filtering, kept 105933 out of a possible 106044 Sites
Run Time = 22.00 seconds

# extract 10K random neutral SNPs
bcftools query -f '%CHROM\t%POS\t%ID\n' genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outlier.recode.vcf > neutral_SNPs_n_105933.list
shuf -n 5000 neutral_SNPs_n_105933.list | cut -f3 -d$'\t' > neutral_SNPs_n_5K_list.txt
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outlier.recode.vcf --snps neutral_SNPs_n_5K_list.txt --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_5K

#format for structure format
# extract the original vcf sample id
bcftools query -l genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_5K.recode.vcf
# format with perl script vcf2genepop.pl
./vcf2genepop_hg.pl vcf=genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_5K.recode.vcf pops=MEW1,MEW2,UMFS,MEH2,LIW1,LIW2,NEH1,NEH2,DBW1,DBW2,DBX1,DBX2,DBX3,NCW1,NCW2,UNC1,UNC2 > genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_10K.gen


The plink.hom file has the following format, one row per identified homozygous region:
     FID      Family ID
     IID      Individual ID
     CHR      Chromosome
     SNP1     SNP at start of region
     SNP2     SNP at end of region
     POS1     Physical position (bp) of SNP1
     POS2     Physical position (bp) of SNP2
     KB       Length of region (kb)
     NSNP     Number of SNPs in run
     DENSITY  Average SNP density (1 SNP per kb)
     PHOM     Proportion of sites homozygous
     PHET     Proportion of sites heterozygous


Het - heterzogisty across SNPs (it is not individual based het)
Nanimals - number of samples in each popualtion
NROH - sum of ROH across individuals
KBROH - mean of ROH across individuals
dummy_length - dummy ROH length assuming all SNPs are homozygous across the genome
FROH - Total Froh as total kb of ROH divided by dummy ROH length (maximal detectable ROH using current settings)
FROH1_2, FROH2_4, FROH4_8, FROH8_16, FROH16 - ROH per individual when only considering runs between 1000 kb and 2000 kb etc.  

Fhat are inbreeding coefficient from Plink ibc command
Fhat1 = Diagonal of GRM (based on variance in additive genetic values): not completely true
Fhat2 = Plink inbreeding based on expected heterozygosity (excess heterozygotes)
Fhat3 = Correlation between uniting gametes
