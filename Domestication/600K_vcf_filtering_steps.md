# exclude inversions and mtDNA
vcftools --vcf genetyped_data_all_samples.vcf --exclude-bed invers.bed --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --recode --recode-INFO-all --out genetyped_data_all_samples_nochr156invers

# invers.bed list the three big inversions that used during filtering
#chrom chromStart  chromEnd
#1 40630000  42400000
#5 62190000  79000000
#6 30500000  43900000

# extract the sample id from vcf file
bcftools query -l genetyped_data_all_samples_nochr156invers.recode.vcf > sample_id_all_842.txt

# exclude LGF, PCs, VC familes, NYH, and CBW populations
vcftools --vcf genetyped_data_all_samples_nochr156invers.recode.vcf --keep sample_id_subset_514.txt --recode --recode-INFO-all --out genetyped_data_n_514

# check the missing ind
./filter_missing_ind.sh genetyped_data_n_514.recode.vcf genetyped_data_n_514_indmiss

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf genetyped_data_n_514.recode.vcf
	--missing-indv
	--out genetyped_data_n_514_indmiss

After filtering, kept 514 out of 514 Individuals
Outputting Individual Missingness
After filtering, kept 276327 out of a possible 276327 Sites
Run Time = 13.00 seconds



                                          Histogram of % missing data per individual
      450 +---------------------------------------------------------------------------------------------------------+
          |                 +                +                 +                 +                +                 |
          |                                                   'totalmissing' using (bin($1,binwidth)):(1.0) ******* |
      400 |-+               ************************************                                                  +-|
          |                 *                                  *                                                    |
      350 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      300 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
      250 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      200 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
      150 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      100 |-+               *                                  *                                                  +-|
          |                 *                                  ************************************                 |
       50 |-+               *                                  *                                  *               +-|
          |******************                                  *                                  *                 |
          |                 *                +                 *                 +                ******************|
        0 +---------------------------------------------------------------------------------------------------------+
        0.025              0.03            0.035              0.04             0.045             0.05             0.055
                                                       % of missing data

The 85% cutoff would be 0.0413025

# maf and genotype rate filtering
vcftools --vcf genetyped_data_n_514_indmiss.recode.vcf --maf 0.05 --max-missing 0.95 --recode --recode-INFO-all --out genetyped_data_n_514_maf05_maxmiss095

# call rate missing filtering in each population
./pop_missing_filter.sh genetyped_data_n_514_maf05_maxmiss095.recode.vcf popmap.txt 0.95 17 genetyped_data_n_514_maf05_maxmiss095_popmiss095

# filter for HWE
./filter_hwe_by_pop.pl -v genetyped_data_n_514_maf05_maxmiss095_popmiss095.recode.vcf -p popmap.txt -h 0.01 -c 0.5 -o genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe
 
#filter_hwe_by_pop.pl -v <vcffile> -p <popmap> [options]
#
#    Options: -v <vcffile> input vcf file -p <popmap> tab-separated file of
#    samples and population designations -h [hwe] minimum Hardy-Weinberg
#    p-value cutoff for SNPs -c [cutoff] proportion of all populations that a
#    locus can be below HWE cutoff without being filtered -o [out] name of outfile

Processing population: DBW1 (32 inds)
Processing population: DBW2 (32 inds)
Processing population: DBX1 (32 inds)
Processing population: DBX2 (31 inds)
Processing population: DBX3 (31 inds)
Processing population: LIW1 (31 inds)
Processing population: LIW2 (30 inds)
Processing population: MEH2 (36 inds)
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
Kept 179798 of a possible 185937 loci (filtered 6139 loci)

# access the missing rate, call rate and allele frequency distribution
vcftools --vcf genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe.recode.vcf --missing-indv --out genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe
vcftools --vcf genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe.recode.vcf --missing-site --out genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe
vcftools --vcf genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe.recode.vcf --freq2 --max-alleles 2 --out genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe

# code for plink formating
/programs/plink-a2.3LM/plink2 --vcf genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe.recode.vcf --allow-extra-chr --make-bed --out genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe

# extract pruned snp
vcftools --vcf genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe.recode.vcf --snps all_pruned_SNP_list.txt --recode --recode-INFO-all --out genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe_pruned

# summary

Table S1.  Summary of data filtering procedures: rows refer to filtering steps; columns refer to statistics for each step.  For columns, ‘sites’ refers to SNPs, and ‘Inds’ refers to individuals.  ‘Start’, ‘End’, and ‘Removed’ refer, respectively, to the number of each unit before the filtering step, the number after the filtering step, and the number removed with the filter.

|     Filter steps                                                            |     Start sites    |     End sites    |     Start Inds    |     End Inds    |     Removed sites    |     Removed Inds    |
|-----------------------------------------------------------------------------|--------------------|------------------|-------------------|-----------------|----------------------|---------------------|
|     Exclude inversion in Chr 1, 5 and 6, and loci in the mtDNA              |     300446         |     276327       |     842           |     842         |     24119            |     0               |
|     Exclude irrelevant populations                                          |     276327         |     276327       |     842           |     514         |     0                |     328             |
|     Filter_missing_ind script; Ind call rate > 0.9                          |     276327         |     276327       |     514           |     514         |     0                |     0               |
|     Minor allele frequency > 0.05, genotype call rate > 0.95                |     276327         |     147773       |     514           |     514         |     128554           |     0               |
|     pop_missing_filter script; call rate > 0.95 in any single population    |     147773         |     147773       |     514           |     514         |     0                |     0               |
|     Hardy-Weinberg equilibrium                                              |     147773         |     141634       |     514           |     514         |     6139             |     0               |
|     LD clumping using bigsnpr                                               |     141634         |     106205       |     514           |     514         |     35429            |     0               |

