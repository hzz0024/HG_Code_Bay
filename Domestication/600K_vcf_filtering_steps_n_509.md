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

# check the missing ind
source /programs/miniconda3/bin/activate dDocent-2.8.13
./filter_missing_ind.sh genetyped_data_n_509.recode.vcf genetyped_data_n_509_indmiss

# maf and genotype rate filtering
vcftools --vcf genetyped_data_n_509_indmiss.recode.vcf --maf 0.05 --max-missing 0.95 --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095

# start population level filtering
# replace the origianl sample id with more clear sample strata (prepared by Excel)
# extract the original vcf sample id
bcftools query -l genetyped_data_n_509_maf05_maxmiss095.recode.vcf > sample_original_name
bcftools reheader -s sample_rename genetyped_data_n_509_maf05_maxmiss095.recode.vcf > genetyped_data_n_509_maf05_maxmiss095.rename.vcf

# call rate missing filtering in each population
./pop_missing_filter.sh genetyped_data_n_509_maf05_maxmiss095.rename.vcf popmap.txt 0.95 17 genetyped_data_n_509_maf05_maxmiss095_popmiss095

# filter for HWE

#filter_hwe_by_pop.pl -v <vcffile> -p <popmap> [options]
#
#    Options: -v <vcffile> input vcf file -p <popmap> tab-separated file of
#    samples and population designations -h [hwe] minimum Hardy-Weinberg
#    p-value cutoff for SNPs -c [cutoff] proportion of all populations that a
#    locus can be below HWE cutoff without being filtered -o [out] name of outfile

./filter_hwe_by_pop.pl -v genetyped_data_n_509_maf05_maxmiss095_popmiss095.recode.vcf -p popmap.txt -h 0.01 -c 0.5 -o genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe

# evaluate the missing rate, call rate and allele frequency distribution
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --missing-indv --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --missing-site --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --freq2 --max-alleles 2 --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe

# code for plink formating
/programs/plink-a2.3LM/plink2 --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --allow-extra-chr --make-bed --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe

# start LD clumping process, see PCAdapt code
# extract pruned snp
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --snps all_pruned_SNP_list.txt --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned
/programs/plink-a2.3LM/plink2 --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned.recode.vcf --allow-extra-chr --make-bed --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned


# Table 1. Summary of data filtering procedures: rows refer to filtering steps; columns refer to statistics for each step.  For columns, ‘sites’ refers to SNPs, and ‘Inds’ refers to individuals.  
# ‘Start’, ‘End’, and ‘Removed’ refer, respectively, to the number of each unit before the filtering step, the number after the filtering step, and the number removed with the filter.

|                             Filter   steps                             | Start sites | End sites | Start Inds | End Inds | Removed sites | Removed Inds |
|:----------------------------------------------------------------------:|:-----------:|:---------:|:----------:|:--------:|:-------------:|:------------:|
|      Exclude inversion in Chr 1, 5   and 6, and loci in the mtDNA      |    300446   |   276327  |     842    |    842   |     24119     |       0      |
|                     Exclude irrelevant populations                     |    276327   |   276327  |     842    |    509   |       0       |      333     |
|            Filter_missing_ind script; Ind   call rate > 0.9            |    276327   |   276327  |     509    |    509   |       0       |       0      |
|       Minor allele frequency >   0.05, genotype call rate > 0.95       |    276327   |   148063  |     509    |    509   |     128264    |       0      |
| pop_missing_filter script;   call rate > 0.95 in any single population |    148063   |   148063  |     509    |    509   |       0       |       0      |
|                       Hardy-Weinberg equilibrium                       |    148063   |   141960  |     509    |    509   |      6103     |       0      |
|                       LD clumping using bigsnpr                        |    141960   |   106456  |     509    |    509   |     35504     |       0      |

