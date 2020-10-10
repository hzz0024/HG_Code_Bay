#!/bin/bash
# Read a string with spaces using for loop
for pop in CHR
do
    #echo -e 'module load angsd/0.931\n###this script will work on bamfiles by population and calculate saf & maf\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /workdir/hz269/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 -minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth 384 -b $'$pop' -out "/workdir/hz269/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask"\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -setMaxDepth 384 -b $'$pop' -out "/workdir/hz269/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask_allvar"\n#main features\n# -P nb of threads\n# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2)\n# -doMajorMinor 1 use the most frequent allele as major\n# -anc provide a ancestral sequence = reference in our case\n# -rf (file with the region written) work on a defined region : OPTIONAL\n# -b (bamlist) input file\n# -out  output file\n\n#main filters\n#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ minimum quality of reads\n#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%\n#filter on allele frequency -minMaf, set to 0.05\n\n#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles\n#order the sites file by chromosome names\n#makes a region file matching the sites files and with same order\n#index sites file\necho "from the maf file, extract a list of SNP chr, positoin, major all, minor all"\ncd /workdir/hz269/04_saf_maf_gl_all\nzcat "$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask.mafs.gz | tail -n +2 > FILE_mask.tmp && mv FILE_mask.tmp "$target"_snplist_mask\nawk '\''{print $1,$2,$3,$4}'\'' "$target"_snplist_mask > "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_mask\nangsd sites index "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_mask' >> formal/'04_saf_maf_al_all_'$pop'_mask.sh'
    echo -e 'module load angsd/0.931\n###this script will work on bamfiles by population and calculate saf & maf\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /workdir/hz269/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 -minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth 384 -b $'$pop' -out "/workdir/hz269/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -setMaxDepth 384 -b $'$pop' -out "/workdir/hz269/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30_allvar"\n#main features\n# -P nb of threads\n# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2)\n# -doMajorMinor 1 use the most frequent allele as major\n# -anc provide a ancestral sequence = reference in our case\n# -rf (file with the region written) work on a defined region : OPTIONAL\n# -b (bamlist) input file\n# -out  output file\n\n#main filters\n#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ minimum quality of reads\n#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%\n#filter on allele frequency -minMaf, set to 0.05\n\n#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles\n#order the sites file by chromosome names\n#makes a region file matching the sites files and with same order\n#index sites file\necho "from the maf file, extract a list of SNP chr, positoin, major all, minor all"\ncd /workdir/hz269/04_saf_maf_gl_all\nzcat "$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30.mafs.gz | tail -n +2 > FILE_cv30.tmp && mv FILE_cv30.tmp "$target"_snplist_cv30\nawk '\''{print $1,$2,$3,$4}'\'' "$target"_snplist_cv30 > "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_cv30\nangsd sites index "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_cv30' >> formal/'04_saf_maf_al_all_'$pop'_cv30.sh'
done
