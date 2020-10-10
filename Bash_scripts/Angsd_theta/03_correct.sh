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