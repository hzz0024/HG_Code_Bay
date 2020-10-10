### ngsLD

1) extract the SNPs from snplist and assign them into different chromosomes with script *extract_chr.sh*

output example:

Number of SNPs in chr1 is:
273241
Number of SNPs in chr2 is:
294774
Number of SNPs in chr3 is:
315491
Number of SNPs in chr4 is:
283496
Number of SNPs in chr5 is:
408592
Number of SNPs in chr6 is:
14784
Number of SNPs in chr7 is:
75924
Number of SNPs in chr8 is:
107832
Number of SNPs in chr9 is:
132474
Number of SNPs in chr10 is:
27430

2) Using script *wt_script_05_ALL.sh* to create beagle files for ngsLD calculation

note change the output folder and location of 01_config.sh for script running

3) Generate the ngsLD output using *extract.sh*

4) 

Rscript --vanilla --slave fit_LDdecay.R --ld_files LD.list --out ALL_chr5_k25_no_edit.jpg --fit_level 0 --max_kb_dist 500 --fit_level 10
Rscript --vanilla --slave fit_LDdecay.R --ld_files LD.list --out ALL_chr5_k25_no_edit1.jpg --fit_level 0 --max_kb_dist 500 --fit_level 10 --plot_size 1.5,2
Rscript --vanilla --slave fit_LDdecay.R --ld_files LD.list --out ALL_chr5_k25_no_edit2.jpg --fit_level 0 --max_kb_dist 500 --fit_level 10 --plot_size 2,2 --plot_axis_scales free_y