Steps of GO enrichment analysis 

1. submit the protein sequences to http://eggnog-mapper.embl.de/ for annotation and mapping
2. Download the protein table from NCBI (a total of 60201 proteins) 
   https://www.ncbi.nlm.nih.gov/genome/398 > click tabular in "Download genome annotation"
3. Using R script 17_adding_LOC_to_eggnog.R to format the eggnog mapper result - "out.emapper.annotations.txt". Remember to delete the headers with # and process information at the bottom. However, uncomment the column header for R processing.
4. Build the GO and KEGG library: Rscript 17_create_orgdb_from_emapper.R out.emapper.annotation.withLOC.uniq.txt
5. Do enrichment test