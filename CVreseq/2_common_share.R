# check common shared SNP
file1 = 'CS_NEH_outlier_SNP.txt' 
h1 = read.delim(file1, header = FALSE, sep='\t')

file2 = 'CS_DEBY_outlier_SNP.txt' 
h2 = read.delim(file2, header = FALSE, sep='\t')

file3 = 'SL_OBOYS2_outlier_SNP.txt' 
h3 = read.delim(file3, header = FALSE, sep='\t')

intersect(h1$V1, h2$V1)

intersect(h3$V1, h1$V1)

intersect(h3$V1, h2$V1)

intersect(h3$V1, intersect(h1$V1, h2$V1))
