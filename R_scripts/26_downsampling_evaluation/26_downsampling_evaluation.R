setwd("~/Dropbox/Mac/Documents/HG/DelBay20_adult/26_downsampling_evaluation")
name1 = "CHR20_all_minq20_minmq30_CV30_masked.original.mafshhhgggg"
dat1 = read.delim(name1, header = TRUE, sep='\t')

name2 = "CHR20_all_minq20_minmq30_CV30_masked.downsampling.mafs"
dat2 = read.delim(name2, header = TRUE, sep='\t')

cor(dat1$knownEM, dat2$knownEM, method = c("pearson"))
cor.test(dat1$knownEM, dat2$knownEM, method = c("pearson"))
