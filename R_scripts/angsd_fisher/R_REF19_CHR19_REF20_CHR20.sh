#!/bin/sh
Rscript Fish_exact.R --ref1 REF19_minq20_minmq30_1x_CV30_masked.mafs --ch1 CHR19_minq20_minmq30_1x_CV30_masked.mafs --ref2 REF20_minq20_minmq30_1x_CV30_masked.mafs --ch2 CHR20_minq20_minmq30_1x_CV30_masked.mafs --alt alt_REF19_CHR19_REF20_CHR20.txt --out1 fish_REF20_CHR20_REF19_CHR19_REF20_CHR20.txt --out2 fish_REF19_CHR19_REF19_CHR19_REF20_CHR20.txt
Rscript Fish_pcomb.R --file1 fish_REF19_CHR19_REF19_CHR19_REF20_CHR20.txt --file2 fish_REF20_CHR20_REF19_CHR19_REF20_CHR20.txt --global 'REF19_CHR19_REF20_CHR20_'

