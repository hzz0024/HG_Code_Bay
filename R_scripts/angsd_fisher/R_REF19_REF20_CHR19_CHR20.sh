#!/bin/sh
Rscript Fish_exact.R --ref1 REF19_minq20_minmq30_1x_CV30_masked.mafs --ch1 REF20_minq20_minmq30_1x_CV30_masked.mafs --ref2 CHR19_minq20_minmq30_1x_CV30_masked.mafs --ch2 CHR20_minq20_minmq30_1x_CV30_masked.mafs --alt alt_REF19_REF20_CHR19_CHR20.txt --out1 fish_CHR19_CHR20_REF19_REF20_CHR19_CHR20.txt --out2 fish_REF19_REF20_REF19_REF20_CHR19_CHR20.txt
Rscript Fish_pcomb.R --file1 fish_REF19_REF20_REF19_REF20_CHR19_CHR20.txt --file2 fish_CHR19_CHR20_REF19_REF20_CHR19_CHR20.txt --global 'REF19_REF20_CHR19_CHR20_'

