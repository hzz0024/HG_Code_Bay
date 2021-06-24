#!/bin/sh
Rscript Fish_exact.R --ref1 REF19_minq20_minmq30_1x_CV30_masked.mafs --ch1 REF20_minq20_minmq30_1x_CV30_masked.mafs --ref2 ARN_minq20_minmq30_1x_CV30_masked.mafs --ch2 COH_minq20_minmq30_1x_CV30_masked.mafs --alt alt_REF19_REF20_ARN_COH.txt --out1 fish_ARN_COH_REF19_REF20_ARN_COH.txt --out2 fish_REF19_REF20_REF19_REF20_ARN_COH.txt
Rscript Fish_pcomb.R --file1 fish_REF19_REF20_REF19_REF20_ARN_COH.txt --file2 fish_ARN_COH_REF19_REF20_ARN_COH.txt --global 'REF19_REF20_ARN_COH_'

