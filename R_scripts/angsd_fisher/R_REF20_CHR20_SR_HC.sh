#!/bin/sh
Rscript Fish_exact.R --ref1 REF20_minq20_minmq30_1x_CV30_masked.mafs --ch1 CHR20_minq20_minmq30_1x_CV30_masked.mafs --ref2 SR_minq20_minmq30_1x_CV30_masked.mafs --ch2 HC_minq20_minmq30_1x_CV30_masked.mafs --alt alt_REF20_CHR20_SR_HC.txt --out1 fish_SR_HC_REF20_CHR20_SR_HC.txt --out2 fish_REF20_CHR20_REF20_CHR20_SR_HC.txt
Rscript Fish_pcomb.R --file1 fish_REF20_CHR20_REF20_CHR20_SR_HC.txt --file2 fish_SR_HC_REF20_CHR20_SR_HC.txt --global 'REF20_CHR20_SR_HC_'

