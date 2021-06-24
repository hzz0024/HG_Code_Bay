#!/bin/sh
Rscript Fish_exact.R --ref1 REF19_minq20_minmq30_1x_CV30_masked.mafs --ch1 CHR19_minq20_minmq30_1x_CV30_masked.mafs --ref2 SR_minq20_minmq30_1x_CV30_masked.mafs --ch2 HC_minq20_minmq30_1x_CV30_masked.mafs --alt alt_REF19_CHR19_SR_HC.txt --out1 fish_SR_HC_REF19_CHR19_SR_HC.txt --out2 fish_REF19_CHR19_REF19_CHR19_SR_HC.txt
Rscript Fish_pcomb.R --file1 fish_REF19_CHR19_REF19_CHR19_SR_HC.txt --file2 fish_SR_HC_REF19_CHR19_SR_HC.txt --global 'REF19_CHR19_SR_HC_'

