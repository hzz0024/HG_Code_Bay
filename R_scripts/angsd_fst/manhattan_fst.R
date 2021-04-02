library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)

plot <- function(fname, outname) {
  #DT <- fread("CHR19_REF19_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_fold.fst")
  DT <- fread(fname)
  print(DT)
  DT$chr <- as.numeric(DT$chr)
  #pdf("Mahattan_ch_ref_1kb_fold.pdf",width=15,height=10)
  par(mfrow=c(1,1)) 
  jpeg(outname, width = 16, height = 9, units = 'in', res = 300)
  manhattan(DT,chr="V1",bp="V2",p="V5",logp=FALSE, cex = 0.2, cex.axis = 0.8, ylim = c(0, 1),
            col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
            ylab="2019 Ref vs 2020 Ref Fst ", cex.lab=1.4) #main = "Chromosome",
  dev.off()
}

plot("CHR19_REF19_minq20_minmq30_1x_CV30_masked_fold.fst", "Mahattan_CHR19_REF19_fold.jpg")
plot("CHR19_CHR20_minq20_minmq30_1x_CV30_masked_fold.fst", "Mahattan_CHR19_CHR20_fold.jpg")
plot("REF19_REF20_minq20_minmq30_1x_CV30_masked_fold.fst", "Mahattan_REF19_REF20_fold.jpg")
