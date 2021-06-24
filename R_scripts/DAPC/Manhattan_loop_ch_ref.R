# manhattan plot of Fst
# install.packages("qqman")
# install.packages("caret")
install.packages("animation")

library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)

require(data.table)

##### script for single-SNP plot (due to difficulty in opening the pdf, I export the jepg plot here)
DT <- fread("ch_ref_doSAF_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#pdf("Mahattan_ch_ref_1kb_fold.pdf",width=15,height=10)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs REF Fst ", cex.lab=1.4) #main = "Chromosome",
dev.off()

DT <- fread("ch_ref_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#pdf("Mahattan_ch_ref_1kb_fold.pdf",width=15,height=10)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_1kb_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs REF 1kb Window Fst ", cex.lab=1.4) #main = "Chromosome",
dev.off()

DT <- fread("ch_ref_5kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_5kb_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs REF 5kb Window Fst ", cex.lab=1.4) #main = "Chromosome",
dev.off()

DT <- fread("ch_ref_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_15kb_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs REF 15kb Window Fst ", cex.lab=1.4) #main = "Chromosome",
dev.off()

#Necessary lines for all plot scripts below
LG.chr = scan("LG.list",what=" ") # empty quotes indicates character string
DT <- fread("ch_ref_no16inv_minI8D8maxD32.fst")
DT_complete <- DT[complete.cases(DT),]
testDT<-DT_complete
##### Leo script for wide plots, loop to make 10 individual LARGE pdf files
##### I couldn't figure out how to combine them onto one page, but they are too large to combine anyway

##### make 10 individual square Fst manhattan plots, png format
i=0
LG.chr = scan("LG.list",what=" ") # empty quotes indicates character string
for (LG in LG.chr){
  i <- i+1
  currDT <- testDT[testDT$chr == i,]  #using the bracket "]" notation which designates the indices of the data set. The first index is for the rows and the second for the columns. leaving the index for the columns blank indicates that we want currDT to contain all the variables (columns) of the original data frame.
  currDT$chr <- as.character(i)  # replace with counter number
  currDT$chr <- as.numeric(currDT$chr)  #convert back to numeric
  
  options(device=png)   # plot to screen (used to work)
  jpeg(paste("Mahattan_ch_ref_no16inv_minI8D8maxD32_fold_", i, ".jpg", sep = ""), width = 16, height = 9, units = 'in', res = 300)
  manhattan(currDT,chr="chr",bp="pos",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 1),
            col=c("#3D3D3D","#B0B0B0"),genomewideline=F, suggestiveline=F,
            main = LG, ylab="ch_ref angsd Fst")  
  #plotname <- sprintf("Mahattan_ch_ref_no16inv_minI8D8maxD32_unfold_%d.png", i) # sprintf = string print function with %s for text and %d for digit
  #dev.copy(png, filename=plotname, width = 16, height = 9, units = 'in', res = 300)  #copy the contents of the graph window to a file without having to re-enter the commands.
  dev.off()
}

