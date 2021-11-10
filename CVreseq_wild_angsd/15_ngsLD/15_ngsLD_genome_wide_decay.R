#########################################
####   process mafs file for Angsd   ####
#########################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/no_shared/format_snp_list")
pname = "All_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked_noinvers.snplist.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
dat = dat[with(dat, order(V1, V2)),]
angsd_list <- paste0(dat$V1, ":" , dat$V2, "-", dat$V2)
write.table(angsd_list, paste0(strsplit(pname, split = ".txt")[[1]][1], ".rf.txt"), row.names=F, col.names = F, quote=F, sep="\t")

# bash script
# for i in {0..9}; do
# j=$[i+1]
# grep 'NC_03578'$i'.1' All_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked_noinvers.snplist.rf.txt > 'All_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked_noinvers.snplist.rf.chr'$j'.txt'
# done

####################################
###   format the ngsLD output    ###
####################################
setwd("~/Dropbox/Mac/Documents/HG/CVreseq_wild_angsd/15_ngsLD")
#ngsLD outputs a TSV file with LD results for all pairs of sites for which LD was calculated, where the first two columns are positions of the SNPs, the third column is the distance (in bp) between the SNPs, and the following 4 columns are the various measures of LD calculated (r^2 from pearson correlation between expected genotypes, D from EM algorithm, D' from EM algorithm, and r^2 from EM algorithm). 

rm(list=ls())
format<-function(ngsLD_output, win){
  #ngsLD_output = "./CLP.chr00.output"
  dat1 = read.delim(ngsLD_output, header = FALSE, sep='\t')
  results = dat1[with(dat1, order(V1, V2)),]
  results <- results[complete.cases(results), ]
  colnames(results)=c('Snp1', 'Snp2', 'Distance', "R1", "D", "Dp", "R2")
  ### Plot the Decay of LD with Distance
  #plot(results$Distance, y=results$R2, xlab="Distance between SNPs (bp)", ylab=expression(R^2), pch=20, col=rgb(0,0,0,alpha=0.2), main="Decay of Linkage Disequilibrium")
  ### Add a smoothed moving average line
  kb.bins=seq(0,max(results$Distance), by=win) # Get all numbers between 1 and the maximum Distance value, in increments of win = 500 for example
  r2=rep(0, length(kb.bins))
  kb.midpt=rep(0, length(kb.bins))# Set up an empty vector to hold the mean r-squared for each "bin"
  LD.averages=data.frame(kb.bins, r2, kb.midpt) # Create the results table (to hold calculations)
  
  ### Go through all of the "bin" values,
  ### For each one, find the subset of data that
  ### falls in that bin, and get the mean r-squared
  for (i in 1:length(kb.bins)) {
    data.interval=subset(results, (results$Distance >= kb.bins[i] & results$Distance < (kb.bins[i]+win))) # Get the subset of data that falls in the ith bin
    LD.averages$r2[i]=mean(data.interval$R2) # Calculate the average R-squared value for this data subset, and save in my "LD.averages" results table
    LD.averages$kb.midpt[i]=(kb.bins[i]+kb.bins[i]+win)/2000
  }
  return(LD.averages)
}

####################################
### save as rdata format (list)  ###
####################################
out <- list()
chr <- seq(0,9)
pop <- c("CLP","CSD","HCD", "HCV")
for (c in chr) {
  for (p in pop){
    name0 = paste0(p, ".chr", c, ".output")
    out[[name0]] = format(paste0(p, ".chr0", c, ".output"), 500)
  }
}

save(out, file = "./output/ngsLD.RData")

####################################
###   start plot for r2 pattern  ###
####################################
library(export)
setwd("~/Dropbox/Mac/Documents/HG/CVreseq_wild_angsd/15_ngsLD")
# We generated mean r2 within 500 bp bins of distances using r-code.
# 3 columns: kb.bins, r2, kb.midpt

rm(list=ls())
library(ggplot2) # cut_interval()
load("./output/ngsLD.RData")
# check the rdata
#LDrdata <- get(load('data/LDanalysis.500bpBINS.rdata'))
pop <- substr(names(out),1,3)
meta <- read.csv('./output/meta_pops_V1.csv', stringsAsFactors=T) # added by @HG
site <- meta$Site.Abb[match(pop,meta$Site.Abb)]

#tiff("output/LDanalysis.makeCurves.500bpBINS.jpg", units="in", width=16, height=12, res=300)
#pdf('output/LDanalysis.makeCurves.500bpBINS.pdf',width=8,height=5)
#quartz(width=8,height=5)
par(mfrow=c(2,2),mar=c(8,2,2,1)) # increase the first two will show x and y-axis labels, last one is controling the width
stats <- c()
cbPalette <- c("#1BA3C6", "#33A65C", "#F8B620", "#E03426", "#EB73B3", "#AEC7E8", "#FF7F0E", "#9EDAE5", "#FFBB78")
for (i in 1:4)
{
  plot(1,1,xlim=c(0,5),ylim=c(0,0.4),xlab="distance (kbp)", ylab="LD (r^2)",type="n")
  tmp <- out[site==levels(site)[i]]
  
  for(j in 1:length(tmp))
  {
    xbar <- tmp[[j]]
    Cstart <- c(C=0.1)
    CDist <- function(n,C,distance)
    {
      ((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance)))
    }
    n=10 ###### change the value to corresponding sample size
    xbar = xbar[which(xbar$r2 != "Inf"),] ###### to get rid of inf r2 values
    modelC = try(nls(r2 ~ CDist(n,C,kb.midpt), data=xbar, start=Cstart, control=nls.control(maxiter=100)))#
    #error = function(e) modelC=print("oops"))
    if(length(modelC)>1)
    {
      xbar$prd <- predict(modelC, newdata = xbar) #@ newdata = xbar added by HG following the https://stackoverflow.com/questions/33309792/r-predict-function-returning-too-many-values 
      halfdecay = (max(xbar$prd))*0.5
      halfdecaydist <- xbar$kb.midpt[which.min(abs(xbar$prd-halfdecay))]
      dist.LD.is.10perc <- xbar$kb.midpt[which.min(abs(xbar$prd-0.1))] 
      stats <- rbind(stats,data.frame(pop=substr(names(tmp)[j],1,3),
                                      chr=substr(names(tmp)[j],5,9),
                                      halfdecay,
                                      halfdecaydist,
                                      dist.LD.is.10perc))
      lines(xbar$kb.midpt, xbar$prd, col=alpha(cbPalette[i],.5), lwd=1)
    }
    else {stats <- rbind(stats,data.frame(pop=substr(names(tmp)[j],1,3),
                                          chr=substr(names(tmp)[j],5,9),
                                          halfdecaydist=NA,
                                          halfdecaydist=NA,
                                          dist.LD.is.10perc=NA))}
  }
}

write.table(stats,"output/ngsLD.stats.500bpBINS.csv",sep=",",row.names=F,quote = F)

graph2ppt(file="output/LD_50k.pptx", width=12, height=12)

#dev.off()

################################################
###   Wilcoxon rank sum test for r2 pattern  ###
################################################

CH1 = stats$dist.LD.is.10perc[which(stats$pop == "CH1")]
RE1 = stats$dist.LD.is.10perc[which(stats$pop == "RE1")]
CH2 = stats$dist.LD.is.10perc[which(stats$pop == "CH2")]
RE2 = stats$dist.LD.is.10perc[which(stats$pop == "RE2")]
DHC = stats$dist.LD.is.10perc[which(stats$pop == "DHC")]
DNB = stats$dist.LD.is.10perc[which(stats$pop == "DNB")]
DSR = stats$dist.LD.is.10perc[which(stats$pop == "DSR")]
wilcox.test(DHC, DNB)
wilcox.test(DHC, DSR)
wilcox.test(CH2, RE2)

