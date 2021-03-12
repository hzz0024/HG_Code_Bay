
##########################
#  rveal the relation-   #
#ship between deltap     #
#and start p, first plot #
##########################

##################### reveal the relationship between deltap and start p, first plot #####################
setwd("/Volumes/cornell/Fisher_exact")
# Delta p plot
library(ggplot2)
# load the neutral datasets
# files <- list.files('.', pattern = "*.txt")
# load the dataset
ch_file = 'CH_maf0.05_pctind0.7_cv30.mafs.extracted'
ref_file = 'REF_maf0.05_pctind0.7_cv30.mafs.extracted'
# obtain the deltap from obs dataset
ch = read.delim(ch_file, header = FALSE, sep = "\t", dec = ".")
ref = read.delim(ref_file, header = FALSE, sep = "\t", dec = ".")
p0 = ref$V6
p1 = ch$V6
DATA = data.frame(MIN=p0, MAX=p1)
DATA = DATA[order(DATA$MIN),]
num_snp = 96 # change this number according to the outlier file
DATA$X = seq(1, num_snp)

sp <- ggplot(DATA, aes(x=X, y=MIN)) +
  geom_point(size=.5)+
  geom_point(aes(x=X, y=MAX),size=.5,color='red')+
  # draws the range bars
  geom_errorbar(data=DATA, aes(ymin=MIN, ymax=MAX), width=.001,color='yellow',alpha=.8)
# add x and y-axis titles
sp + scale_x_continuous(name="SNP", limits=c(0, 100)) +
  scale_y_continuous(name="Allele frequency", limits=c(0, 1)) +
  labs(title = "Allele frequency changes for 96 potential outliers identified from Fisher's exact test",
       subtitle = "ref allele = black, ch allele = red, actual deltap = yellow")

##################### reveal the relationship between deltap and start p, second plot #####################
library(ggplot2)
# load the dataset
ch_file = 'CH_maf0.05_pctind0.7_cv30.mafs.extracted'
ref_file = 'REF_maf0.05_pctind0.7_cv30.mafs.extracted'
#deltap_file = 'obs_deltap.output'
# obtain the deltap from obs dataset
ch = read.delim(ch_file, header = FALSE, sep = "\t", dec = ".")
ref = read.delim(ref_file, header = FALSE, sep = "\t", dec = ".")
deltap = ch$V6 - ref$V6
p0 = ref$V6
p1 = ch$V6
DATA = data.frame(p=p0, delta_p=deltap)
DATA = DATA[order(DATA$p),]
num_snp = 96
sp <- ggplot(DATA, aes(x=p, y=delta_p)) +
  geom_point(size=.5)
# add x and y-axis titles
sp + scale_x_continuous(name="p", limits=c(min(DATA$p), max(DATA$p))) +
  scale_y_continuous(name="Deltap (absolute values)", limits=c(min(DATA$delta_p), max(DATA$delta_p))) +
  labs(title = "Deltap against reference allele p for Fisher's exact result (96 outliers after FDR correction)",
       subtitle = "note deltap ranges from 0-0.5")

## check the delta_p pattern

library(ggplot2)
ch_file = 'SR_maf0.05_pctind0.7_cv30.mafs'
ref_file = 'NB_maf0.05_pctind0.7_cv30.mafs'
# obtain the deltap from obs dataset
ch = read.delim(ch_file, header = TRUE, sep = "\t", dec = ".")
ref = read.delim(ref_file, header = TRUE, sep = "\t", dec = ".")
p0 = ref$knownEM
p1 = ch$knownEM

#filename='/Users/ryan/Downloads/softwareEstMAF_20110510/test/out1/test_tabgeno_fl0.0'
#kimdat = read.delim(filename, header=FALSE,sep=' ', dec='.')
#kimdat = read.delim(text = gsub("\t", " ", readLines(filename)),header=FALSE,sep=' ', dec='.')
#p0 = kimdat$V2
#p1 = kimdat$V3
deltap = p1-p0

p1 = p1[abs(deltap)> 0.2]
p0 = p0[abs(deltap)> 0.2]
deltap = deltap[abs(deltap) > 0.2]

DATA = data.frame(p=p0, delta_p=deltap)

DATA = DATA[order(DATA$p),]
num_snp = length(p1)
sp <- ggplot(DATA, aes(x=p, y=delta_p)) +
  geom_point(size=.5)
# add x and y-axis titles
sp + scale_x_continuous(name="p", limits=c(min(DATA$p), max(DATA$p))) +
  scale_y_continuous(name="Deltap", limits=c(min(DATA$delta_p), max(DATA$delta_p))) +
  labs(title = "Deltap against reference allele p for Fisher's exact result (96 outliers after FDR correction)",
       subtitle = "note deltap ranges from 0-0.5")

d=data.frame(x=p0,y=p1-p0)
p <- ggplot(d, aes(x=x, y=y)) + geom_violin()

