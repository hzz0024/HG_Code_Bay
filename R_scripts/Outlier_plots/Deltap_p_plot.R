##########################
#  eveal the relation-   #
#ship between deltap     #
#and start p, first plot #
##########################

##################### reveal the relationship between deltap and start p, first plot #####################
#setwd("/Volumes/cornell/DelBay19_Hopper/permutation/4_deltap_plot/deltap_vs_p_plot1")
# Delta p plot
library(ggplot2)
# load the neutral datasets
files <- list.files('.', pattern = "*.txt")
# load the dataset
ch_file = 'CH_ref_98_ch_doMAF_filter.mafs.extracted'
ref_file = 'REF_ref_98_ref_doMAF_filter.mafs.extracted'
# obtain the deltap from obs dataset
ch = read.delim(ch_file, header = TRUE, sep = "\t", dec = ".")
ref = read.delim(ref_file, header = TRUE, sep = "\t", dec = ".")
p0 = ref$knownEM
p1 = ch$knownEM
DATA = data.frame(MIN=p0, MAX=p1)
DATA = DATA[order(DATA$MIN),]
num_snp = 3664
DATA$X = seq(1, num_snp)

sp <- ggplot(DATA, aes(x=X, y=MIN)) +
  geom_point(size=.5)+
  geom_point(aes(x=X, y=MAX),size=.5,color='red')+
  # draws the range bars
  geom_errorbar(data=DATA, aes(ymin=MIN, ymax=MAX), width=.001,color='yellow',alpha=.8)
# add x and y-axis titles
sp + scale_x_continuous(name="SNP", limits=c(0, 400)) +
  scale_y_continuous(name="Allele frequency", limits=c(0, 1)) +
  labs(title = "Minor allele changes for 386 SNP outliers",
       subtitle = "ref allele = black, ch allele = red, actual deltap = yellow")

##################### reveal the relationship between deltap and start p, second plot #####################
#setwd("/Volumes/cornell/DelBay19_Hopper/permutation/4_deltap_plot/deltap_vs_p_plot2")
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
  labs(title = "Deltap against reference allele p for the observation data",
       subtitle = "note deltap ranges from 0-0.5")


