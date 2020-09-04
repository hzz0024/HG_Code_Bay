##########################
#  reveal the relation-   #
#ship between deltap     #
#and start p, first plot #
##########################
library(ggplot2)
# load the maf files
setwd("/Volumes/cornell/Fisher_exact")
CH_file = 'CH_maf0.05_pctind0.7_cv30.mafs'
REF_file = 'REF_maf0.05_pctind0.7_cv30.mafs'
HC_file = 'HC_maf0.05_pctind0.7_cv30.mafs'
ARN_file = 'ARN_maf0.05_pctind0.7_cv30.mafs'
COH_file = 'COH_maf0.05_pctind0.7_cv30.mafs'
SR_file = 'SR_maf0.05_pctind0.7_cv30.mafs'
NB_file = 'NB_maf0.05_pctind0.7_cv30.mafs'

CH = read.delim(CH_file, header = TRUE, sep = "\t", dec = ".")
REF = read.delim(REF_file, header = TRUE, sep = "\t", dec = ".")
HC = read.delim(HC_file, header = TRUE, sep = "\t", dec = ".")
ARN = read.delim(ARN_file, header = TRUE, sep = "\t", dec = ".")
COH = read.delim(COH_file, header = TRUE, sep = "\t", dec = ".")
SR = read.delim(SR_file, header = TRUE, sep = "\t", dec = ".")
NB = read.delim(NB_file, header = TRUE, sep = "\t", dec = ".")

# load the outlier files
setwd("/Volumes/cornell/Fisher_exact/results")
CH_REF_name = 'Fish_CH_REF_fdr2.txt'
CH_REF <- read.delim(CH_REF_name, header = TRUE, sep=' ')
HC_ARN_name = 'Fish_HC_ARN_fdr2.txt'
HC_ARN <- read.delim(HC_ARN_name, header = TRUE, sep=' ')
HC_COH_name = 'Fish_HC_COH_fdr2.txt'
HC_COH <- read.delim(HC_COH_name, header = TRUE, sep=' ')
HC_SR_name = 'Fish_HC_SR_fdr2.txt'
HC_SR <- read.delim(HC_SR_name, header = TRUE, sep=' ')
HC_NB_name = 'Fish_HC_NB_fdr2.txt'
HC_NB <- read.delim(HC_NB_name, header = TRUE, sep=' ')
ARN_COH_name = 'Fish_ARN_COH_fdr2.txt'
ARN_COH <- read.delim(ARN_COH_name, header = TRUE, sep=' ')
ARN_SR_name = 'Fish_ARN_SR_fdr2.txt'
ARN_SR <- read.delim(ARN_SR_name, header = TRUE, sep=' ')
ARN_NB_name = 'Fish_ARN_NB_fdr2.txt'
ARN_NB <- read.delim(ARN_NB_name, header = TRUE, sep=' ')
COH_SR_name = 'Fish_COH_SR_fdr2.txt'
COH_SR <- read.delim(COH_SR_name, header = TRUE, sep=' ')
COH_NB_name = 'Fish_COH_NB_fdr2.txt'
COH_NB <- read.delim(COH_NB_name, header = TRUE, sep=' ')
SR_NB_name = 'Fish_SR_NB_fdr2.txt'
SR_NB <- read.delim(SR_NB_name, header = TRUE, sep=' ')

##################### reveal the relationship between deltap and start p, first plot #####################

FIRST = REF
SECOND = CH
ALL = CH_REF
pop_name = "CH-REF"
pop0 = "REF"

#==================fixed=======================
p0 = FIRST$knownEM
p1 = SECOND$knownEM
names = paste0(FIRST$chromo,'_',SECOND$position)
target = ALL
target_names = paste0(target$chr,'_',target$pos)
p0 = p0[names %in% target_names]
p1 = p1[names %in% target_names]
deltap = p1-p0 
DATA = data.frame(X=p0, Y=deltap)
DATA = DATA[order(DATA$X),]
sp <- ggplot(DATA, aes(x=X, y=Y)) +
  geom_point(size=.5)
# add x and y-axis titles
sp + scale_x_continuous(name="p0", limits=c(min(DATA$X), max(DATA$X))) +
  scale_y_continuous(name="Delta_p", limits=c(min(DATA$Y), max(DATA$Y))) +
  labs(title = paste0(pop_name, " delta_p against p0, p0 is allele frequency in population ", pop0)) + 
  geom_abline(slope = -1, intercept = 0) +
  geom_text(aes(x = 0.05, y = 0, label = "y= -x", color = "red")) +
  theme(legend.position="none")
# save as jpg
ggsave(paste0(pop_name, '_fdr02.jpg'), width = 20, height = 16, units = "cm")
print(sp)
dev.off()

DATA0 = DATA

#================sample========================
p0 = FIRST$knownEM
p1 = SECOND$knownEM
names = paste0(FIRST$chromo,'_',SECOND$position)
target_names =sample(names, length(target_names), replace = FALSE)
p0 = p0[names %in% target_names]
p1 = p1[names %in% target_names]
DATA = data.frame(X=p0, Y=p1-p0)
#DATA = DATA[order(DATA$X),]
sp <- ggplot(DATA, aes(x=X, y=Y)) +
  geom_point(size=.5)
# add x and y-axis titles
sp + scale_x_continuous(name="p0", limits=c(min(DATA0$X), max(DATA0$X))) +
  scale_y_continuous(name="Delta_p", limits=c(min(DATA0$Y), max(DATA0$Y))) +
  labs(title = paste0(pop_name, " random sampling delta_p against p0, p0 is allele frequency in population ", pop0)) + 
  geom_abline(slope = -2, intercept = 0.1) +
  geom_text(aes(x = 0.18, y = -0.375, label = "y= -2x + 0.1", color = "red")) + 
  geom_abline(slope = -2, intercept = 1) +
  geom_text(aes(x = 0.68, y = -0.25, label = "y= -2x + 1", color = "red")) + 
  geom_abline(slope = -1, intercept = 0) +
  geom_text(aes(x = 0.35, y = -0.375, label = "y= -x", color = "red")) +
  theme(legend.position="none")

# save as jpg
ggsave(paste0(pop_name, '_fdr02_random_sample.jpg'), width = 20, height = 16, units = "cm")
print(sp)
dev.off()

##################### reveal the relationship between deltap and start p, second plot #####################
FIRST = REF # FIRST is p0
SECOND = CH # SECOND is p1
ALL = CH_REF
pop_name = "CH_REF"

# obtain the deltap from obs dataset
p0 = FIRST$knownEM
p1 = SECOND$knownEM
names = paste0(FIRST$chromo,'_',SECOND$position)
target = ALL
target_names = paste0(target$chr,'_',target$pos)
p0 = p0[names %in% target_names]
p1 = p1[names %in% target_names]
deltap = p1-p0 
DATA = data.frame(MIN=p0, MAX=p1)
DATA = DATA[order(DATA$MIN),]
num_snp = length(deltap) # change this number according to the outlier file
DATA$X = seq(1, num_snp)

sp <- ggplot(DATA, aes(x=X, y=MIN)) +
  geom_point(size=.5)+
  geom_point(aes(x=X, y=MAX),size=.5,color='red')+
  # draws the range bars
  geom_errorbar(data=DATA, aes(ymin=MIN, ymax=MAX), width=.001,color='yellow',alpha=.8)
# add x and y-axis titles
sp + scale_x_continuous(name="SNPs sorted based on p0", limits=c(0, num_snp)) +
  scale_y_continuous(name="Allele frequency", limits=c(0, 1)) +
  labs(title = paste0("Outlier p1 relative to p0 in population pair ", pop_name),
       subtitle = "p0 = black, p1 = red, delta_p = yellow")

# save as jpg
ggsave(paste0(pop_name, '_fdr02_p0_p1.jpg'), width = 20, height = 16, units = "cm")
print(sp)
dev.off()


#================sample========================

p0 = FIRST$knownEM
p1 = SECOND$knownEM
names = paste0(FIRST$chromo,'_',SECOND$position)
target_names =sample(names, length(target_names), replace = FALSE)
p0 = p0[names %in% target_names]
p1 = p1[names %in% target_names]
deltap = p1-p0 
DATA = data.frame(MIN=p0, MAX=p1)
DATA = DATA[order(DATA$MIN),]
num_snp = length(deltap) # change this number according to the outlier file
DATA$X = seq(1, num_snp)

sp <- ggplot(DATA, aes(x=X, y=MIN)) +
  geom_point(size=.5)+
  geom_point(aes(x=X, y=MAX),size=.5,color='red')+
  # draws the range bars
  geom_errorbar(data=DATA, aes(ymin=MIN, ymax=MAX), width=.001,color='yellow',alpha=.8)
# add x and y-axis titles
sp + scale_x_continuous(name="SNPs sorted based on p0", limits=c(0, num_snp)) +
  scale_y_continuous(name="Allele frequency", limits=c(0, 1)) +
  labs(title = paste0("Random SNPs p1 relative to p0 in population pair ", pop_name),
       subtitle = "p0 = black, p1 = red, delta_p = yellow")

# save as jpg
ggsave(paste0(pop_name, '_fdr02_p0_p1_random_sample.jpg'), width = 20, height = 16, units = "cm")
print(sp)
dev.off()
################################################## Below code is test purpose, useful to check the delta_p pattern

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

