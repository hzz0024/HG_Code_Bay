##########################
#  reveal the relation-  #
#ship between deltap     #
#and start p, first plot #
##########################
library(ggplot2)
setwd("~/Documents/Ryan_workplace/DelBay/Selection_fish_1/combined_ps")

out_name = "REF-CH-SR-HC_out_0.1.txt"
outlier <- read.delim(out_name, header = FALSE, sep='\t')
#c <- make_id(out_name)

setwd("~/Documents/Ryan_workplace/DelBay/Selection_fish_1/combined_ps/sep_ps")
pv_file1 = 'fish_REF_CH_REF-CH-SR-HC.txt'
pv_file2 = 'fish_SR_HC_REF-CH-SR-HC.txt'
A = read.delim(pv_file1, header = FALSE, sep='\t')
B = read.delim(pv_file2, header = FALSE, sep='\t')

names = paste0(A$V1,'_',A$V2)
A_pv = A[names %in% a,]$V7
B_pv = B[names %in% a,]$V7

# shared outliers and output as table
#tmp1 = A[names %in% c,]
#tmp2 = B[names %in% c,]
#cmp_table = data.frame(tmp1$V1, tmp1$V2, tmp1$V7, tmp2$V7 )
#names(cmp_table) = c('chr','pos','REF_CH','SR_HC')
#write.table(cmp_table, 'REF-CH-SR-HC_pv.txt', quote=F, col.names = T, row.names = F, sep='\t')

##################################################################
# reveal the relationship between deltap and start p, first plot #
##################################################################

setwd("~/Documents/Ryan_workplace/DelBay/Selection_fish_1/combined_ps")

out_name = "REF-CH-NB-HC_out_0.1.txt"
outlier <- read.delim(out_name, header = FALSE, sep='\t')

file1 = 'HC_maf0.05_pctind0.7_cv30.mafs' #p1
file2 = 'NB_maf0.05_pctind0.7_cv30.mafs' #p0
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat2 = read.delim(file2, header = TRUE, sep='\t')

file3 =  'CH_maf0.05_pctind0.7_cv30.mafs' #p1
file4 =  'REF_maf0.05_pctind0.7_cv30.mafs' #p0
dat3 = read.delim(file3, header = TRUE, sep='\t')
dat4 = read.delim(file4, header = TRUE, sep='\t')

#REF1 = dat2 #p0
#CH1 = dat1 #p1
#ALL = outlier
pop_name = "NB_HC"
pop0 = "NB"

######## test if randomly picked SNP has a sum delta_p > 0.3 or not
a = dat1$knownEM - dat2$knownEM
b = dat3$knownEM - dat4$knownEM
sum( (abs(b) + abs(a)) > 0.3)
sum( (abs(b) + abs(a)) > 0.3)/length(a)

sum( (abs(deltap) + abs(deltap2)) > 0.3)

sample_id = sample(seq(1,length(p0)),11)
sum((abs(b[sample_id]) + abs(a[sample_id])) > 0.3)

##################################################################
#                      start first plot                          #
##################################################################
#==================fixed=======================
p0 = dat2$knownEM
p1 = dat1$knownEM
names = paste0(dat2$chromo,'_',dat2$position)
target = outlier
target_names = paste0(target$V1,'_',target$V2)
#target_names = c
p0 = p0[names %in% target_names]
p1 = p1[names %in% target_names]
deltap = p1-p0 

p3 = dat4$knownEM
p4 = dat3$knownEM
p3 = p3[names %in% target_names]
p4 = p4[names %in% target_names]
deltap2 = p4 - p3 

colors = vector(mode="character", length=length(deltap))
filter_idx1 =  abs(deltap) < 0.1 
filter_idx2 =  abs(deltap) >= 0.1 & abs(deltap) < 0.2
filter_idx3 =  abs(deltap) >= 0.2
colors[filter_idx1] = 'blue'
colors[filter_idx2] = 'black'
colors[filter_idx3] = 'gray'

DATA = data.frame(X=p0, Y=deltap)
sp <- ggplot(DATA, aes(x=X, y=Y)) +
  geom_point(size=1.5, color = colors) 
# add x and y-axis titles
sp + scale_x_continuous(name="p0", limits=c(min(DATA$X), max(DATA$X))) +
  scale_y_continuous(name="Delta_p", limits=c(min(DATA$Y), max(DATA$Y))) +
  labs(title = paste0(pop_name, " delta_p against p0, p0 is allele frequency in population ", pop0)) + 
  #geom_abline(slope = -1, intercept = 0) +
  #geom_text(aes(x = 0.05, y = 0, label = "y= -x", color = "red")) +
  theme(legend.position="none")
# save as jpg
ggsave(paste0(pop_name, '_fisher1.jpg'), width = 20, height = 16, units = "cm")
print(sp)
dev.off()

DATA = data.frame(X=p3, Y=deltap2)
#DATA = DATA[order(DATA$X),]
sp <- ggplot(DATA, aes(x=X, y=Y)) +
  geom_point(size=1.5, color = colors) 
# add x and y-axis titles
sp + scale_x_continuous(name="p0", limits=c(min(DATA$X), max(DATA$X))) +
  scale_y_continuous(name="Delta_p", limits=c(min(DATA$Y), max(DATA$Y))) +
  labs(title = paste0('CH_REF', " delta_p against p0, p0 is allele frequency in population ", "REF")) + 
  #geom_abline(slope = -1, intercept = 0) +
  #geom_text(aes(x = 0.05, y = 0, label = "y= -x", color = "red")) +
  theme(legend.position="none")
# save as jpg
ggsave(paste0(pop_name, '_fisher2.jpg'), width = 20, height = 16, units = "cm")
print(sp)
dev.off()

DATA0 = DATA
#================sample========================
p0 = dat2$knownEM
p1 = dat1$knownEM
names = paste0(dat2$chromo,'_',dat2$position)
target_names = sample(names, length(target_names), replace = FALSE)
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
  #geom_abline(slope = -2, intercept = 0.1) +
  #geom_text(aes(x = 0.18, y = -0.375, label = "y= -2x + 0.1", color = "red")) + 
  #geom_abline(slope = -2, intercept = 1) +
  #geom_text(aes(x = 0.68, y = -0.25, label = "y= -2x + 1", color = "red")) + 
  geom_abline(slope = -1, intercept = 0) +
  geom_text(aes(x = 0.35, y = -0.375, label = "y= -x", color = "red")) +
  theme(legend.position="none")

# save as jpg
ggsave(paste0(pop_name, '_fdr01_random_sample.jpg'), width = 20, height = 16, units = "cm")
print(sp)
dev.off()

##################### reveal the relationship between deltap and start p, second plot #####################
#REF1 = dat2 # REF1 is p0
#CH1 = dat1 # CH1 is p1
#ALL = outlier
pop_name = "NB_HC"

# obtain the deltap from obs dataset
p0 = dat2$knownEM
p1 = dat1$knownEM
names = paste0(dat2$chromo,'_',dat2$position)
target = outlier
target_names = paste0(target$V1,'_',target$V2)
#target_names = c
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
       subtitle = paste0("p0 = black, p1 = red, delta_p = yellow, n = ", length(target$V3)))
# save as jpg
ggsave(paste0(pop_name, '_fdr1_p0_p1.jpg'), width = 20, height = 16, units = "cm")
print(sp)
dev.off()

# determine the population name
pop_name = "REF_CH"
# obtain the deltap from obs dataset
p0 = dat4$knownEM
p1 = dat3$knownEM
names = paste0(dat4$chromo,'_',dat4$position)
target = outlier
target_names = paste0(target$V1,'_',target$V2)
#target_names = c
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
       subtitle = paste0("p0 = black, p1 = red, delta_p = yellow, n = ", length(target$V3)))
# save as jpg
ggsave(paste0(pop_name, '_fdr1_p0_p1.jpg'), width = 20, height = 16, units = "cm")
print(sp)
dev.off()

#================random sample========================

p0 = dat2$knownEM
p1 = dat1$knownEM
names = paste0(dat2$chromo,'_',dat1$position)
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
ggsave(paste0(pop_name, '_fdr01_p0_p1_random_sample.jpg'), width = 20, height = 16, units = "cm")
print(sp)
dev.off()



################################################## Below code is test purpose, useful to check the delta_p pattern

#names1 = paste0(dat1$chromo,'_',dat1$position)
#names2 = paste0(dat2$chromo,'_',dat2$position)
#datt1 = dat1[names1 %in% c,]
#datt2 = dat2[names2 %in% c,]

# load the maf files
#setwd("/Volumes/cornell/Fisher_exact")
#CH_file = 'CH_maf0.05_pctind0.7_cv30.mafs'
#REF_file = 'REF_maf0.05_pctind0.7_cv30.mafs'

#CH = read.delim(CH_file, header = TRUE, sep = "\t", dec = ".")
#REF = read.delim(REF_file, header = TRUE, sep = "\t", dec = ".")

# load the outlier files
#setwd("/Volumes/cornell/Fisher_exact/results")
#CH_REF_name = 'Fish_CH_REF_fdr1.txt'
#CH_REF <- read.delim(CH_REF_name, header = TRUE, sep=' ')

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