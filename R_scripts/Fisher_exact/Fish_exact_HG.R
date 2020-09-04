#######################
#  Fisher's exact test#
#######################
##############################################################
# Fisher's exact test for domajorminor 3 and domaf 1 outputs #
##############################################################
get.dir <- function(ref_name1, ch_name1, ref_name2, ch_name2, output_name){
  dat1 <- read.delim(ref_name1, header = TRUE, sep='\t')
  dat2 <- read.delim(ch_name1, header = TRUE, sep='\t')
  dat3 <- read.delim(ref_name2, header = TRUE, sep='\t')
  dat4 <- read.delim(ch_name2, header = TRUE, sep='\t')
  outputfile = output_name
  sink(outputfile)
  idxs <- seq(1, length(dat1$knownEM))
  for(i in idxs){
  #for(i in seq(1,10000)){
    # print the running process
    s = paste0(i,'/',dim(dat1)[1])
    message(s,"\r",appendLF=FALSE)
    chr <- dat1$chromo[i]
    pos <- dat1$position[i]
    # Calculate delta_p in replicate study 1
    p1 <- dat1$knownEM[i]
    p2 <- dat2$knownEM[i]
    dp1 <- p2 - p1
    # Calculate delta_p in replicate study 2
    p3 <- dat3$knownEM[i]
    p4 <- dat4$knownEM[i]
    dp2 <- p4 - p3
    dp_mean = (dp1+dp2)/2
    #get T/F values for whether change in allele frequency is positive or negative
    dp_mean.greater = dp_mean <= 0
    alternative = dp_mean.greater
    alternative[alternative == TRUE] <- 'greater'
    alternative[alternative == FALSE] <- 'less'
    cat(paste0(chr,'\t',pos,'\t',dp1,'\t',dp2,'\t',dp_mean,'\t',alternative,'\n'))
  }
  sink()
}

fish.test <- function(ref_name, ch_name, alt_name, output_name){
  dat1 <- read.delim(ref_name, header = TRUE, sep='\t')
  dat1_total = dat1$nInd*2
  dat1_refCount = round(dat1$knownEM*dat1_total)
  dat1_altCount = dat1_total - round(dat1$knownEM*dat1_total)

  # load challenge file with header in it
  dat2 <- read.delim(ch_name, header = TRUE, sep='\t')
  dat2_total = dat2$nInd*2
  dat2_refCount = round(dat2$knownEM*dat2_total)
  dat2_altCount = dat2_total - round(dat2$knownEM*dat2_total)
  
  # load the alternative for Fisher one-side test
  alt <- read.delim(alt_name, header = FALSE, sep='\t')
  alt_vector <- alt$V6
  
  outputfile = output_name
  sink(outputfile)
  idxs <- seq(1, length(dat1_total))
  for(i in idxs){
    # print the running process
    s = paste0(i,'/',dim(dat1)[1])
    message(s,"\r",appendLF=FALSE)
    # create the matrix for Fisher's exact test
    chr <- dat1$chromo[i]
    pos <- dat1$position[i]
    c.a1 = dat1_refCount[i]
    s.a1 = dat2_refCount[i]
    c.a2 = dat1_altCount[i]
    s.a2 = dat2_altCount[i]
    st = dat2_total[i]
    ct = dat1_total[i]
    survival = c(s.a1, s.a2)
    control = c(c.a1, c.a2)
    M = as.table(cbind(control, survival))
    rownames(M) = c("ref","alt")
    colnames(M) = c("control", "survival")
    # Calculate delta_p
    p0 <- dat1$knownEM[i]
    p1 <- dat2$knownEM[i]
    delta_p <- p1 - p0
    #print(delta_p)
    alternative= toString(alt_vector[i])
    # perform fisher's exact test
    if(class(alternative) != 'character'){
      message('ERROR!NOT STRING!')
      message(alternative)
    } 
    p = fisher.test(M, alternative=alternative)$p.value
    p2 = p*2 ##multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
    p2[p2 > 1] <- 1
    #p_value <- fisher.test(M)$p
    cat(paste0(chr,'\t',pos,'\t',p1,'\t',p0,'\t',delta_p,'\t',p, '\t',p2,'\n'))
  }
  sink()
#  return(M)
}


get.dir("REF_maf0.05_pctind0.7_cv30.mafs", "CH_maf0.05_pctind0.7_cv30.mafs",
        "COH_maf0.05_pctind0.7_cv30.mafs", "HC_maf0.05_pctind0.7_cv30.mafs", "alt_CH_REF_HC_COH.txt")
fish.test("REF_maf0.05_pctind0.7_cv30.mafs", "CH_maf0.05_pctind0.7_cv30.mafs", "alt_CH_REF_HC_COH.txt", "fish_CH_REF.txt")
fish.test("COH_maf0.05_pctind0.7_cv30.mafs", "HC_maf0.05_pctind0.7_cv30.mafs", "alt_CH_REF_HC_COH.txt", "fish_HC_COH.txt")

##############################################################
# Fisher's exact test for domajorminor 5 and domaf 2 outputs #
##############################################################
get_anc.dir <- function(ref_name1, ch_name1, ref_name2, ch_name2, output_name){
  dat1 <- read.delim(ref_name1, header = TRUE, sep='\t')
  dat2 <- read.delim(ch_name1, header = TRUE, sep='\t')
  dat3 <- read.delim(ref_name2, header = TRUE, sep='\t')
  dat4 <- read.delim(ch_name2, header = TRUE, sep='\t')
  outputfile = output_name
  sink(outputfile)
  idxs <- seq(1, length(dat1$unknownEM))
  for(i in idxs){
    #for(i in seq(1,10000)){
    # print the running process
    s = paste0(i,'/',dim(dat1)[1])
    message(s,"\r",appendLF=FALSE)
    chr <- dat1$chromo[i]
    pos <- dat1$position[i]
    # Calculate delta_p in replicate study 1
    p1 <- dat1$unknownEM[i]
    p2 <- dat2$unknownEM[i]
    dp1 <- p2 - p1
    # Calculate delta_p in replicate study 2
    p3 <- dat3$unknownEM[i]
    p4 <- dat4$unknownEM[i]
    dp2 <- p4 - p3
    dp_mean = (dp1+dp2)/2
    #get T/F values for whether change in allele frequency is positive or negative
    dp_mean.greater = dp_mean <= 0
    alternative = dp_mean.greater
    alternative[alternative == TRUE] <- 'greater'
    alternative[alternative == FALSE] <- 'less'
    cat(paste0(chr,'\t',pos,'\t',dp1,'\t',dp2,'\t',dp_mean,'\t',alternative,'\n'))
  }
  sink()
}

fish_anc.test <- function(ref_name, ch_name, alt_name, output_name){
  dat1 <- read.delim(ref_name, header = TRUE, sep='\t')
  dat1_total = dat1$nInd*2
  dat1_refCount = round(dat1$unknownEM*dat1_total)
  dat1_altCount = dat1_total - round(dat1$unknownEM*dat1_total)
  
  # load challenge file with header in it
  dat2 <- read.delim(ch_name, header = TRUE, sep='\t')
  dat2_total = dat2$nInd*2
  dat2_refCount = round(dat2$unknownEM*dat2_total)
  dat2_altCount = dat2_total - round(dat2$unknownEM*dat2_total)
  
  # load the alternative for Fisher one-side test
  alt <- read.delim(alt_name, header = FALSE, sep='\t')
  alt_vector <- alt$V6
  
  outputfile = output_name
  sink(outputfile)
  idxs <- seq(1, length(dat1_total))
  for(i in idxs){
    # print the running process
    s = paste0(i,'/',dim(dat1)[1])
    message(s,"\r",appendLF=FALSE)
    # create the matrix for Fisher's exact test
    chr <- dat1$chromo[i]
    pos <- dat1$position[i]
    c.a1 = dat1_refCount[i]
    s.a1 = dat2_refCount[i]
    c.a2 = dat1_altCount[i]
    s.a2 = dat2_altCount[i]
    st = dat2_total[i]
    ct = dat1_total[i]
    survival = c(s.a1, s.a2)
    control = c(c.a1, c.a2)
    M = as.table(cbind(control, survival))
    rownames(M) = c("ref","alt")
    colnames(M) = c("control", "survival")
    # Calculate delta_p
    p0 <- dat1$unknownEM[i]
    p1 <- dat2$unknownEM[i]
    delta_p <- p1 - p0
    #print(delta_p)
    alternative= alt_vector[i]
    # perform fisher's exact test
    if(class(alternative) != 'character'){
      message('ERROR!NOT STRING!')
    } 
    p = fisher.test(M, alternative=alternative)$p.value
    p2 = p*2 ##multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
    p2[p2 > 1] <- 1
    #p_value <- fisher.test(M)$p
    cat(paste0(chr,'\t',pos,'\t',p1,'\t',p0,'\t',delta_p,'\t',p, '\t',p2,'\n'))
  }
  sink()
  #  return(M)
}


get_anc.dir("REF_maf0.05_pctind0.7_cv30_anc.mafs", "CH_maf0.05_pctind0.7_cv30_anc.mafs",
        "COH_maf0.05_pctind0.7_cv30_anc.mafs", "HC_maf0.05_pctind0.7_cv30_anc.mafs", "alt_CH_REF_HC_COH_anc.txt")
fish_anc.test("REF_maf0.05_pctind0.7_cv30_anc.mafs", "CH_maf0.05_pctind0.7_cv30_anc.mafs", "alt_CH_REF_HC_COH_anc.txt", "fish_CH_REF_anc.txt")
fish_anc.test("COH_maf0.05_pctind0.7_cv30_anc.mafs", "HC_maf0.05_pctind0.7_cv30_anc.mafs", "alt_CH_REF_HC_COH_anc.txt", "fish_HC_COH_anc.txt")


#######################
#  Adjust the p-value #
#######################

library(qvalue)
write_table <- function(res_name, output_fdr2, output_fdr1, output_fdr5){
  #res_name = deparse(substitute(dat))
  res <- read.delim(res_name, header = FALSE, sep='\t')
  p1 <- res$V3
  p0 <- res$V4
  p_value <- res$V7
  length(which(p_value<0.05))
  hist(p_value)
  qobj <- qvalue(p = p_value)
  fdr <- qobj$qvalues
  length(fdr[fdr<0.2])
  outlier <- as.data.frame(cbind(res$V1[fdr<0.2], res$V2[fdr<0.2], res$V4[fdr<0.2], res$V5[fdr<0.2], p_value[fdr<0.2], fdr[fdr<0.2]))
  colnames(outlier) <- c("chr","pos", "p0", "delta_p", "p_value","fdr")
  # export the dataframe to txt file
  write.table(outlier,output_fdr2, quote = FALSE, row.names = FALSE)
  length(fdr[fdr<0.1])
  outlier <- as.data.frame(cbind(res$V1[fdr<0.1], res$V2[fdr<0.1], res$V4[fdr<0.1], res$V5[fdr<0.1], p_value[fdr<0.1], fdr[fdr<0.1]))
  colnames(outlier) <- c("chr","pos", "p0", "delta_p", "p_value","fdr")
  # export the dataframe to txt file
  write.table(outlier,output_fdr1, quote = FALSE, row.names = FALSE)
  length(fdr[fdr<0.05])
  outlier <- as.data.frame(cbind(res$V1[fdr<0.05], res$V2[fdr<0.05], res$V4[fdr<0.05], res$V5[fdr<0.05], p_value[fdr<0.05], fdr[fdr<0.05]))
  colnames(outlier) <- c("chr","pos", "p0", "delta_p", "p_value","fdr")
  # export the dataframe to txt file
  write.table(outlier,output_fdr5, quote = FALSE, row.names = FALSE)
  print('Done...')
}
#------------------------------------------------------------
# result by domajorminor 3 and domaf 1
write_table("fish_CH_REF_pvalue.txt", "./results/Fish_CH_REF_fdr2.txt" , "./results/Fish_CH_REF_fdr1.txt" , "./results/Fish_CH_REF_fdr5.txt")
write_table("fish_HC_COH_pvalue.txt", "./results/Fish_HC_COH_fdr2.txt" , "./results/Fish_HC_COH_fdr1.txt" , "./results/Fish_HC_COH_fdr5.txt")
write_table("fish_HC_NB_pvalue.txt", "./results/Fish_HC_NB_fdr2.txt" , "./results/Fish_HC_NB_fdr1.txt" , "./results/Fish_HC_NB_fdr5.txt")
# result by domajorminor 5 and domaf 2
write_table("fish_CH_REF_pvalue_anc.txt", "./results/Fish_CH_REF_fdr2_anc.txt" , "./results/Fish_CH_REF_fdr1_anc.txt" , "./results/Fish_CH_REF_fdr5_anc.txt")
write_table("fish_HC_COH_pvalue_anc.txt", "./results/Fish_HC_COH_fdr2_anc.txt" , "./results/Fish_HC_COH_fdr1_anc.txt" , "./results/Fish_HC_COH_fdr5_anc.txt")
write_table("fish_HC_NB_pvalue_anc.txt", "./results/Fish_HC_NB_fdr2_anc.txt" , "./results/Fish_HC_NB_fdr1_anc.txt" , "./results/Fish_HC_NB_fdr5_anc.txt")

##########################
#  reveal the relation-  #
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
CH_REF_name = 'Fish_CH_REF_fdr1.txt'
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
ggsave(paste0(pop_name, '_fdr01.jpg'), width = 20, height = 16, units = "cm")
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
ggsave(paste0(pop_name, '_fdr01_p0_p1.jpg'), width = 20, height = 16, units = "cm")
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
ggsave(paste0(pop_name, '_fdr01_p0_p1_random_sample.jpg'), width = 20, height = 16, units = "cm")
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

