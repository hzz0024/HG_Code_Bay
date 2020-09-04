#######################
#  Fisher's exact test#
#######################
# mafs produced by domajorminor 3 and domaf 1
# load reference file with header in it
dat1_name = 'HC_maf0.05_pctind0.7_cv30.mafs'
dat1 <- read.delim(dat1_name, header = TRUE, sep='\t')
dat1_total = dat1$nInd*2
dat1_refCount = round(dat1$knownEM*dat1_total)
dat1_altCount = dat1_total - round(dat1$knownEM*dat1_total)

# load challenge file with header in it
dat2_name = 'COH_maf0.05_pctind0.7_cv30.mafs'
dat2 <- read.delim(dat2_name, header = TRUE, sep='\t')
dat2_total = dat2$nInd*2
dat2_refCount = round(dat2$knownEM*dat2_total)
dat2_altCount = dat2_total - round(dat2$knownEM*dat2_total)

outputfile='fish_HC_COH_pvalue.txt'
sink(outputfile)
idxs <- seq(1, length(dat1_total))
delta_p <- c()
p_values <- c()
alternative <- c()
for(i in idxs){
  #for(i in seq(1,10000)){
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
  #print(M)
  # Calculate delta_p
  p0 <- dat1$knownEM[i]
  p1 <- dat2$knownEM[i]
  delta_p <- p1 - p0
  #print(delta_p)
  #get T/F values for whether change in allele frequency is positive or negative 
  alternative = delta_p <= 0
  alternative[alternative == TRUE] <- 'greater'
  alternative[alternative == FALSE] <- 'less'
  # perform fisher's exact test
  p = fisher.test(M, alternative=alternative)$p.value
  p2 = p*2 ##multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
  p2[p2 > 1] <- 1
  #p_value <- fisher.test(M)$p
  cat(paste0(chr,'\t',pos,'\t',p1,'\t',p0,'\t',delta_p,'\t',p, '\t',p2,'\n'))
}
sink()

#------------------------------------------------------------
# mafs produced by domajorminor 5 and domaf 2
# load reference file with header in it
dat1_name = 'HC_maf0.05_pctind0.7_cv30_anc.mafs'
dat1 <- read.delim(dat1_name, header = TRUE, sep='\t')
dat1_total = dat1$nInd*2
dat1_refCount = round(dat1$unknownEM*dat1_total)
dat1_altCount = dat1_total - round(dat1$unknownEM*dat1_total)

# load challenge file with header in it
dat2_name = 'COH_maf0.05_pctind0.7_cv30_anc.mafs'
dat2 <- read.delim(dat2_name, header = TRUE, sep='\t')
dat2_total = dat2$nInd*2
dat2_refCount = round(dat2$unknownEM*dat2_total)
dat2_altCount = dat2_total - round(dat2$unknownEM*dat2_total)

outputfile='fish_HC_COH_pvalue_anc.txt'
sink(outputfile)
idxs <- seq(1, length(dat1_total))
delta_p <- c()
p_values <- c()
alternative <- c()
for(i in idxs){
  #for(i in seq(1,10000)){
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
  #print(M)
  # Calculate delta_p
  p0 <- dat1$unknownEM[i]
  p1 <- dat2$unknownEM[i]
  delta_p <- p1 - p0
  #print(delta_p)
  #get T/F values for whether change in allele frequency is positive or negative 
  alternative = delta_p <= 0
  alternative[alternative == TRUE] <- 'greater'
  alternative[alternative == FALSE] <- 'less'
  # perform fisher's exact test
  p = fisher.test(M, alternative=alternative)$p.value
  p2 = p*2 ##multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
  p2[p2 > 1] <- 1
  #p_value <- fisher.test(M)$p
  cat(paste0(chr,'\t',pos,'\t',p1,'\t',p0,'\t',delta_p,'\t',p, '\t',p2,'\n'))
}
sink()
