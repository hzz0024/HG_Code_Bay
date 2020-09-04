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


get.dir("REF_maf0.05_pctind0.7_cv30.mafs", "CH_maf0.05_pctind0.7_cv30.mafs", "SR_maf0.05_pctind0.7_cv30.mafs", "HC_maf0.05_pctind0.7_cv30.mafs", "alt_CH_REF_HC_SR.txt")
fish.test("REF_maf0.05_pctind0.7_cv30.mafs", "CH_maf0.05_pctind0.7_cv30.mafs", "alt_CH_REF_HC_SR.txt", "fish_CH_REF_HCSR.txt")
fish.test("SR_maf0.05_pctind0.7_cv30.mafs", "HC_maf0.05_pctind0.7_cv30.mafs", "alt_CH_REF_HC_SR.txt", "fish_HC_SR.txt")

