#######################
#  Fisher's exact test#
#######################
##############################################################
# Fisher's exact test for domajorminor 3 and domaf 1 outputs #
##############################################################
library("Exact")
fish.test <- function(ref_name, ch_name, output_name){
  #ref_name = "REF19_minq20_minmq30_1x_CV30_masked.mafs"
  dat1 <- read.delim(ref_name, header = TRUE, sep='\t')
  dat1_total = dat1$nInd*2
  dat1_refCount = round(dat1$knownEM*dat1_total)
  dat1_altCount = dat1_total - round(dat1$knownEM*dat1_total)
  
  # load challenge file with header in it
  #ch_name = "CHR19_minq20_minmq30_1x_CV30_masked.mafs"
  dat2 <- read.delim(ch_name, header = TRUE, sep='\t')
  dat2_total = dat2$nInd*2
  dat2_refCount = round(dat2$knownEM*dat2_total)
  dat2_altCount = dat2_total - round(dat2$knownEM*dat2_total)
  
  #output_name = "test.txt"
  outputfile = output_name
  if (file.exists(outputfile)) {
  #Delete file if it exists
  file.remove(outputfile)
  }
  #sink(outputfile)
  idxs <- seq(1, length(dat1_total))
  #idxs = seq(1,100)
  fisher_pvs = c()
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
    M <- matrix(c(c.a1, s.a1, c.a2, s.a2), nrow = 2, dimnames = list(c("ref", "surv"), c("allele1", "allele2")))
    # Calculate delta_p
    p0 <- dat1$knownEM[i]
    p1 <- dat2$knownEM[i]
    delta_p <- p1 - p0
    # https://rdrr.io/cran/Exact/man/exact.test.html for barnard's exact test
    #p = fisher.test(M, alternative=alternative)$p.value
    t = exact.test(M, alternative="two.sided", method="z-pooled", )
    p=t$p.value
    cat(paste0(chr,'\t',pos,'\t',p1,'\t',p0,'\t',delta_p,'\t',p, '\n'), file=outputfile, append=TRUE)
  }
  #sink()
  #  return(M)
}


library("optparse")
 
option_list = list(
    make_option(c("-a", "--ref1"), type="character", default=NULL, 
              help="ref1", metavar="character"),
    make_option(c("-b", "--ch1"), type="character", default=NULL,
              help="ch1", metavar="character"),
    make_option(c("-o", "--out1"), type="character", default=NULL, 
              help="output1", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

fish.test(opt$ref1, opt$ch1, opt$out1)

