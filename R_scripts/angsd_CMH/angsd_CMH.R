
# The null hypothesis is that the relative proportions of one variable are independent of the other variable within the repeats; in other words, there is no consistent difference in proportions in the 2×2 tables.

# For our imaginary legwarmers experiment, the null hypothesis would be that the proportion of people feeling pain was the same for legwarmer-wearers and non-legwarmer wearers, 
# after controlling for the time of year. 

# The alternative hypothesis is that the proportion of people feeling pain was different for legwarmer and non-legwarmer wearers.

# If some repeats have a big difference in proportion in one direction, and other repeats have a big difference in proportions but in the opposite direction, 
# the Cochran-Mantel-Haenszel test may give a non-significant result.

########################
#######  CMH test ######
########################

library(dplyr)
library(vcd)

# data entered as a data frame
CMH.test <- function(ref_name1, ch_name1, ref_name2, ch_name2, ref_name3, ch_name3, output_name){
  # ref_name1 = "REF19_minq20_minmq30_1x_CV30_masked.mafs"
  # ch_name1 = "CHR19_minq20_minmq30_1x_CV30_masked.mafs" 
  # ref_name2 = "REF20_minq20_minmq30_1x_CV30_masked.mafs"
  # ch_name2 = "CHR20_minq20_minmq30_1x_CV30_masked.mafs"
  # ref_name3 = "NB_minq20_minmq30_1x_CV30_masked.mafs"
  # ch_name3 = "HC_minq20_minmq30_1x_CV30_masked.mafs"
  dat1 <- read.delim(ref_name1, header = TRUE, sep='\t')
  dat1_total = dat1$nInd*2
  dat1_refCount = round(dat1$knownEM*dat1_total)
  dat1_altCount = dat1_total - round(dat1$knownEM*dat1_total)
  # load challenge file with header in it
  dat2 <- read.delim(ch_name1, header = TRUE, sep='\t')
  dat2_total = dat2$nInd*2
  dat2_refCount = round(dat2$knownEM*dat2_total)
  dat2_altCount = dat2_total - round(dat2$knownEM*dat2_total)
  
  dat3 <- read.delim(ref_name2, header = TRUE, sep='\t')
  dat3_total = dat3$nInd*2
  dat3_refCount = round(dat3$knownEM*dat3_total)
  dat3_altCount = dat3_total - round(dat3$knownEM*dat3_total)
  # load challenge file with header in it
  dat4 <- read.delim(ch_name2, header = TRUE, sep='\t')
  dat4_total = dat4$nInd*2
  dat4_refCount = round(dat4$knownEM*dat4_total)
  dat4_altCount = dat4_total - round(dat4$knownEM*dat4_total)
  
  dat5 <- read.delim(ref_name3, header = TRUE, sep='\t')
  dat5_total = dat5$nInd*2
  dat5_refCount = round(dat5$knownEM*dat5_total)
  dat5_altCount = dat5_total - round(dat5$knownEM*dat5_total)
  # load challenge file with header in it
  dat6 <- read.delim(ch_name3, header = TRUE, sep='\t')
  dat6_total = dat6$nInd*2
  dat6_refCount = round(dat6$knownEM*dat6_total)
  dat6_altCount = dat6_total - round(dat6$knownEM*dat6_total)
  
  outputfile = output_name
  if (file.exists(outputfile)) {
    #Delete file if it exists
    file.remove(outputfile)
  }
  idxs <- seq(1, length(dat1_total))
  #idxs = seq(1,10)
  fisher_pvs = c()
  for(i in idxs){
    # print the running process
    s = paste0(i,'/',dim(dat1)[1])
    message(s,"\r",appendLF=FALSE)
    # chromosome and position information
    chr <- dat1$chromo[i]
    pos <- dat1$position[i]
    
    cnt_column = c(dat1_refCount[i], dat2_refCount[i], dat1_altCount[i], dat2_altCount[i], 
                   dat3_refCount[i], dat4_refCount[i], dat3_altCount[i], dat4_altCount[i],
                   dat5_refCount[i], dat6_refCount[i], dat5_altCount[i], dat6_altCount[i])
    Input =("
   contrast  pop         Allele   
       2019  reference        1     
       2019  chalelnge        1     
       2019  reference        2     
       2019  chalelnge        2     
       2020  reference        1     
       2020  chalelnge        1     
       2020  reference        2     
       2020  chalelnge        2     
       wild  reference        1     
       wild  chalelnge        1     
       wild  reference        2     
       wild  chalelnge        2     
    ")
    Data = read.table(textConnection(Input),header=TRUE)
    Data$Count = cnt_column
    ### Specify the order of factor levels
    ### Otherwise, R will alphabetize them
    Data =
      mutate(Data,
             contrast = factor(contrast, levels=unique(contrast)),
             pop = factor(pop, levels=unique(pop)),
             Allele = factor(Allele, levels=unique(Allele))
      )
    ### Cross-tabulate the data
    ###   Note here, Location is stratum variable (is last)
    ###              Habitat x Allele are 2 x 2 tables
    Data.xtabs = xtabs(Count ~ Allele + pop + contrast,
                       data=Data)
    ftable(Data.xtabs)
    # Cochran–Mantel–Haenszel test
    p = mantelhaen.test(Data.xtabs)$p.value
    # Woolf test
    oddsratio(Data.xtabs, log=TRUE) 
    wp = woolf_test(Data.xtabs)$p.value
    # Individual Fisher exact tests
    # n = dim(Data.xtabs)[3]
    # 
    # for(i in 1:n){
    #   Name = dimnames(Data.xtabs)[3]$Location[i]
    #   P.value = fisher.test(Data.xtabs[,,i])$p.value
    #   cat(Name, "\n")
    #   cat("Fisher test p-value: ", P.value, "\n")
    #   cat("\n")
    # }
    p1 <- dat1$knownEM[i]
    p2 <- dat2$knownEM[i]
    delta_p1 <- p2 - p1
    p3 <- dat3$knownEM[i]
    p4 <- dat4$knownEM[i]
    delta_p2 <- p4 - p3
    p5 <- dat5$knownEM[i]
    p6 <- dat6$knownEM[i]
    delta_p3 <- p6 - p5
    cat(paste0(chr,'\t',pos,'\t',delta_p1, '\t',delta_p2,'\t',delta_p3,'\t',p,'\t', wp,'\n'), file=output_name, append=TRUE)
  }
  #  return(M)
}

CMH.test("REF19_minq20_minmq30_1x_CV30_masked.mafs", "CHR19_minq20_minmq30_1x_CV30_masked.mafs", 
         "REF20_minq20_minmq30_1x_CV30_masked.mafs", "CHR20_minq20_minmq30_1x_CV30_masked.mafs", 
         "NB_minq20_minmq30_1x_CV30_masked.mafs", "HC_minq20_minmq30_1x_CV30_masked.mafs", 
         "CMH_test.txt")

library("optparse")

option_list = list(
  make_option(c("-a", "--ref_name1"), type="character", default=NULL, 
              help="ref_name1", metavar="character"),
  make_option(c("-b", "--ch_name1"), type="character", default=NULL,
              help="ch_name1", metavar="character"),
  make_option(c("-c", "--ref_name2"), type="character", default=NULL,
              help="ref_name2", metavar="character"),
  make_option(c("-d", "--ch_name2"), type="character", default=NULL,
              help="ch_name2", metavar="character"),
  make_option(c("-c", "--ref_name3"), type="character", default=NULL,
              help="ref_name3", metavar="character"),
  make_option(c("-d", "--ch_name3"), type="character", default=NULL,
              help="ch_name3", metavar="character"),
  make_option(c("-l", "--output_name"), type="character", default=NULL,
              help="output_name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

CMH.test(opt$ref_name1, opt$ch_name1, opt$ref_name2, opt$ch_name2, opt$ref_name3, opt$ch_name3, opt$output_name)

