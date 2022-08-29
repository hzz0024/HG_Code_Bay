################################
#  process outlier for B-score #
################################

######################################
# Step 1: generate the spectrum file # 
######################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/20_balancing_selection/sfs/")

#function to normalize
norm <- function(x) x/sum(x)

spec_format <- function(sfs_name){
  sfs<-(scan(sfs_name))
  sfs_<-norm(sfs[-c(1,length(sfs))])
  barplot(sfs_) #plot variable sites 
  df <- data.frame(seq(1:length(sfs_)), length(sfs_)+1, sfs_)
  write.table(df, file = paste0("./spect_results/",sfs_name, ".nosub_spect.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

spec_format("NB_18.neutral.sfs")
spec_format("NB_19.neutral.sfs")
spec_format("NB_21.neutral.sfs")
spec_format("Ref_19.sfs")
spec_format("Ref_20.sfs")
spec_format("HC_18.neutral.sfs")
spec_format("HC_19.neutral.sfs")
spec_format("HC_21.neutral.sfs")
spec_format("Sur_19.sfs")
spec_format("Sur_20.sfs")

#######################################
# Step 2: convert mafs to B0maf format 
#######################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/20_balancing_selection/outliers")
format_mafs <- function(headname, pop_n){
  dat <- read.table(paste0(headname,"_shared_outliers.mafs"), header = T)
  # do formatting for each chromosome
  for (j in c( 'NC_035780.1', 'NC_035781.1', 'NC_035782.1', 'NC_035783.1', 'NC_035784.1', 'NC_035785.1', 'NC_035786.1', 'NC_035787.1', 'NC_035788.1', 'NC_035789.1')){ 
    DT = dat[which(dat$chromo %in% j), ]
    # filter out SNPs with maf < 0.05 or maf > 0.95 in each population
    chr_ = c()
    pos_ = c()
    genPos_ = c()
    g_a_ = c()
    g_tot_ = c()
    #for(i in seq(10)){
    for (i in seq(length(DT[, 1]))) {
      chr = DT$chromo[i]
      pos = DT$position[i]
      genPos = "NA"
      sub = "0"
      g_tot = pop_n * 2
      g_a1 = round(DT$knownEM[i] * pop_n * 2)
      g_a2 = g_tot - g_a1
      if (DT$knownEM[i] < 0.5) {
        pos_ = c(pos_, pos)
        genPos_ = c(genPos_, genPos)
        g_a_ = c(g_a_, g_a1)
        g_tot_ = c(g_tot_, g_tot)
      } else {
        pos_ = c(pos_, pos)
        genPos_ = c(genPos_, genPos)
        g_a_ = c(g_a_, g_a2)
        g_tot_ = c(g_tot_, g_tot)
      }
    }
    output = data.frame(pos_, genPos_, g_a_, g_tot_)
    output = output[order(output$pos_), ]
    print(length(output$pos_))
    write.table(output, file = paste0(headname, "_", j , ".mafs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

format_mafs("HC_18", 48)
format_mafs("HC_19", 48)
format_mafs("HC_21", 50)
format_mafs("Sur_19", 50)
format_mafs("Sur_20", 48)

format_mafs("NB_18", 49)
format_mafs("NB_19", 48)
format_mafs("NB_21", 50)
format_mafs("Ref_19", 47)
format_mafs("Ref_20", 247)




########### old trash format function ############

format_mafs <- function(headname){
  dat <- read.table(paste0(headname,"_shared_outliers.mafs"), header = T)
  # do formatting for each chromosome
  for (j in c( 'NC_035780.1', 'NC_035781.1', 'NC_035782.1', 'NC_035783.1', 'NC_035784.1', 'NC_035785.1', 'NC_035786.1', 'NC_035787.1', 'NC_035788.1', 'NC_035789.1')){ 
    DT = dat[which(dat$chromo %in% j), ]
    # filter out SNPs with maf < 0.05 or maf > 0.95 in each population
    chr_ = c()
    pos_ = c()
    genPos_ = c()
    g_a_ = c()
    g_tot_ = c()
    #for(i in seq(10)){
    for (i in seq(length(DT[, 1]))) {
      chr = DT$chromo[i]
      pos = DT$position[i]
      genPos = "NA"
      sub = "0"
      g_tot = DT$nInd[i] * 2
      g_a1 = round(DT$knownEM[i] * DT$nInd[i] * 2)
      g_a2 = g_tot - g_a1
      if (DT$knownEM[i] < 0.5) {
        pos_ = c(pos_, pos)
        genPos_ = c(genPos_, genPos)
        g_a_ = c(g_a_, g_a1)
        g_tot_ = c(g_tot_, g_tot)
      } else {
        pos_ = c(pos_, pos)
        genPos_ = c(genPos_, genPos)
        g_a_ = c(g_a_, g_a2)
        g_tot_ = c(g_tot_, g_tot)
      }
    }
    output = data.frame(pos_, genPos_, g_a_, g_tot_)
    output = output[order(output$pos_), ]
    print(length(output$pos_))
    write.table(output, file = paste0(headname, "_", j , ".mafs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}
