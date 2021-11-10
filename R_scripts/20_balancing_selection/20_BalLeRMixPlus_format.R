#################################################
# format the mafs file for BalLeRMixPlus        #
#################################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/20_Balancing_selection/format")
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/CHR19_REF19/")

format_mafs <- function(headname){
  dat <- read.table(paste0(headname,"_all_minq20_minmq30_CV30_masked.mafs"), header = T)
  #dat <- read.table(paste0(headname,"_minmapq30_minq20_CV30_masked_noinvers.mafs"), header = T)
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

format_dafs <- function(headname){
  dat <- read.table(paste0(headname,"_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs"), header = T)
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
      if (DT$major[i] == DT$anc[i]) {
        pos_ = c(pos_, pos)
        genPos_ = c(genPos_, genPos)
        g_a_ = c(g_a_, g_a1)
        g_tot_ = c(g_tot_, g_tot)
      } else if (DT$minor[i] == DT$anc[i]) {
        pos_ = c(pos_, pos)
        genPos_ = c(genPos_, genPos)
        g_a_ = c(g_a_, g_a2)
        g_tot_ = c(g_tot_, g_tot)
      } else {
        next # this is useful for B0 and B2 derived allele, which allows to remove the site when major or minor do not match the ancestral allele 
        
        # uncommand below when sites of substitution are needed for input
        # pos_ = c(pos_, pos) 
        # genPos_ = c(genPos_, genPos)
        # g_a1_ = c(g_a1_, sub)
        # g_tot_ = c(g_tot_, g_tot)
      }
    }
    output = data.frame(pos_, genPos_, g_a_, g_tot_)
    output = output[order(output$pos_), ]
    print(length(output$pos_))
    write.table(output, file = paste0(headname, "_", j , ".dafs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

format_mafs("NB")
format_mafs("COH")
format_mafs("REF20")
format_mafs("REF19")
format_mafs("SR")
format_mafs("CHR19")
format_mafs("CHR20")
format_mafs("HC")
format_mafs("ARN")
