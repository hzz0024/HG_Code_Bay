###############################################
#######  process random snps for B-score ######
###############################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/20_balancing_selection/random")

#######################################
# Step 1: convert mafs to B0maf format 
#######################################

format_mafs <- function(headname, k, pop_n){
  dat <- read.table(paste0(headname,"_random_", k, ".mafs"), header = T)
  # do formatting for each chromosome
  for (j in c( 'NC_035780.1', 'NC_035781.1', 'NC_035782.1', 'NC_035783.1', 'NC_035784.1', 'NC_035785.1', 'NC_035786.1', 'NC_035787.1', 'NC_035788.1', 'NC_035789.1')){ 
    DT = dat[which(dat$chromo %in% j), ]
    if(dim(DT)[1]==0) next
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
    write.table(output, file = paste0(headname,"_",k, "_",j, ".mafs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}
  


setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/20_balancing_selection/random")

for (k in seq(1:100)){
  format_mafs("HC_18", k, 48)
} 

for (k in seq(1:100)){
  format_mafs("HC_19", k, 48)
} 

for (k in seq(1:100)){
  format_mafs("HC_21", k, 50)
} 

for (k in seq(1:100)){
  format_mafs("Sur_19", k, 50)
} 

for (k in seq(1:100)){
  format_mafs("Sur_20", k, 48)
} 

for (k in seq(1:100)){
  format_mafs("NB_18", k, 49)
}

for (k in seq(1:100)){
  format_mafs("NB_19", k, 48)
}

for (k in seq(1:100)){
  format_mafs("NB_21", k, 50)
}

for (k in seq(1:100)){
  format_mafs("Ref_19", k, 47)
}

for (k in seq(1:100)){
  format_mafs("Ref_20", k, 247)
}






