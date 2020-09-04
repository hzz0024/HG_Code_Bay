########### Combine two p-values ##############
combine_p <- function(p1_name, p2_name){
  dat_1 <- read.delim(p1_name, header = FALSE, sep='\t')
  dat_2 <- read.delim(p2_name, header = FALSE, sep='\t')
  p1 <- dat_1$V7
  p2 <- dat_2$V7
  ps <- data.frame(p1, p2)
  return(ps)
}

HC_COH.ps <- combine_p("fish_CH_REF_HCCOH_anc.txt", "fish_HC_COH_anc.txt")
HC_COH.ps <- combine_p("fish_CH_REF_HCCOH.txt", "fish_HC_COH.txt")


p = get.combined(HC_COH.ps[1:5000,])
qobj <- qvalue(p = p)
fdr <- qobj$qvalues
length(fdr[fdr<0.05])

p = HC_COH.ps
p = get.combined(p) 
adj = p.adjust(p, method = 'BH')
length(p[adj<0.05])


