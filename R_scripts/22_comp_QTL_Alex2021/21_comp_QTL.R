setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/21_comp_QTL/")

pname = "ps_Del19_HC_NB.txt"
dat = read.delim(pname, header = F, sep='\t')
dat$id = paste0(dat[,1],'_',dat[,2])
head(dat)

pname = "HC_SR_FDR_outlier.list"
dat = read.delim(pname, header = T, sep='\t')

qname = "QLT_outlier.txt"
QTL = read.delim(qname, header = TRUE, sep='\t')
QTL$id = paste0(QTL$chr,'_',QTL$pos)
head(QTL)

intersect(QTL$id, dat$id)

# no shared QTL
# only 6 SNP were included in my whole snp list
# Trait	Chromosome	Position (bp)	-log10(p)	R2	Gene
# 1	26229830	7.39	0.097	rho GTPase-activating protein 190-like
# 1	27205754	6.80	0.089	uncharacterized
# 1	25822872	6.65	0.157	nuclear receptor coactivator 2-like
# 1	26221696	6.16	0.128	rho GTPase-activating protein 190-like	
# 1	27372715	6.14	0.124	transient receptor potential cation channel subfamily M member 1-like	
# 1	26221697	5.79	0.119	rho GTPase-activating protein 190-like


format_bed <- function(pname, distance){
  #pname = "challenge_FDR_outlier.list"
  dat = read.delim(pname, header = T, sep='\t')
  message(paste0("total number of outliers is ", length(dat$chr)))
  dat_ = dat[with(dat, order(chr, pos)),]
  bed_list <- paste0(dat_$chr, "\t" , dat_$pos-distance, "\t", dat_$pos+distance)
  print(head(bed_list))
  write.table(bed_list, paste0(pname, ".bed"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_bed("challenge_FDR_outlier.list", 0)
format_bed("HC_NB_FDR_outlier.list", 0)

