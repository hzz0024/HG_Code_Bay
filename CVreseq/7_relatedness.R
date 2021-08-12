### Calculate Relatedness ####
setwd("~/Documents/Ryan_workplace/CVreseq_pcadapt")
vcf <- read.vcfR("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.thin.recode.vcf.gz")
geno <- vcf@gt[,-1] 
genosub <- geno
G_relate <- matrix(NA, nrow=ncol(genosub), ncol=nrow(genosub)*2)
abob <- t(substr(genosub,start = 1, stop=1))
class(abob) <- "numeric"
abob = abob + 1
bbob <- t(substr(genosub,start = 3, stop=3))
class(bbob) <- "numeric"
bbob = bbob + 1
odd <- seq(1, ncol(G_relate), by=2)
G_relate[,odd] <- abob
G_relate[,odd+1] <- bbob
rownames(G_relate) <- rownames(abob)
position = getPOS(vcf)
colnames(G_relate) <- rep(position, each=2)
class(G_relate) <- "numeric"
G_relate0 = G_relate

ind_group = rep(c(1,2,3,4,5),each=6)
G_relate <- data.frame(Individual = rownames(G_relate), population=ind_group, G_relate)
head(G_relate[,1:11])
G_relate2 <- G_relate[,-2]

head(G_relate2[,1:11])
G_relate2$Individual <- as.character(G_relate$Individual)
#G_relate2[,2:ncol(G_relate2)] <- G_relate2[,2:ncol(G_relate2)]+1 # 0s read as missing data
head(G_relate2[,1:11])
library('related')
poprelate2 <- coancestry(G_relate2,  lynchrd =1 )
(rem_ind <- unique(poprelate2$relatedness$ind2.id[which(poprelate2$relatedness$lynchrd>0.5)]))
ind_keep <- which(!(G_relate2$Individual %in% rem_ind))
print(c("final sample size:", length(ind_keep)))
# relatedness plotting
relateness = poprelate2$relatedness
write.table(relateness, file = "relatedness.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
filename = '/out.relatedness2'
dat = read.csv(filename, sep='\t')
relateness  = matrix(relateness, nrow=30)
rownames(relateness) =  dat[,2][1:30]
colnames(relateness) =  dat[,2][1:30]
#heatmap(relateness, Colv = NA, Rowv = NA)
library('glots')
heatmap.2(relateness, trace="none", density.info='none')