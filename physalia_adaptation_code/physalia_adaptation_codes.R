# Libraries
library(dplyr)
library(magrittr)
library(tibble)
library(gplots)
library(RColorBrewer)
library(corrplot)
# --------------
makeSymm <- function(m, position) {
  # Add symetrical triangle matrix (upper or lower)
  if (position == 'upper'){
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  if (position == 'lower'){
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
}

#-------------- fst matrix for all SNPs --------------------
fst.mat <- read.table("03_fst/canada_fst_matrix.txt")
#use the function to fill the full matrix
fst.all.mat<- fst.mat %>%
  as.matrix(.) %>%
  makeSymm(., 'upper')
fst.all.mat[is.na(fst.all.mat)] <-  0 #replace NAs by 0 (NAs unaccepted for the heatmap function)
fst.all.mat[1:10,1:10] #check the fst_matrix

#visualise values
corrplot(fst.all.mat, is.corr = FALSE, method="number", addgrid.col = FALSE, diag=F, type="lower", number.digits = 3, number.cex=0.7)

#Visualize pairwise FST through a heatmap
gplots::heatmap.2(fst.all.mat, trace = 'none',
                  col= colorRampPalette(brewer.pal(9, "Reds"))(15),
                  key.xlab='FST')



library(SoDA)
library(reshape2)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)

#import information about populations
info_pop <- read.table("02_data/info_pop_geo_eco.txt", header=T)
head(info_pop)
#calculate geogrpahic (euclidian) distances between all pairs of populations
distance <- dist(SoDA::geoXY(info_pop$latitude, info_pop$longitude)) %>%
  as.matrix(.)
dimnames(distance) <- list(info_pop$pop,info_pop$pop)
distance

#prepare datasets
#linearize the distance matrix
dist.melt <- reshape2::melt(distance) %>%
  set_colnames(., c('pop1', 'pop2','distance'))
head(dist.melt)

#linearize the fst matrix
fst.melt <- reshape2::melt(fst.all.mat) %>%
  set_colnames(., c('pop1', 'pop2','FST'))

#join the distance and fst
IBD.df <- left_join(dist.melt, fst.melt, by=c('pop1','pop2')) %>%
  filter(., distance > 0)
head(IBD.df)

#test association with FST
cor.test(log(IBD.df$distance), IBD.df$FST/(1-IBD.df$FST))

#plot IBD
ggplot(IBD.df) + aes(x=log(distance), y=FST/(1-FST)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw()

info_ind<-read.table("02_data/info_samples_canada.txt", header=T)
head(info_ind)
table (info_ind$pop,info_ind$sex)



#1st intall the librairies

if (!("devtools" %in% installed.packages())){install.packages(devtools)}
library(devtools)
if (!("qvalue" %in% installed.packages())){BiocManager::install("qvalue")}
if (!("vcfR" %in% installed.packages())){install.packages("vcfR")} 
devtools::install_github("whitlock/OutFLANK")

#we use the library vcfR to convert the vcf into the OutFLANK format
library(OutFLANK)
library(vcfR)
library(ggplot2)
obj.vcfR <- read.vcfR("02_data/canada.vcf")

#extract useful informations about snp id and position
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP

#gather this info in a dataframe
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) # save info about id, chr, position
#R is sometimes not good at categorizing columns and here i had a problem that position was a factor... 
#this is an easy way to transform into a numeric
chr_pos$position<-as.numeric(as.character(chr_pos$position)) 

#we expect that it will be useful for subsequent analysis to have a file with snp id and position so let's write it in our folder 02_data
write.table(chr_pos, "02_data/SNP_pos.txt", sep="\t", quote=F, row.names=F)

#extract and format genotype matrix
geno <- extract.gt(obj.vcfR) # Character matrix containing the genotypes

#an empty matrix, (9 stands for missing data)
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno))

#that we fill with genotypes
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

#an overview of our data and its first 10 rows/10 columns
table(as.vector(G))
dim(G)
G[1:10,1:10]

#As it will be useful later, I suggest that we export it now
write.table(G, "02_data/geno_matrix.txt", sep="\t", col.names=F, row.names=F)

#import pop info
info_samples_canada<- read.table("02_data/info_samples_canada.txt", header=T)
head(info_samples_canada)
pop_vector<- info_samples_canada$pop

# FST matrix with OutFLANK
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector)


#import pruned info
id_snp_pruned<-read.table("02_data/canada.prune.in")

#those are SNPs id, we need to know at which position they are
#this can be done with the %in% function
lines_trim<-which(id_snp %in% id_snp_pruned[,1]) 
head(lines_trim)
length(lines_trim)

#run outFLANK on pruned SNPs
#numberOfSamples is the number of populations
#qthreshold is the false discovery rate
out_trim <- OutFLANK(my_fst[which(id_snp %in% id_snp_pruned[,1]),], NumberOfSamples=20, qthreshold = 0.05, Hmin = 0.1)
str(out_trim)

#have a look at the results
#the jpeg line allow to output an image in your folder that you can later download to have a look at
jpeg("04_outflank/outflank_prunedSNP_fst.jpeg")
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL)
dev.off()

jpeg("04_outflank/outflank_prunedSNP_pvalues.jpeg")
hist(out_trim$results$pvaluesRightTail)
dev.off()


P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
head(P1)

#we need to add the chromosome/position info for plotting. the left_join function in dplyr is super useful to match differebnt table
library(dplyr)
P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))

#We can have a look at the results by exporting the figures
#we can look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
#And we can look along the genome with our manhattan plots

jpeg("04_outflank/outflank_outlier_fst_He.jpeg")
ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag))+ 
  geom_point()+
  theme_classic()
dev.off()

#note that we divide here position by 1000000 so the scale is in MB
jpeg("04_outflank/outflank_outlier_fst.jpeg")
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ 
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x")+
  labs(  x = "position (in MB)")
dev.off()

#It may also be easier to export the matrix and play in Rstudio
write.table(P1_pos, "04_outflank/outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)










args = commandArgs(trailingOnly=TRUE)

#source functions from Gautier file
source("baypass_utils.R")

#load omega matrix
omega = as.matrix(read.table("05_baypass/prunedsnps.output_mat_omega.out"))

#load beta parameters
pi.beta.coef=read.table("05_baypass/allsnps.output_summary_beta_params.out", h=T)$Mean

#create the POD - simulate
simu.bta <- simulate.baypass(omega.mat=omega, nsnp=500,beta.pi=pi.beta.coef, suffix="simulates_pods")





# run model for 5 repetitions
# for i in 1 2 3 4 5
# do
# seed=$((1000 + RANDOM % 9999)) #set random seed
# ./g_baypass -npop $N_pop -gfile $geno -efile $env -omegafile $omega_mat -outprefix output"$i" -nthreads $N_CPU -npilot 25 -burnin 5000 -seed $seed
# done



#load xtx values
xtx_allsnps<-read.table("05_baypass/allsnps.output_summary_pi_xtx.out", header = T)
head(xtx_allsnps)
#we will mostly work with M_WtX

#load position info about the SNPs.
SNP_pos<-read.table("02_data/SNP_pos.txt", header=T)
#should be same nb of rows
dim(xtx_allsnps)
dim(SNP_pos)

xtx_pos<-cbind(SNP_pos, xtx_allsnps)

ggplot(xtx_pos, aes(x=position, y=M_XtX, colour=chromosome))+ 
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x")




#load xtx values from simulatd data
xtx_simu<-read.table("05_baypass/simulate.output_summary_pi_xtx.out", header=T)
head(xtx_simu)

#calculate the threshold
threshold_fdr0.01 = quantile(xtx_simu$M_XtX,probs=0.99)
threshold_fdr0.05 = quantile(xtx_simu$M_XtX,probs=0.95)

#add it on the plot
ggplot(xtx_pos, aes(x=position, y=M_XtX, colour=chromosome))+ 
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x")+
  geom_hline(aes(yintercept =threshold_fdr0.05), linetype="dotted", size=1, col="red", show.legend = FALSE)+
  geom_hline(aes(yintercept =threshold_fdr0.01), linetype="dotted", size=1, show.legend = FALSE)

#output outliers
xtx_pos[xtx_pos$M_XtX>= threshold_fdr0.05,]








#read file
info_pop<-read.table("02_data/info_pop_geo_eco.txt", header =T)
head(info_pop)
#env vector and its transpose form
info_pop$temperature
t(info_pop$temperature)
write.table (t(info_pop$temperature), "05_baypass/env.txt", sep="\t", quote=F, row.names=F, col.names=F)







#load bf values
bf_allsnps<-read.table("05_baypass/allsnps_env.output_summary_betai_reg.out", header = T)
bf_pos<-cbind(SNP_pos, bf_allsnps)

#load bf values from simulatd data
bf_simu<-read.table("05_baypass/simulate_env.output_summary_betai_reg.out", header=T)

#calculate the threshold from simu (or you cna use BF = 10)
threshold_fdr0.01 = quantile(bf_simu$BF.dB,probs=0.99)
threshold_fdr0.05 = quantile(bf_simu$BF.dB,probs=0.95)

outliers<-bf_pos[bf_pos$BF.dB>=threshold_fdr0.05,]
write.table(outliers, "05_baypass/outlier_temp_bp.txt", row.names=F, quote=F, sep="\t")





geno<-read.table("02_data/geno_matrix.txt")
SNP_pos<-read.table("02_data/SNP_pos.txt", header=T)

#transpose data and give meaningful colnames
gen<-t(geno)
colnames(gen)<-paste(SNP_pos$chromosome, SNP_pos$position, sep="_")
gen [1:10,1:10]

#replace 9 by NA
gen[which(gen=="9")]<- NA
#evaluate % of missing
sum(is.na(gen))/(dim(gen)[1]*dim(gen)[2]) # about 1.5% of missing data
#impute missing with the most common geno
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs





info<-read.table("02_data/info_samples_canada.txt", header=T)
head(info)




library(vegan)
temp.rda <- rda(gen.imp ~ info$temperature, scale=T)
temp.rda


RsquareAdj(temp.rda)

temp.signif.full <- anova.cca(temp.rda, parallel=getOption("mc.cores")) # default is permutation=999
temp.signif.full




jpeg("06_rda/rda1_triplot.jpeg")
plot(temp.rda, scaling=3) 
points(temp.rda, display="sites", pch=20, cex=1.3, col=info$pop, scaling=3)
dev.off()







load.temp.rda <- scores(temp.rda, choices=c(1), display="species") 

#load info about snps
SNP_pos<-read.table("02_data/SNP_pos.txt", header=T)
load.temp.rda.pos<-cbind(SNP_pos,load.temp.rda )
head(load.temp.rda.pos)

#plot distribution
jpeg("06_rda/rda1_loading_hist.jpeg")
hist(load.temp.rda.pos$RDA1, main="Loadings on RDA1")
dev.off()

#chose the sd limit
z=2
lim_min<- mean(load.temp.rda.pos$RDA1) - ( z * sd(load.temp.rda.pos$RDA1) )
lim_max<- mean(load.temp.rda.pos$RDA1) + ( z * sd(load.temp.rda.pos$RDA1) )

#outliers
outlier_temp <- load.temp.rda.pos[load.temp.rda.pos$RDA1 >=lim_max | load.temp.rda.pos$RDA1 <= lim_min ,]
outlier_temp

#export them
write.table(outlier_temp, "06_rda/outlier_temp_rda.txt", row.names=F, quote=F, sep="\t")





library(ggplot2)

#plot
jpeg("06_rda/loading_temp_manhattanplot.jpeg")
ggplot(load.temp.rda.pos, aes(x=position, y=RDA1, colour=chromosome))+ 
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x")+
  geom_hline(aes(yintercept =lim_min), linetype="dotted", size=1,  show.legend = FALSE)+
  geom_hline(aes(yintercept =lim_max), linetype="dotted", size=1, show.legend = FALSE)
dev.off()





geo.rda <- rda(gen.imp ~ info$lat + info$long, scale=T)
geo.rda
RsquareAdj(geo.rda)

library('network')
jpeg("06_rda/rda1_triplot_geo.jpeg")
plot(geo.rda, scaling=3) 
points(geo.rda, display="sites", pch=20, cex=1.3, col=as.color(info$pop), scaling=3)
dev.off()



temp.geo.rda <- rda(X= gen.imp, Y= info$temp, Z=cbind(info$lat,info$long), scale=T)
temp.geo.rda
RsquareAdj(temp.rda)
RsquareAdj(temp.geo.rda)





library(cluster)
#Create a geographic matrix with dbmem
coorgeo<-info[,6:7]
dist_eucl<-daisy(coorgeo)
geo_pcnm<-pcnm(dist_eucl)
geo_mat<-geo_pcnm$vectors
head(geo_mat)

#perform rda
geo.rda<-rda(gen.imp~geo_mat[,1]+geo_mat[,2]+geo_mat[,3]+geo_mat[,4]+geo_mat[,5]+geo_mat[,6]+geo_mat[,7]+geo_mat[,8])
RsquareAdj(geo.rda)
#select model
anova(geo.rda, step=1000, by= "margin")
ordistep(geo.rda)


#----------------------

sex.rda <- rda(gen.imp ~ info$sex, scale=T)

RsquareAdj(sex.rda)

sex.signif.full <- anova.cca(sex.rda, parallel=getOption("mc.cores")) # default is permutation=999
sex.signif.full

load.sex.rda <- scores(sex.rda, choices=c(1), display="species") 
load.sex.rda.pos<-cbind(load.sex.rda, SNP_pos)

#plot
jpeg("06_rda/loading_sex_manhattanplot.jpeg")
ggplot(load.sex.rda.pos, aes(x=position, y=RDA1, colour=chromosome))+ 
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x")
dev.off()

#---------------------

temp_gps.rda <- rda(gen.imp ~ info$temperature+ info$lat + info$long, scale=T)
temp_gps.rda

RsquareAdj(temp_gps.rda)

signif.axis <- anova.cca(temp_gps.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
vif.cca(temp_gps.rda)

#--------------------


library(dplyr)

#load outliers tables
outlier_temp_rda<-read.table("06_rda/outlier_temp_rda.txt", header=T)
head(outlier_temp_rda)
nRDA<-dim(outlier_temp_rda)[1]
nRDA #how many outliers?

outlier_temp_bp<-read.table("05_baypass/outlier_temp_bp.txt", header =T)
head(outlier_temp_bp)
outlier_temp_bp<-outlier_temp_bp[,c(1,2,3,8)] #we keep snp id, chr, pos and BF
dim(outlier_temp_bp)
nBP<-dim(outlier_temp_bp)[1]
nBP #how many outliers?

#join outliers keeping positions present in either the 1st or the 2nd database (or both)
outlier_temp_fulljoin<-full_join(outlier_temp_rda,outlier_temp_bp)
head(outlier_temp_fulljoin)
nALL<-dim(outlier_temp_fulljoin)[1]
nALL # how many in total?

#join outliers keeping positions present in either the 1st or the 2nd database (or both)
outlier_temp_innerjoin<-inner_join(outlier_temp_rda,outlier_temp_bp)
head(outlier_temp_innerjoin)
dim(outlier_temp_innerjoin)
nboth<-dim(outlier_temp_innerjoin)[1]
nboth #how many joint outliers?

#visualize
library(ggVennDiagram)
ggVennDiagram(list(rda = 1:nRDA, BP = (nRDA+1-nboth):(nRDA-nboth+nBP)))


