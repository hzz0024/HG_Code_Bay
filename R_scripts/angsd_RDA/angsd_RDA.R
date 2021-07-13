# Detecting multilocus adaptation using Redundancy Analysis (RDA)
#install.packages(c("psych","vegan"), dependencies=TRUE)

# Load packages
# -------------
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA

# below is the script scriptRDA.r in https://github.com/clairemerot/selection_analysis/blob/master/01_scripts/Rscripts/scriptRDA.r
# code R
library(vegan)

argv <- commandArgs(T)
MAF_FILE <- argv[1]
OUTPUT_FOLDER <- argv[2]
ENV<-argv[3]
ENV_FOLDER<-argv[4]



data_maf<-read.table(MAF_FILE, header=TRUE)
head(data_maf)


data_freq<-t(data_maf)
data_env<-read.table(paste0(ENV_FOLDER, "/",ENV), header=TRUE)

###rda
#perform rda
pop6rda<-rda(data_freq~., data=data_env, scale=T)
pop6rda

#look at R² and important axis
RsquareAdj(pop6rda)
summary(pop6rda)$concont

write.table(RsquareAdj(pop6rda), paste0(OUTPUT_FOLDER,ENV,"_adjR2.txt"), row.names=FALSE, quote=FALSE)
pdf(file=paste0(OUTPUT_FOLDER,ENV,"_screeplot.pdf"))
screeplot(pop6rda)
dev.off()


#plot RDA SNP are red at the cener of the plot
#population are black #vector are env predictorss
jpeg(file=paste0(OUTPUT_FOLDER,ENV,"_rda_1-2_1-3.jpg"))
par(mfrow=c(1,2))
plot(pop6rda, scaling=3)
text(pop6rda, display="sites", col=1, scaling=3)
plot(pop6rda, choices=c(1,3), scaling=3)#axis 1 and 3
text(pop6rda, display="sites", col=1, scaling=3)
dev.off()

#extract significative snp
load.rda<-summary(pop6rda)$species[,1:3]
pdf(file=paste0(OUTPUT_FOLDER,ENV,"_hist_loadings.pdf"))
par(mfrow=c(1,3))
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")
dev.off()

write.table (load.rda[,1:3], paste0(OUTPUT_FOLDER,ENV,"_rda_loading.txt"), quote=FALSE)

#function
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#applyfuntcion to outliers with rda 1 and rda 2
#z is the number of sd 3 is standard, 2 less conservative, 
cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 
#cand3 <- outliers(load.rda[,3],3) 
ncand<-length(cand1) + length ( cand2) #+ length ( cand3)#total nb of snps
lim<-3

if (ncand==0){
  cand1 <- outliers(load.rda[,1],2.5) 
  cand2 <- outliers(load.rda[,2],2.5) 
  #cand3 <- outliers(load.rda[,3],2) 
  ncand<-length(cand1) + length ( cand2) #+ length ( cand3)#total nb of snps
  lim<-2.5
  print("no SNP outlier with 3 sd, re-run with 2.5sd")
}


#Next, we'll organize our results by making one data frame with the axis, SNP name, loading, & correlation with each predictor:
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
#cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1)=colnames(cand2)=c("axis", "snp", "loading")#=colnames(cand3)
cand <- rbind(cand1, cand2)#, cand3)
cand$snp <- as.character(cand$snp)
head(cand)

#Let's add in the correlations of each candidate SNP with the  environmental predictors:
n<-dim(data_env)[2]
foo <- matrix(nrow=ncand, ncol=n)  # n columns for n predictors
colnames(foo) <- colnames(data_env)

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- data_freq[,nam]
  foo[i,] <- apply(data_env,2,function(x) cor(x,snp.gen))
}
#table of candidate snp with loading on each axis and correlation with env predictors
cand <- cbind.data.frame(cand,foo)  
head(cand)

#investigate the candidates
#are there snp associated with several axis? if yes remove them
n_dupli<-length(cand$snp[duplicated(cand$snp)])
n_dupli
if (n_dupli>=1){cand<-cand[!duplicated(cand$snp)]}

#Next, we'll see which of the predictors each candidate SNP is most strongly correlated with:
n<-dim(cand)[2]
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,n+1] <- names(which.max(abs(bar[4:n]))) # gives the variable
  cand[i,n+2] <- max(abs(bar[4:n]))              # gives the correlation
}

colnames(cand)[n+1] <- "predictor"
colnames(cand)[n+2] <- "correlation"
head(cand)
table(cand$predictor) 
write.table(cand, paste0(OUTPUT_FOLDER,ENV,"_candidate_SNP_",lim,"_sd.txt"), quote=FALSE, sep=" ", row.names=FALSE)

#Let's look at RDA plots again, but this time focus in on the SNPs in the ordination space. 
#We'll color code the SNPs based on the predictor variable that they are most strongly correlated with. 

sel <- cand$snp
env <- cand$predictor
variables<-levels(as.factor(cand$predictor))
for (i in 1: length(variables))
{env[env==variables[i]] <- i} #associates a color with each env predictor

# color by predictor:
col.pred <- rownames(pop6rda$CCA$v) # pull the SNP names and put them in a vector of colour

for (i in 1:length(sel)) {           # 
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("LG",col.pred)] <- '#f1eef6' # non-candidate SNPs make them transparent
col.pred[grep("scaffold",col.pred)] <- '#f1eef6' # non-candidate SNPs make them transparent
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c(1: length(variables))

#Now we're ready to plot the SNPs:
# axes 1 & 2
jpeg(paste0(OUTPUT_FOLDER,ENV,"_rda_SNP.jpg"))
#par(mfrow=c(1,2))
plot(pop6rda, type="n", scaling=1, xlim=c(-1,1), ylim=c(-1,1))
points(pop6rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=1)
points(pop6rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=1)
text(pop6rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=variables, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
#plot(pop6rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
#points(pop6rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
#points(pop6rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
#text(pop6rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
#legend("bottomright", legend=variables, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()

#significance of the model with permutation
signif.full<-anova.cca(pop6rda, parallel=4)
signif.full
write.table(signif.full, paste0(OUTPUT_FOLDER,ENV,"_signifmodel.txt"), quote=FALSE, row.names=FALSE)

#significance of each axis (take time)
signif.axis<-anova.cca(pop6rda, by="axis",parallel=4)
signif.axis
write.table(signif.axis, paste0(OUTPUT_FOLDER,ENV,"_signifaxisrda.txt"), quote=FALSE, row.names=FALSE)

#variance inflexion factor
vif.cca(pop6rda)

{"mode":"full","isActive":false}