# install.packages(c("psych","vegan"), dependencies=TRUE)
# website https://popgen.nescent.org/2018-03-27_RDA_GEA.html
# Load packages
# -------------
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA

setwd("~/Downloads")
datzip <- ("data/wolf_geno_samp_10000.zip") 
zipd <- tempdir()
unzip(datzip, exdir=zipd)
gen <- read.csv(paste0(zipd,"/wolf_geno_samp_10000.csv"), row.names=1)
dim(gen)

sum(is.na(gen)) # 27,987 NAs in the matrix (~3% missing data)

gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs

env <- read.csv("data/wolf_env.csv")
str(env) # Look at the structure of the data frame

env$individual <- as.character(env$individual) # Make individual names characters (not factors)
env$land_cover <- as.factor(env$land_cover)    # Make land cover a factor (not an integer)

# Confirm that genotypes and environmental data are in the same order
identical(rownames(gen.imp), env[,1]) 

pred <- subset(env, select=-c(precip_coldest_quarter, max_temp_warmest_month, min_temp_coldest_month))
pred$land_cover

pred <- subset(pred, select=-c(land_cover))

pred <- pred[,5:12]
colnames(pred) <- c("AMT","MDR","sdT","AP","cvP","NDVI","Elev","Tree")
pairs.panels(pred, scale=T)

wolf.rda <- rda(gen.imp ~ ., data=pred, scale=T)
wolf.rda

RsquareAdj(wolf.rda)

summary(eigenvals(wolf.rda, model = "constrained"))

screeplot(wolf.rda)

signif.full <- anova.cca(wolf.rda, parallel=getOption("mc.cores")) # default is permutation=999

signif.axis <- anova.cca(wolf.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis

# > signif.axis
# Permutation test for rda under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = gen.imp ~ AMT + MDR + sdT + AP + cvP + NDVI + Elev + Tree, data = pred, scale = T)
# Df Variance      F Pr(>F)    
# RDA1      1    281.2 3.3226  0.001 ***
# RDA2      1    216.9 2.5625  0.001 ***
# RDA3      1    179.5 2.1213  0.003 ** 
# RDA4      1    110.8 1.3097  0.076 .  
# RDA5      1     89.8 1.0609  0.870    
# RDA6      1     87.1 1.0298  0.876    
# RDA7      1     75.9 0.8971  0.996    
# RDA8      1     72.2 0.8532  0.911    
# Residual 85   7193.5                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


levels(env$ecotype) <- c("Western Forest","Boreal Forest","Arctic","High Arctic","British Columbia","Atlantic Forest")
eco <- env$ecotype
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c") # 6 nice colors for our ecotypes

# axes 1 & 2
plot(wolf.rda, type="n", scaling=3)
points(wolf.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(wolf.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the wolves
text(wolf.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

load.rda <- scores(wolf.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) # 38
cand2 <- outliers(load.rda[,2],3) # 69
cand3 <- outliers(load.rda[,3],3) # 34

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=8)  # 8 columns for 8 predictors
colnames(foo) <- c("AMT","MDR","sdT","AP","cvP","NDVI","Elev","Tree")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$snp[duplicated(cand$snp)])  # 7 duplicate detections

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #  7 duplicates on axis 2
table(foo[foo[,1]==3,2]) # no duplicates on axis 3
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,12] <- names(which.max(abs(bar[4:11]))) # gives the variable
  cand[i,13] <- max(abs(bar[4:11]))              # gives the correlation
}

colnames(cand)[12] <- "predictor"
colnames(cand)[13] <- "correlation"

table(cand$predictor) 

sel <- cand$snp
env <- cand$predictor
env[env=="AP"] <- '#1f78b4'
env[env=="cvP"] <- '#a6cee3'
env[env=="MDR"] <- '#6a3d9a'
env[env=="AMT"] <- '#e31a1c'
env[env=="NDVI"] <- '#33a02c'
env[env=="Elev"] <- '#ffff33'
env[env=="sdT"] <- '#fb9a99'
env[env=="Tree"] <- '#b2df8a'

# color by predictor:
col.pred <- rownames(wolf.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99','#b2df8a')


# axes 1 & 2
plot(wolf.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(wolf.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(wolf.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(wolf.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("AP","cvP","MDR","AMT","NDVI","Elev","sdT","Tree"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
