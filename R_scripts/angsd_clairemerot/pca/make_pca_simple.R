#this R script performs a PCA on the covariance matrix and export it labelling rows with bamfile name and columns as PC
#it also makes a first visualisation in a pdf file with PC1-2 and PC3-4

argv <- commandArgs(T)
INPUT <- argv[1]
BAM <- argv[2]

INPUT = "All_0p7x_minq20_minmq30_CV30_masked_noinvers_shared_sites.covMat"
#perform a pca on covariance matrix
print(paste("read cov matrix", INPUT))
cov_mat<-as.matrix(read.table(INPUT), header=F)
pca<-eigen(cov_mat)
     
pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

#add column names
nPC<-dim(pca$vectors)[2]
col_PC<-vector(length=nPC)
for (i in 1 : nPC) {col_PC[i]<-paste0("PC",i)}
colnames(pca.mat)<-c(col_PC)

#add rownames
bam_names<-read.table(BAM,header=F)
rownames(pca.mat)<-bam_names$V1

#calculate varsum(eigen_mats$values[eigen_mats$values>=0]
var1<-round(pca$values[1]*100/sum(pca$values[pca$values>=0]),2)
var2<-round(pca$values[2]*100/sum(pca$values[pca$values>=0]),2)
var3<-round(pca$values[3]*100/sum(pca$values[pca$values>=0]),2)
var4<-round(pca$values[4]*100/sum(pca$values[pca$values>=0]),2)

#make kmeans for 3 groups on PC1
kmeans_res<-kmeans(as.matrix(pca.mat[,1]), c(min(pca.mat[,1]), median(pca.mat[,1]), max(pca.mat[,1])))
k_ss<-round(kmeans_res$betweenss/kmeans_res$totss,3)

#save 4PCS eigenvalues and k means SS
write.table(pca.mat[,1:4], paste0(INPUT,".pca"), quote=F)
write.table(c(var1,var2,var3,var4,k_ss), paste0(INPUT,".eig"), quote=F)

#plot pca
jpeg(file=paste0(INPUT,".pca.jpg"))
par(mfrow=c(1,1))
plot(pca.mat[,1], pca.mat[,2], pch=20, ylab=paste("PC2", var2), xlab=paste("PC1", var1),col=kmeans_res$cluster, main=paste("k_SS",k_ss))
#plot(pca.mat[,3], pca.mat[,4], pch=20, ylab=paste("PC4", var4), xlab=paste("PC3",var3))
dev.off()



#this R script performs a first visualisation in a pdf file with PC1-2 and PC3-4, and colours corresponding to groups
library(dplyr)
library(ade4)
argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERCENT_IND <- argv[2]


#read info file
x_info<-read.table("02_info/info.txt", header=T)
pca.mat<-read.table(paste0("04_pca/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"cov.pca"), header=T)

if (dim(x_info)[1]!=dim(pca.mat)[1]){print ("warning : not the same number of indivuals in the pca and the info file")}

#just in case order is different, use rownames of the pca (coming from bam.list) to order the info matrix
bam_name<-as.data.frame(rownames(pca.mat))
colnames(bam_name)<-colnames(x_info)[1]
bam_x_info<-left_join(bam_name,x_info)
head(bam_x_info)

#for each var, make a pdf with PC1-2 & PC3-4 coloured 
nvar<-dim(x_info)[2]-2 #we assume that col1 is bamfile name & col 2 is ind id

for (i in 1:nvar)
{
  #plot pca
  group<-as.factor(bam_x_info[,i+2])
  group_name<-colnames(bam_x_info)[i+2]
  pdf(file=paste0("04_pca/all_maf",MIN_MAF,"_pctind",PERCENT_IND,".pca.",group_name,".pdf"))
  par(mfrow=c(1,2))
  s.class(pca.mat, group,  xax=1, yax=2, cellipse=T,grid=F, col=as.numeric(group), cstar=F,sub="PC1-2", pch=20)
  s.class(pca.mat, group,  xax=3, yax=4, cellipse=T,grid=F, col=as.numeric(group), cstar=F,sub="PC3-4", pch=20)
  dev.off()
}

