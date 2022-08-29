
library(corrplot)
library(ade4)
library(export)
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/00_env_pca_cor/")
#env_16pop<-read.table("selection/02_raw_data/env_16pop.txt", header=T)
env_18pop<-read.table("info_pop_env.txt", header=T,fill=TRUE)

PCA_pearson(env_18pop)
graph2ppt(file="Env_cor_pop18",width=10,height=10)

var_clim<-env_18pop[,-1]
PCA_subset(var_subset=var_clim, categ="clim", out_folder="pca_18pop")

env_16pop<-read.table("evo_16pop.txt", header=T,fill=TRUE)
env_16pop
var_clim<-env_16pop[,-1]
PCA_pearson(env_16pop)
PCA_subset(var_subset=var_clim, categ="clim", out_folder="env_16pop")

var_clim<-cbind(env_16pop$GpsNS,env_16pop$MeanSST,env_16pop$MeanTemp, env_16pop$MeanPrec)
colnames(var_clim)<-c("GPSNS","MeanSST", "MeanTemp", "MeanPrec")
PCA_subset(var_subset=var_clim, categ="clim", out_folder="")

var_clim2<-cbind(env_16pop$MeanSST,env_16pop$MeanTemp, env_16pop$MeanPrec)
colnames(var_clim2)<-c("MeanSST", "MeanTemp", "MeanPrec")
PCA_subset(var_subset=var_clim2, categ="clim2", out_folder="selection/02_raw_data/pca_env_var/")

var_temp<-cbind(env_16pop$GpsNS,env_16pop$MeanSST,env_16pop$MeanTemp)
colnames(var_temp)<-c("GPSNS","MeanSST", "MeanTemp")
PCA_subset(var_subset=var_temp, categ="temp", out_folder="selection/02_raw_data/pca_env_var/")

var_sali<-cbind(env_16pop$GpsEW,env_16pop$MeanSSS,env_16pop$MeanTA)
colnames(var_sali)<-c("GPSEW","MeanSSS", "MeanTA")
PCA_subset(var_subset=var_sali, categ="sali", out_folder="selection/02_raw_data/pca_env_var/")

var_sata<-cbind(env_16pop$MeanSSS,env_16pop$MeanTA)
colnames(var_sata)<-c("MeanSSS", "MeanTA")
PCA_subset(var_subset=var_sata, categ="sata", out_folder="selection/02_raw_data/pca_env_var/")

var_lafu<-cbind(env_16pop$Lamin,env_16pop$Fucus)
colnames(var_lafu)<-c("Lamin","Fucus")
PCA_subset(var_subset=var_lafu, categ="lafu", out_folder="selection/02_raw_data/pca_env_var/")

var_comp<-cbind(env_16pop$Lamin,env_16pop$Fucus, env_16pop$OtherSeaweed, env_16pop$Zoostera, env_16pop$OtherMat)
colnames(var_comp)<-c("Lamin","Fucus", "OS","zoo","OM")
PCA_subset(var_subset=var_comp, categ="comp", out_folder="selection/02_raw_data/pca_env_var/")

var_beda<-cbind(env_16pop$BedD,env_16pop$BedS,env_16pop$BedT)
colnames(var_beda)<-c("BedD","BedS", "BedT")
PCA_subset(var_subset=var_beda, categ="beda", out_folder="selection/02_raw_data/pca_env_var/")

var_best<-cbind(env_16pop$BedS,env_16pop$BedT)
colnames(var_best)<-c("BedS", "BedT")
PCA_subset(var_subset=var_best, categ="best", out_folder="selection/02_raw_data/pca_env_var/")



#useful functions

PCA_subset<-function(var_subset, categ, out_folder)
{
  x_noNA<-var_subset
  na_posit<-which(is.na(x_noNA), TRUE)
  if (dim(na_posit)[1]>=1)
  {
    for (i in 1:dim(na_posit)[1] )
    {
      x_noNA[na_posit[i,1],na_posit[i,2]]<-mean(var_subset[,na_posit[i,2]], na.rm=T)
    }}
  pca_subset<-dudi.pca(x_noNA, scannf=FALSE, nf=20)
  percent_PC1<-round(pca_subset$eig[1]/sum(pca_subset$eig),2)*100
  percent_PC2<-round(pca_subset$eig[2]/sum(pca_subset$eig),2)*100
  jpeg( paste0(out_folder, "pca_", categ, ".jpg"))
  s.arrow(pca_subset$co, xax=1, yax=2, sub= paste("PC1", percent_PC1,"%", " PC2",percent_PC2,"%"))
  dev.off()
  evplot(pca_subset$eig)
  s.arrow(pca_subset$co, xax=1, yax=2, sub= paste("PC1", percent_PC1,"%", " PC2",percent_PC2,"%"), grid=F)
  
  pca_var<-as.data.frame(pca_subset$li[,1])
  colnames(pca_var)<-c(paste0("pca_",categ))#, paste0("pca_",categ,"_2"))
  head(pca_var)
  
  write.table(pca_var,  paste0(out_folder, "pca_", categ, ".txt"), row.names=FALSE, quote=FALSE)
  write.table(pca_subset$co, paste(out_folder, "pca_co_", categ, ".txt", sep=""), quote=FALSE)
  if (percent_PC1<60)
  {
    pca_var<-as.data.frame(pca_subset$li[,1:2])
    colnames(pca_var)<-c(paste0("pca_",categ), paste0("pca_",categ,"_2"))
    head(pca_var)
    write.table(pca_var, paste(out_folder, "pca_", categ, ".txt", sep=""), row.names=FALSE, quote=FALSE)
  }
}


PCA_pearson <- function(x)
{
  ##display variables on a pca
  x_noNA<-x
  na_posit<-which(is.na(x_noNA), TRUE)
  if (dim(na_posit)[1]>=1)
  {
    for (i in 1:dim(na_posit)[1] )
    {
      x_noNA[na_posit[i,1],na_posit[i,2]]<-mean(x[,na_posit[i,2]], na.rm=T)
    }
  }
  data_var_cs_noNA<<-cbind(x[,1],scale(x_noNA[,2:dim(x)[2]])) #center and scale the matrix
  var_pca<<-dudi.pca(x_noNA[,2:dim(x)[2]], scannf=FALSE, nf=10)
  percent_PC1<-round(var_pca$eig[1]/sum(var_pca$eig),2)*100
  percent_PC2<-round(var_pca$eig[2]/sum(var_pca$eig),2)*100
  
  s.class(var_pca$li, xax=1, yax=2, fac=as.factor(x[,1]))
  cat(paste("PC1", percent_PC1,"%"),paste("PC2",percent_PC2,"%"))
  s.arrow(var_pca$co, xax=1, yax=2)
  evplot(var_pca$eig)
  
  ##test pearson correlation
  variable<-c(2:dim(x)[2])
  c2<-combn(variable,2)
  cor_info<-matrix(ncol=4)
  cor_variables<-matrix(ncol=4)
  cor_all<-matrix(ncol=4)
  for (j in 1: dim (c2)[2])
  {
    pearson<-cor.test(x[,c2[1,j]],x[,c2[2,j]])
    cor_info<-rbind(cor_info,c((colnames(x)[c2[1,j]]),(colnames(x)[c2[2,j]]), pearson$p.value, pearson$estimate))
    cor_all<-rbind(cor_all, c((colnames(x)[c2[1,j]]),(colnames(x)[c2[2,j]]), pearson$p.value, pearson$estimate) )
    if (pearson$p.value<=0.05){cor_variables<-rbind(cor_variables, c((colnames(x)[c2[1,j]]),(colnames(x)[c2[2,j]]), pearson$p.value, pearson$estimate) )}
    
  }
  
  colnames(cor_variables)<-c("var1", "var2", "p_value", "cor")
  colnames(cor_all)<-c("var1", "var2", "p_value", "cor")
  print(cor_variables[2:dim(cor_variables)[1],])
  
  cor_all<-as.data.frame(cor_all[2:dim(cor_all)[1],])
  print
  cor_all$p_value<-as.numeric(as.character(cor_all$p_value))
  L<-length(cor_all$p_value)
  w=which(sort(cor_all$p_value)<(0.05*c(1:L)/L))
  plim=sort(cor_all$p_value)[length(w)]
  print(plim)
  # Benjamini–Hochberg or FDR correction
  
  cormat<<-(signif(cor(x[, 2: (dim(x)[2])]),2))
  p_cor<-cor.test.p (x[, 2: (dim(x)[2])])
  col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                             "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                             "#D6604D", "#B2182B", "#67001F")) (200)
  #corrplot(abs(cormat), diag=F, type="lower", is.corr=F,col=col3, tl.col=1,cl.pos="b", cl.lim=c(0,1))
  corrplot((cormat),   method="square",col=col2, tl.col=1, cl.pos="b", p.mat=p_cor, insig="label_sig", pch.cex=1, sig.level = plim)
}

evplot <- function(ev)
{
  # Broken stick model (MacArthur 1957)
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

