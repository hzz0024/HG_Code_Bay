library(dplyr)
outlier_1<-read.table("challenge_FDR_outlier.list", header=T)
head(outlier_1)
n1<-dim(outlier_1)[1]

outlier_2<-read.table("HC_NB_FDR_outlier.list", header=T)
head(outlier_2)
n2<-dim(outlier_2)[1]

outlier_temp_fulljoin<-full_join(outlier_1,outlier_2)
head(outlier_temp_fulljoin)
nALL<-dim(outlier_temp_fulljoin)[1]

outlier_temp_innerjoin<-intersect(outlier_1$id,outlier_2$id)
head(outlier_temp_innerjoin)
nboth<-length(outlier_temp_innerjoin)

library(ggVennDiagram)
ggVennDiagram(list(Del19_challenge = 1:n1, Del19_HC_NB = (n1+1-nboth):(n1-nboth+n2)), category.names = c("Surv. vs Ref","HC vs NB"))+ scale_fill_gradient(low="dimgray",high = "gainsboro")
graph2ppt(file="SGS",width=6,height=6)
############ SGS ############
outlier_1<-read.table("SGS_Del19.txt", header=T)
head(outlier_1)
n1<-dim(outlier_1)[1]

outlier_2<-read.table("SGS_HC_NB.txt", header=T)
head(outlier_2)
n2<-dim(outlier_2)[1]

outlier_3<-read.table("SGS_HC_SR.txt", header=T)
head(outlier_3)
n3<-dim(outlier_3)[1]

x = list( SGS_Surv19_Ref19 = outlier_1$ID,
          SGS_HC_NB = outlier_2$ID,
          SGS_HC_SR = outlier_3$ID)

ggVennDiagram(x, label_alpha = 0)


ggVennDiagram(list(Del2019 = 1:n1, CF_Surv19_Ref19_HC_SR = (n1+1-nboth):(n1-nboth+n2)), label_geom = geom_label)