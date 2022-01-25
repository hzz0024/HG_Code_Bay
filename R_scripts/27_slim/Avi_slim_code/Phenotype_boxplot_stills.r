library(ggplot2)
H1=read.csv("/Users/avi/Desktop/10MigXSel/HMHS_2_0516.csv", F)
H1=H1[1482:1782, ]
colnames(H1)=c("Gen","Pop","Par_fit","Juv_fit","Par_phe","Juv_phe","Co_A")

H2=read.csv("/Users/avi/Desktop/10MigXSel/HMMS_2_0516.csv", F)
H2=H2[1482:1782, ]
colnames(H2)=c("Gen","Pop","Par_fit","Juv_fit","Par_phe","Juv_phe","Co_A")

H3=read.csv("/Users/avi/Desktop/10MigXSel/HMLS_2_0516.csv", F)
H3=H3[1482:1782, ]
colnames(H3)=c("Gen","Pop","Par_fit","Juv_fit","Par_phe","Juv_phe","Co_A")

H4=read.csv("/Users/avi/Desktop/10MigXSel/MMHS_2_0516.csv", F)
H4=H4[1482:1782, ]
colnames(H4)=c("Gen","Pop","Par_fit","Juv_fit","Par_phe","Juv_phe","Co_A")

H5=read.csv("/Users/avi/Desktop/10MigXSel/MMMS_2_0516.csv", F)
H5=H5[1482:1782, ]
colnames(H5)=c("Gen","Pop","Par_fit","Juv_fit","Par_phe","Juv_phe","Co_A")

H6=read.csv("/Users/avi/Desktop/10MigXSel/MMLS_2_0516.csv", F)
H6=H6[1482:1782, ]
colnames(H6)=c("Gen","Pop","Par_fit","Juv_fit","Par_phe","Juv_phe","Co_A")

H7=read.csv("/Users/avi/Desktop/10MigXSel/LMHS_2_0516.csv", F)
H7=H7[1482:1782, ]
colnames(H7)=c("Gen","Pop","Par_fit","Juv_fit","Par_phe","Juv_phe","Co_A")

################################################################################

#H1

H1p=ggplot(H1, aes(x=Pop, y=Par_phe)) +
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#0099CC")

H1j=ggplot(H1A, aes(x=Pop, y=Juv_phe))+
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#66CC00")
################################################################################
#H2

H2p=ggplot(H2, aes(x=Pop, y=Par_phe)) +
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#0099CC")

H2j=ggplot(H2, aes(x=Pop, y=Juv_phe))+
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#66CC00")
################################################################################

#H3

H3p=ggplot(H3, aes(x=Pop, y=Par_phe)) +
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#0099CC")

H3j=ggplot(H3, aes(x=Pop, y=Juv_phe))+
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#66CC00")
################################################################################

#H4

H4p=ggplot(H4, aes(x=Pop, y=Par_phe)) +
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#0099CC")

H4j=ggplot(H4, aes(x=Pop, y=Juv_phe))+
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#66CC00")

################################################################################

#H5

H5p=ggplot(H5, aes(x=Pop, y=Par_phe)) +
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#0099CC")

H5j=ggplot(H5, aes(x=Pop, y=Juv_phe))+
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#66CC00")
#######################################################################
#H6

H6p=ggplot(H6, aes(x=Pop, y=Par_phe)) +
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#0099CC")

H6j=ggplot(H6, aes(x=Pop, y=Juv_phe))+
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#66CC00")
################################################################################




#H7

H7p=ggplot(H7, aes(x=Pop, y=Par_phe)) +
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#0099CC")

H7j=ggplot(H7, aes(x=Pop, y=Juv_phe))+
  theme(legend.position="right") +
  ylim(c(0,1)) +
  ylab("Phenotype")+
  xlab("Population")+
  theme(text = element_text(size=40),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#66CC00")

