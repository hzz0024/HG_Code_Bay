H1=read.csv("/Users/avi/Desktop/10MigXSel/HMHS_2_0516.csv", F)
H1=H1[1482:1782, ]
colnames(H1)=c("Gen","Pop","Par_fit","Juv_fit","Par_phe","Juv_phe","Co_A")
H1$Pop=as.factor(H1$Pop)

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

H8=read.csv("/Users/avi/Desktop/10MigXSel/LMMS_2_0516.csv", F)
H8=H8[1482:1782, ]
colnames(H8)=c("Gen","Pop","Par_fit","Juv_fit","Par_phe","Juv_phe","Co_A")

H9=read.csv("/Users/avi/Desktop/10MigXSel/LMLS_2_0516.csv", F)
H9=H9[1482:1782, ]
colnames(H9)=c("Gen","Pop","Par_fit","Juv_fit","Par_phe","Juv_phe","Co_A")

ggplot(H2, aes(x=Pop, y=Juv_phe)) +
  theme(legend.position="none")+
  ylim(c(0,1))+
  xlab("Population")+
  ylab("Phenotype")+
  theme(axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="#ff66ff")
  
Calypso=ggplot(H2, aes(x=Pop, y=Par_phe)) +
  theme(legend.position="none")+
  ylim(c(0,1))+
  xlab("Population")+
  ylab("Phenotype")+
  theme(axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"),
        axis.text.x = element_text(colour = "black", face="bold"),
        axis.text.y = element_text(colour = "black", face="bold"))+
  geom_boxplot(color="blue")

#Next steps: save both r plots and code animation to flip between them 
  
  