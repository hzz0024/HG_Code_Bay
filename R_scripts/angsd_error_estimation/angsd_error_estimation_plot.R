library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)

a <- read.delim("Del19_error.error.chunkunordered", header = FALSE, sep='\t')
a$V17 <- NULL
a <- as.vector(colMeans(a, dim=1))
b <- read.delim("Del20_error.error.chunkunordered", header = FALSE, sep='\t')
b$V17 <- NULL
b <- as.vector(colMeans(b, dim=1))

meacDNA<-c("A","C","G","T")
type <- as.vector((sapply(1:4,function(x) paste(DNA[x],DNA,sep="-to-"))))
df <- data.frame(Base_substitution=rep(type, 2), Frequency=c(a, b), Batch=rep(c("2019", "2020"), each=16))
#barplot(a)
#barplot(b)
ggplot(df, aes(x=Base_substitution, y=Frequency, group=Batch)) +
  geom_point(aes(shape=Batch, color=Batch), size=5)+
  geom_line(aes(color=Batch), size=2)+
  theme(axis.text.x = element_text(color = "grey20", size = 24, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 24, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 26, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 26, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(legend.title = element_text(colour="grey20", size=26))+
  theme(legend.text = element_text(colour="grey20", size=24))

graph2ppt(file="Base_substitution",width=18,height=9)