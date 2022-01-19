
########################
##### Mantel test ######
########################
library(ade4)

setwd("~/Dropbox/CVreseq/IBD/IBD_input")
m1 <- read.table("neutral_fst_matrix.txt",header=TRUE,sep="\t")
m2 <- read.table("aquatic_dis_matrix.txt",header=TRUE,sep="\t")

gen <- quasieuclid(as.dist(m1))
geo <- quasieuclid(as.dist(m2))
plot(r1 <- mantel.randtest(geo,gen, nrepet = 10000), main = "Mantel test")
r1
# Monte-Carlo test
# Call: mantel.randtest(m1 = geo, m2 = gen, nrepet = 10000)
# 
# Observation: 0.5004774 
# 
# Based on 10000 replicates
# Simulated p-value: 0.01059894 
# Alternative hypothesis: greater 
# 
# Std.Obs   Expectation      Variance 
# 2.0880824275 -0.0002230122  0.0574990352 

#####################
##### IBD plot ######
#####################
setwd("~/Dropbox/CVreseq/IBD/IBD_input")
library(ggplot2)
library(cowplot)
library(ggrepel)
library(export)
library(ggpmisc)
library(lmodel2)

IBD <- read.csv("fst_neutral.csv", header = TRUE)

IBD$Fst
my.formula <- y ~ x
neutral = IBD[which(IBD$group == "Neutral"),]
mod_neu = lmodel2(Fst ~ Distance, data=neutral, "interval", "interval", 95)
reg_neu = mod_neu$regression.results[which(mod_neu$regression.results$Method == "SMA"),]
names(reg_neu) = c("method", "intercept", "slope", "angle", "p-value")
cbPalette <- c("#A71B4B", "#E97302", "#EAC728", "#0BC0B3", "#4461A8", "#F5191C", "#7A7A7A")
my.formula <- y ~ x
reg_neu
ggplot(data = IBD, aes(x = Distance, y = Fst, color=group, shape=group)) +
  geom_point(size=5)+
  #scale_colour_manual(values=cbPalette, breaks=c("CLP", "CS", "HC", "HCVA", "HI"))+
  #geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
  theme_bw() +
  theme(axis.title.y=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  scale_y_continuous(name="Fst/(1-Fst)") +
  cowplot::theme_cowplot()+
  theme(panel.grid.major = element_line(color = "lightgrey",
                                        size = 0.5,
                                        linetype = "dashed"))+
  geom_abline(data = reg_neu, aes(intercept = intercept, slope = slope))

graph2ppt(file="SM_IBD", width=8, height=6) 




ggplot(data = IBD, aes(x = Distance, y = Fst)) +
  geom_point(aes(color=factor(pop2)), shape=20,  size=8)+
  scale_colour_manual(values=cbPalette, breaks=c("CLP", "CS", "HC", "HC_VA", "HI", "SM"))+
  geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
  theme_bw() +
  theme(axis.title.y=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  scale_y_continuous(name="Fst/(1-Fst)") +
  scale_x_continuous(name="Distance(km)", limits=c(0, 1500))+
  geom_smooth(method = "lm", se = TRUE, fullrange=TRUE)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) 


