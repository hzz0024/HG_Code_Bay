########################################
##### Extract the pruned snp list ######
########################################

full_list = read.delim("Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", header = FALSE, sep='\t')
pruned_list = read.delim("WILD_LD_prunning_snp_n_563625.list", header = FALSE, sep=':')
full_list$SNP = paste0(full_list$V1,'_',full_list$V2)
pruned_list$SNP = paste0(pruned_list$V1,'_',pruned_list$V2)
head(full_list)
head(pruned_list)
pruned_final = full_list[which(full_list$SNP %in% pruned_list$SNP),][,1:4]
head(pruned_final)
dim(pruned_final)
write.table(pruned_final, "Del19_pruned_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", row.names=F, col.names = F, quote=F, sep="\t")

########################
##### Mantel test ######
########################

library(ade4)
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/21_IBD")
m1 <- read.table("purned_fst.csv",header=TRUE,sep=",")
m1 <- read.table("global_fst.csv",header=TRUE,sep=",")
m1 <- read.table("SGS_fst.csv",header=TRUE,sep=",")
m1 <- read.table("fisher_fst.csv",header=TRUE,sep=",")
m1 <- read.table("genome_no_SGS_outlier_fst.csv",header=TRUE,sep=",")

m1 <- read.table("global_fst.csv",header=TRUE,sep=",")
m2 <- read.table("matrix.csv",header=TRUE,sep="\t")

gen <- quasieuclid(as.dist(m1))
geo <- quasieuclid(as.dist(m2))
plot(r1 <- mantel.randtest(geo,gen, nrepet = 10000), main = "Mantel test")
r1

#####################
##### IBD plot ######
#####################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/21_IBD")
# https://boulter.com/gps/distance/?from=39.2935%09-75.3435&to=39.248%09-75.248&units=k
library(ggplot2)
library(cowplot)
library(ggrepel)
library(export)
library(ggpmisc)
IBD <- read.csv("IBD_dis_plot_pruned.csv", header = TRUE)

IBD$Fst

cbPalette <- c("#A71B4B", "#E97302", "#EAC728", "#0BC0B3", "#4461A8", "#F5191C", "#7A7A7A")

my.formula <- y ~ x
ggplot(data = IBD, aes(x = Distance, y = Fst)) +
  geom_point(aes(color=factor(pop2)), shape=20,  size=8)+
  scale_colour_manual(values=cbPalette, breaks=c("ARN", "COH", "SR", "NB"))+
  geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
  theme_bw() +
  theme(axis.title.y=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  scale_y_continuous(name="Fst/(1-Fst)") +
  scale_x_continuous(name="Distance(km)", limits=c(0, 32))+
  geom_smooth(method = "lm", se = TRUE, fullrange=TRUE)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) 

graph2ppt(file="Wild_IBD_pruned1", width=8, height=6)   

###############################
# pruned without NB population#
###############################

IBD <- read.csv("IBD_dis_plot_pruned_noNB.csv", header = TRUE)
IBD$Fst

cbPalette <- c("#A71B4B", "#E97302", "#EAC728", "#0BC0B3", "#4461A8", "#F5191C", "#7A7A7A")

my.formula <- y ~ x
ggplot(data = IBD, aes(x = Distance, y = Fst)) +
  geom_point(aes(color=factor(pop2)), shape=20,  size=8)+
  scale_colour_manual(values=cbPalette, breaks=c("HC", "ARN", "COH", "SR"))+
  geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
  theme_bw() +
  theme(axis.title.y=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  scale_y_continuous(name="Fst/(1-Fst)") +
  scale_x_continuous(name="Distance(km)", limits=c(0, 32))+
  geom_smooth(method = "lm", se = TRUE, fullrange=TRUE)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) 

graph2ppt(file="Wild_IBD_pruned", width=8, height=6)

###################################################
# combine the pruned snps with outlier candidates #
###################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/21_IBD")
library(lmodel2)

IBD <- read.csv("IBD_dis_plot_combined_all.csv", header = TRUE)

whole = IBD[which(IBD$group == "Whole"),]
neutral = IBD[which(IBD$group == "Neutral"),]
Fisher = IBD[which(IBD$group == "Fisher"),]
SGS_wild = IBD[which(IBD$group == "SGS_wild"),]
SGS_challenge = IBD[which(IBD$group == "SGS_challenge"),]
SGS_both = IBD[which(IBD$group == "SGS_both"),]

mod_whole = lmodel2(Fst ~ Distance, data=whole, "interval", "interval", 95)
mod_neutral = lmodel2(Fst ~ Distance, data=neutral, "interval", "interval", 95)
mod_Fisher = lmodel2(Fst ~ Distance, data=Fisher, "interval", "interval", 95)
mod_SGS_wild = lmodel2(Fst ~ Distance, data=SGS_wild, "interval", "interval", 95)
mod_SGS_challenge = lmodel2(Fst ~ Distance, data=SGS_challenge, "interval", "interval", 95)
mod_SGS_both = lmodel2(Fst ~ Distance, data=SGS_both, "interval", "interval", 95)

reg_whole = mod_whole$regression.results[which(mod_whole$regression.results$Method == "SMA"),]
reg_neutral = mod_neutral$regression.results[which(mod_neutral$regression.results$Method == "SMA"),]
reg_Fisher = mod_Fisher$regression.results[which(mod_Fisher$regression.results$Method == "SMA"),]
reg_SGS_wild = mod_SGS_wild$regression.results[which(mod_SGS_wild$regression.results$Method == "SMA"),]
reg_SGS_challenge = mod_SGS_challenge$regression.results[which(mod_SGS_challenge$regression.results$Method == "SMA"),]
reg_SGS_both = mod_SGS_both$regression.results[which(mod_SGS_both$regression.results$Method == "SMA"),]

names(reg_whole) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_neutral) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_Fisher) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_SGS_wild) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_SGS_challenge) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_SGS_both) = c("method", "intercept", "slope", "angle", "p-value")

cbPalette <- c("#A71B4B", "#E97302", "#EAC728", "#0BC0B3", "#4461A8", "#F5191C", "#7A7A7A")
my.formula <- y ~ x

ggplot(data = IBD, aes(x = Distance, y = Fst, color=group, shape=group)) +
  geom_point(size=5)+
  #scale_colour_manual(values=cbPalette, breaks=c("HC", "ARN", "COH", "SR", "NB"))+
  #geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
  theme_bw() +
  #geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
  theme(axis.title.y=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  scale_y_continuous(name="Fst/(1-Fst)") +
  scale_x_continuous(name="Distance(km)", limits=c(0, 32))+
  cowplot::theme_cowplot()+
  theme(panel.grid.major = element_line(color = "lightgrey",
                                        size = 0.5,
                                        linetype = "dashed"))+
  geom_abline(data = reg_whole, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_neutral, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_Fisher, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_SGS_wild, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_SGS_challenge, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_SGS_both, aes(intercept = intercept, slope = slope))+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~~")), 
               parse = TRUE)  # need to replace the lm intercept and slope with RMA ones. Change them in ppt.


# plot the p-values of mantel test
ps <- read.csv("p_values.csv", header = TRUE)
ggplot(ps, aes(x=group, y=ps, color=group)) + 
  geom_violin(trim=FALSE, fill="gray")+
  labs(title="p-values from Mantel tests (10,000 permutations)",x=NULL, y = "p-values")+
  geom_boxplot(width=0.1)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4, color="red")+
  theme_bw() +
  theme(axis.title.y=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  cowplot::theme_cowplot()+
  theme(panel.grid.major = element_line(color = "lightgrey",
                                        size = 0.5,
                                        linetype = "dashed"))
graph2ppt(file="ps", width=6, height=6)

# my.formula <- y ~ x
# ggplot(data = IBD, aes(x = Distance, y = Fst, color=group, shape=group)) +
#   geom_point(size=5)+
#   #scale_colour_manual(values=cbPalette, breaks=c("HC", "ARN", "COH", "SR", "NB"))+
#   #geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
#   theme_bw() +
#   theme(axis.title.y=element_text(size=15))+
#   theme(axis.title.x=element_text(size=15))+
#   theme(axis.text=element_text(size=15))+
#   scale_y_continuous(name="Fst/(1-Fst)") +
#   scale_x_continuous(name="Distance(km)", limits=c(0, 32))+
#   geom_smooth(aes(group=group),
#               method = "lm", se = TRUE, fullrange=TRUE)+
#   cowplot::theme_cowplot()+
#   theme(panel.grid.major = element_line(color = "lightgrey",
#                                         size = 0.5,
#                                         linetype = "dashed"))+
#   facet_grid(group ~ ., scales = 'free_y') +
#   stat_poly_eq(formula = my.formula, 
#                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                parse = TRUE) 



graph2ppt(file="IBD_summary_whole_vs_neutral", width=6, height=6)
graph2ppt(file="IBD_summary_nofree_y", width=6, height=6)
graph2ppt(file="IBD_whole", width=8, height=4)


###################################################
# combine the outlier candidates #
###################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/21_IBD")
library(lmodel2)

IBD <- read.csv("IBD_dis_plot_combined_outlier.csv", header = TRUE)

Fisher = IBD[which(IBD$group == "Fisher"),]
SGS_wild = IBD[which(IBD$group == "SGS_wild"),]
SGS_challenge = IBD[which(IBD$group == "SGS_challenge"),]
SGS_both = IBD[which(IBD$group == "SGS_both"),]

mod_Fisher = lmodel2(Fst ~ Distance, data=Fisher, "interval", "interval", 95)
mod_SGS_wild = lmodel2(Fst ~ Distance, data=SGS_wild, "interval", "interval", 95)
mod_SGS_challenge = lmodel2(Fst ~ Distance, data=SGS_challenge, "interval", "interval", 95)
mod_SGS_both = lmodel2(Fst ~ Distance, data=SGS_both, "interval", "interval", 95)

reg_Fisher = mod_Fisher$regression.results[which(mod_Fisher$regression.results$Method == "SMA"),]
reg_SGS_wild = mod_SGS_wild$regression.results[which(mod_SGS_wild$regression.results$Method == "SMA"),]
reg_SGS_challenge = mod_SGS_challenge$regression.results[which(mod_SGS_challenge$regression.results$Method == "SMA"),]
reg_SGS_both = mod_SGS_both$regression.results[which(mod_SGS_both$regression.results$Method == "SMA"),]

names(reg_Fisher) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_SGS_wild) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_SGS_challenge) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_SGS_both) = c("method", "intercept", "slope", "angle", "p-value")

cbPalette <- c("#A71B4B", "#E97302", "#EAC728", "#0BC0B3", "#4461A8", "#F5191C", "#7A7A7A")
my.formula <- y ~ x

ggplot(data = IBD, aes(x = Distance, y = Fst, color=group, shape=group)) +
  geom_point(size=5)+
  #scale_colour_manual(values=cbPalette, breaks=c("HC", "ARN", "COH", "SR", "NB"))+
  #geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
  theme_bw() +
  #geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
  theme(axis.title.y=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  scale_y_continuous(name="Fst/(1-Fst)") +
  scale_x_continuous(name="Distance(km)", limits=c(0, 32))+
  cowplot::theme_cowplot()+
  theme(panel.grid.major = element_line(color = "lightgrey",
                                        size = 0.5,
                                        linetype = "dashed"))+
  geom_abline(data = reg_Fisher, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_SGS_wild, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_SGS_challenge, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_SGS_both, aes(intercept = intercept, slope = slope))+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~~")), 
               parse = TRUE)  # need to replace the lm intercept and slope with RMA ones. Change them in ppt.

graph2ppt(file="IBD_outlier", width=8, height=4)

###################################################
# combine the pruned snps with outlier candidates #
###################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/21_IBD")
library(lmodel2)

IBD <- read.csv("IBD_dis_plot_combined_whole_neutral_challenge.csv", header = TRUE)

whole = IBD[which(IBD$group == "Whole"),]
neutral = IBD[which(IBD$group == "Neutral"),]
SGS_challenge = IBD[which(IBD$group == "SGS_challenge"),]

mod_whole = lmodel2(Fst ~ Distance, data=whole, "interval", "interval", 95)
mod_neutral = lmodel2(Fst ~ Distance, data=neutral, "interval", "interval", 95)
mod_SGS_challenge = lmodel2(Fst ~ Distance, data=SGS_challenge, "interval", "interval", 95)

reg_whole = mod_whole$regression.results[which(mod_whole$regression.results$Method == "SMA"),]
reg_neutral = mod_neutral$regression.results[which(mod_neutral$regression.results$Method == "SMA"),]
reg_SGS_challenge = mod_SGS_challenge$regression.results[which(mod_SGS_challenge$regression.results$Method == "SMA"),]

names(reg_whole) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_neutral) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_SGS_challenge) = c("method", "intercept", "slope", "angle", "p-value")

cbPalette <- c("#A71B4B", "#E97302", "#EAC728", "#0BC0B3", "#4461A8", "#F5191C", "#7A7A7A")
my.formula <- y ~ x

ggplot(data = IBD, aes(x = Distance, y = Fst, color=group, shape=group)) +
  geom_point(size=5)+
  #scale_colour_manual(values=cbPalette, breaks=c("HC", "ARN", "COH", "SR", "NB"))+
  #geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
  theme_bw() +
  #geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
  theme(axis.title.y=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  scale_y_continuous(name="Fst/(1-Fst)") +
  scale_x_continuous(name="Distance(km)", limits=c(0, 32))+
  cowplot::theme_cowplot()+
  theme(panel.grid.major = element_line(color = "lightgrey",
                                        size = 0.5,
                                        linetype = "dashed"))+
  geom_abline(data = reg_whole, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_neutral, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_Fisher, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_SGS_wild, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_SGS_challenge, aes(intercept = intercept, slope = slope))+
  geom_abline(data = reg_SGS_both, aes(intercept = intercept, slope = slope))+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~~")), 
               parse = TRUE)  # need to replace the lm intercept and slope with RMA ones. Change them in ppt.


# plot the p-values of mantel test
ps <- read.csv("p_values.csv", header = TRUE)
ggplot(ps, aes(x=group, y=ps, color=group)) + 
  geom_violin(trim=FALSE, fill="gray")+
  labs(title="p-values from Mantel tests (10,000 permutations)",x=NULL, y = "p-values")+
  geom_boxplot(width=0.1)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4, color="red")+
  theme_bw() +
  theme(axis.title.y=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  cowplot::theme_cowplot()+
  theme(panel.grid.major = element_line(color = "lightgrey",
                                        size = 0.5,
                                        linetype = "dashed"))
graph2ppt(file="ps", width=6, height=6)

# my.formula <- y ~ x
# ggplot(data = IBD, aes(x = Distance, y = Fst, color=group, shape=group)) +
#   geom_point(size=5)+
#   #scale_colour_manual(values=cbPalette, breaks=c("HC", "ARN", "COH", "SR", "NB"))+
#   #geom_text_repel(data = IBD, mapping=aes(x = Distance, y = Fst, label=Fst), box.padding = unit(0.5, "point"))+
#   theme_bw() +
#   theme(axis.title.y=element_text(size=15))+
#   theme(axis.title.x=element_text(size=15))+
#   theme(axis.text=element_text(size=15))+
#   scale_y_continuous(name="Fst/(1-Fst)") +
#   scale_x_continuous(name="Distance(km)", limits=c(0, 32))+
#   geom_smooth(aes(group=group),
#               method = "lm", se = TRUE, fullrange=TRUE)+
#   cowplot::theme_cowplot()+
#   theme(panel.grid.major = element_line(color = "lightgrey",
#                                         size = 0.5,
#                                         linetype = "dashed"))+
#   facet_grid(group ~ ., scales = 'free_y') +
#   stat_poly_eq(formula = my.formula, 
#                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                parse = TRUE) 



graph2ppt(file="IBD_summary_whole_vs_neutral_vs_challenge", width=8, height=4)
graph2ppt(file="IBD_summary_nofree_y", width=6, height=6)
graph2ppt(file="IBD_whole", width=8, height=4)