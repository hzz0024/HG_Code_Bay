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


########################################
##### Extract the GEA snp list ######
########################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/06_pairwise_fst_IBD/00_extract_snps")
snp_table <- read.table('All_sfs_all_sites_noparalogs.snplist.txt',  sep = '\t', col.names = c('Chromosome', 'Position', 'Major', 'minor'))
#genomic_ranges <- read.table('CD5_lfmm_rda_outliers_n_115_10K.txt', sep = '\t', col.names = c('Chromosome', 'Start', 'End'))
genomic_ranges <- read.table('MAX10_lfmm_rda_outliers_n_45_10K_list.txt', sep = '\t', col.names = c('Chromosome', 'Start', 'End'))
snp_list = data.frame()
for (i in 1:(dim(genomic_ranges)[1])){
  range = genomic_ranges[i,]
  snp_list_ <- snp_table[snp_table$Chromosome == range$Chromosome & snp_table$Position > as.numeric(range$Start) & snp_table$Position < as.numeric(range$End),]
  snp_list <- rbind(snp_list, snp_list_)
}

how_many_in_range <- function(coords){
  coords = as.vector(coords)
  sum(snp_table$Chromosome == coords[1] & snp_table$Position > as.numeric(coords[2]) & snp_table$Position < as.numeric(coords[3]))
}

genomic_ranges$number_of_snps <- apply(genomic_ranges, 1, how_many_in_range)

write.table(snp_list[,1:2], "MAX10_lfmm_rda_outliers_n_45_10K_snp_list.txt", row.names=F, col.names = F, quote=F, sep="\t")

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
m1 <- read.table("SGS_wild_fst.csv",header=TRUE,sep=",")
m2 <- read.table("matrix.csv",header=TRUE,sep=",")

gen <- quasieuclid(as.dist(m1))
geo <- quasieuclid(as.dist(m2))
plot(r1 <- mantel.randtest(geo,gen, nrepet = 10000), main = "Mantel test")
r1

##########################
# smatr package tutorial #
##########################

# Load leaf lifetime dataset:
data(leaflife)

### One sample analyses ###
# Extract only low-nutrient, low-rainfall data:
leaf.low <- subset(leaflife, soilp == 'low' & rain == 'low')

# Fit a MA for log(leaf longevity) vs log(leaf mass per area):
ma(longev ~ lma, log='xy', data=leaflife)

# Test if the MA slope is not significantly different from 1:
ma.test <- ma(longev ~ lma, log='xy', slope.test=1, data=leaflife)
summary(ma.test)

# Construct a residual plot to check assumptions:
plot(ma.test,type="residual")

### Several sample analyses ###

# Now consider low-nutrient sites (high and low rainfall):
leaf.low.soilp <- subset(leaflife, soilp == 'low')

# Fit SMA's separately at each of high and low rainfall sites,
# and test for common slope:
com.test <- sma(longev~lma*rain, log="xy", data=leaf.low.soilp)
com.test

# Plot longevity vs LMA separately for each group:
plot(com.test)

# Fit SMA's separately at each of high and low rainfall sites,
# and test if there is a common slope equal to 1:
sma(longev~lma*rain, log="xy", slope.test=1, data=leaf.low.soilp)

# Fit SMA's with common slope across each of high and low rainfall sites, 
# and test for common elevation:
sma(longev~lma+rain, log="xy", data=leaf.low.soilp)

# Fit SMA's with common slope across each of high and low rainfall sites, 
# and test for no shift along common SMA:
sma(longev~lma+rain, log="xy", type="shift", data=leaf.low.soilp)

Examples
# Load leaf lifetime dataset:
data(leaflife)

# Only consider low-nutrient sites:
leaf.low.soilp <- subset(leaflife, soilp == 'low')

# Fit SMA's separately at each of high and low 
# rainfall sites and test for common slope:
ft <- sma(longev~lma*rain, data=leaf.low.soilp, log="xy")

# Plot leaf longevity (longev) vs leaf mass per area (lma) 
# separately for each of high and low rainfall:
plot(ft)

# As before but add lines which have a common slope:
plot(ft, use.null=TRUE)

#As above, but adding the common slope lines to an existing plot
plot(ft, type='p', col="black")
plot(ft, use.null=TRUE, add=TRUE, type='l')

# Plot with equally spaced tick marks:
plot(ft, xaxis=defineAxis(major.ticks=c(40,80,160,320,640)), 
     yaxis=defineAxis(major.ticks=c(0.5,1,2,4,8)) )

# Produce a residual plot to check assumptions:
plot(ft,which="res")

# Produce a normal quantile plot:
plot(ft,which="qq")

###############################
# plot for outlier candidates #
###############################
rm(list=ls())
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/06_pairwise_fst_IBD/01_GEA_outlier_fst/")
library(lmodel2)
library(cowplot)
library(ggrepel)
library(export)
library(ggpmisc)
library(smatr)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

IBD <- read.csv("IBD_dis_plot_GAE.csv", header = TRUE)

GEA_CD5 = IBD[which(IBD$Group == "GEA_CD5"),]
neutral = IBD[which(IBD$Group == "Neutral"),]

# calculate for slope and intercept
GEA_res <- sma(Fst ~ Distance, data=GEA_CD5, method = "SMA")
GEA_res
Neutral_res <- sma(Fst ~ Distance, data=neutral, method = "SMA")
Neutral_res

plot(GEA_res)
plot(Neutral_res)
# transfer into scientific format
a <- c(GEA_res$groupsummary$Slope, GEA_res$groupsummary$Slope_lowCI, GEA_res$groupsummary$Slope_highCI, GEA_res$groupsummary$Int ,GEA_res$groupsummary$Int_lowCI , GEA_res$groupsummary$Int_highCI)
formatC(a, format = "e", digits = 2)
b <- c(Neutral_res$groupsummary$Slope, Neutral_res$groupsummary$Slope_lowCI, Neutral_res$groupsummary$Slope_highCI, Neutral_res$groupsummary$Int ,Neutral_res$groupsummary$Int_lowCI , Neutral_res$groupsummary$Int_highCI)
formatC(b, format = "e", digits = 2)
# test for common slope
sma(Fst ~ Distance*Group, data=IBD)
# test for common intercept
sma(Fst ~ Distance + Group, data=IBD)
# test for neutral slope
sma(Fst ~ Distance, data=GEA_CD5, slope.test=0)
sma(Fst ~ Distance, data=neutral, slope.test=0)

############ GEA outliers 95CI###############
IBD_GEA <- IBD[which(IBD$Group == "GEA_CD5"),]
# fit sma
mod <- sma(Fst ~ Distance, data=IBD_GEA, method="SMA")
# plot model
plot(mod)
# create new data set of Distance at a high resolution (200 points from min to max)
preds_GEA <- data.frame(expand.grid(Distance = seq(min(IBD_GEA$Distance, na.rm = T), max(IBD_GEA$Distance, na.rm = T), length.out = 1000), stringsAsFactors = FALSE))
# bootstrap data and get predictions
preds_GEA <- IBD_GEA %>%
  # create new bootstrapped data sets
  modelr::bootstrap(n = 1000, id = 'boot_num') %>%
  # fit sma to every bootstrap
  group_by(boot_num) %>%
  mutate(., fit = map(strap, ~ sma(Fst ~ Distance, data=data.frame(.), method="SMA"))) %>%
  ungroup() %>%
  # extract intercept and slope from each fit
  mutate(., intercept = map_dbl(fit, ~coef(.x)[1]),
         slope = map_dbl(fit, ~coef(.x)[2])) %>%
  select(., -fit) %>%
  # get fitted values for each bootstrapped model
  # uses the preds_GEA dataframe we made earlier
  group_by(boot_num) %>%
  do(data.frame(fitted = .$intercept + .$slope*preds_GEA$Distance, 
                Distance = preds_GEA$Distance)) %>%
  ungroup() %>%
  # calculate the 2.5% and 97.5% quantiles at each Distance value
  group_by(., Distance) %>%
  dplyr::summarise(., conf_low = quantile(fitted, 0.025),
                   conf_high = quantile(fitted, 0.975)) %>%
  ungroup() %>%
  # add fitted value of actual unbootstrapped model
  mutate(., Fst = coef(mod)[1] + coef(mod)[2]*Distance)

############ neutral 95CI ###############
IBD_neutral <- IBD[which(IBD$Group == "Neutral"),]
# fit sma
mod <- sma(Fst ~ Distance, data=IBD_neutral, method="SMA")
# plot model
plot(mod)
# create new data set of Distance at a high resolution (200 points from min to max)
preds_neutral <- data.frame(expand.grid(Distance = seq(min(IBD_neutral$Distance, na.rm = T), max(IBD_neutral$Distance, na.rm = T), length.out = 1000), stringsAsFactors = FALSE))
# bootstrap data and get predictions
preds_neutral <- IBD_neutral %>%
  # create new bootstrapped data sets
  modelr::bootstrap(n = 1000, id = 'boot_num') %>%
  # fit sma to every bootstrap
  group_by(boot_num) %>%
  mutate(., fit = map(strap, ~ sma(Fst ~ Distance, data=data.frame(.), method="SMA"))) %>%
  ungroup() %>%
  # extract intercept and slope from each fit
  mutate(., intercept = map_dbl(fit, ~coef(.x)[1]),
         slope = map_dbl(fit, ~coef(.x)[2])) %>%
  select(., -fit) %>%
  # get fitted values for each bootstrapped model
  # uses the preds_neutral dataframe we made earlier
  group_by(boot_num) %>%
  do(data.frame(fitted = .$intercept + .$slope*preds_neutral$Distance, 
                Distance = preds_neutral$Distance)) %>%
  ungroup() %>%
  # calculate the 2.5% and 97.5% quantiles at each Distance value
  group_by(., Distance) %>%
  dplyr::summarise(., conf_low = quantile(fitted, 0.025),
                   conf_high = quantile(fitted, 0.975)) %>%
  ungroup() %>%
  # add fitted value of actual unbootstrapped model
  mutate(., Fst = coef(mod)[1] + coef(mod)[2]*Distance)

preds_neutral$Group <- "Neutral"
preds_GEA$Group <- "GEA_CD5"

preds <- rbind(preds_neutral, preds_GEA)
#preds <-preds_neutral
# plot with ggplot

cbPalette <- c("lightcoral", "lightskyblue", "darkseagreen", "black")

ggplot(IBD, aes(Distance, Fst, color=Group, shape=Group, fill=Group)) +
  geom_point(size=2, alpha = 0.6)+
  geom_line(data = preds) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.1, preds, colour = NA) +
  scale_colour_manual(values=cbPalette, breaks=c("GEA_CD5", "Neutral"))+
  scale_shape_manual(values=c(15, 16, 17, 18), breaks=c("GEA_CD5", "Neutral"))+
  scale_y_continuous(name=expression(italic(F)[ST]~"/"~textstyle(group("(", 1-italic(F)[ST], ")")))) +
  #scale_x_continuous(name="Distance(km)", limits=c(0, 30))+
  #cowplot::theme_cowplot()+
  theme(legend.position="right",
        legend.title = element_text(size = 20),
        legend.text=element_text(size=14)) + 
  theme(text=element_text(family="Times New Roman",  size=12, colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = "Black"))+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=rel(1.5)))+
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.5, linetype = "dashed"))

graph2ppt(file="IBD_GEA_neutral3", width=6, height=6)

########################### old SMA plot ################################
# formal plot
mod_SGS_wild = lmodel2(Fst ~ Distance, data=SGS_wild, "interval", "interval", 95)
mod_neutral = lmodel2(Fst ~ Distance, data=neutral, "interval", "interval", 95)

reg_neutral = mod_neutral$regression.results[which(mod_neutral$regression.results$Method == "SMA"),]
reg_SGS_wild = mod_SGS_wild$regression.results[which(mod_SGS_wild$regression.results$Method == "SMA"),]

names(reg_neutral) = c("method", "intercept", "slope", "angle", "p-value")
names(reg_SGS_wild) = c("method", "intercept", "slope", "angle", "p-value")

cbPalette <- c("lightcoral", "lightskyblue", "darkseagreen", "black")
#my.formula <- y ~ x

ggplot(data = IBD, aes(x = Distance, y = Fst, color=Group, shape=Group)) +
  geom_point(size=3, alpha = 0.6)+
  scale_colour_manual(values=cbPalette, breaks=c("SGS_wild", "Neutral"))+
  scale_shape_manual(values=c(15, 16, 17, 18), breaks=c("SGS_wild", "Neutral"))+
  scale_y_continuous(name=expression(italic(F)[ST]~"/"~textstyle(group("(", 1-italic(F)[ST], ")")))) +
  scale_x_continuous(name="Distance(km)", limits=c(0, 32))+
  cowplot::theme_cowplot()+
  theme(legend.position="right",
        legend.title = element_text(size = 20),
        legend.text=element_text(size=14)) + 
  theme(text=element_text(family="Times New Roman", face="bold", size=12, colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = "Black"))+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=rel(1.5)))+
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.5, linetype = "dashed"))+
  geom_abline(data = reg_SGS_wild, aes(intercept = intercept, slope = slope), colour = "lightcoral")+
  geom_abline(data = reg_neutral, aes(intercept = intercept, slope = slope), colour = "lightskyblue")
#stat_poly_eq(formula = my.formula, 
#             aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~~")), 
#             parse = TRUE)  # need to replace the lm intercept and slope with RMA ones. Change them in ppt.

graph2ppt(file="IBD_SGS_neutral", width=12, height=8)

########################
## scatter plot for fst
########################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/06_pairwise_fst_IBD")
df1 = read.delim("fst_comp_SGS.txt", header = T, sep='\t')

quartz(width=6, height=6)
par(mfrow=c(2,2))

ggplot(df1, aes(x=Cohort_2018, y=Cohort_2021)) +
  geom_point(size=2, shape=23)+
  theme(legend.position="right",
        legend.title = element_text(size = 12),
        legend.text=element_text(size=14)) + 
  theme(text=element_text(family="Times New Roman", face="bold", size=12, colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = "Black"))+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=rel(1.5)))+
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.5, linetype = "dashed"))+
  coord_cartesian(xlim = c(0,6e-2), ylim = c(0,6e-2))+
  geom_abline(intercept = 0, slope = 1, colour = "lightcoral")


