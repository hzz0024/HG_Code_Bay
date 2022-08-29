#setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/trash/trial_with_modified_code/finished")

########################################
# compare SGS and permutation outliers #
########################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Challenge_20//")
SGS_outlier = read.delim("20_SGS_Sur_Ref_ps_outlier.list", header = TRUE, sep='\t') 
SGS_outlier_id <- SGS_outlier$id

df = read.delim("out_selection_index.txt", header = TRUE, sep='\t')
df$SNP <- paste0(df$chromo, "_", df$position)
hist(df$pvalue)
outlier <- df[which(df$pvalue < 0.005),]

head(outlier)

# Check shared SNPs
x <- list(
  A = outlier$SNP, 
  B = SGS_outlier_id
)

library("ggvenn")
names(x) <- c("Wild_21_permutation_outliers", "Wild_21_SGS_outliers")

ggvenn(
  x, columns = c("Wild_21_permutation_outliers", "Wild_21_SGS_outliers"),
  fill_color = c("#EFC000FF", "#bc5090"),
  stroke_size = 0.5, set_name_size = 4
)                 

length(intersect(outlier$SNP, SGS_Sur_Ref_20))    

# ######################
# # single snp example #
# ######################
# setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/single_snp_test_example/")
# example1_null = read.delim("example1.txt", header = FALSE, sep='\t')
# example2_null = read.delim("example2.txt", header = FALSE, sep='\t')
# e1 <- ggplot(example1_null, aes(x=V1)) + 
#   geom_histogram(aes(y=..density..), colour="black", fill="white")+
#   geom_density(alpha=.2, fill="#FF6666") +
#   labs(title = "Example SNP 1")+
#   theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
#   xlab("Delta_p") +
#   ylab("Density") + 
#   theme_classic()+
#   theme(legend.title = element_text(size = 12),
#         legend.text=element_text(size=12)) + 
#   theme(text=element_text(family="Times New Roman", size=12, colour="black"))+
#   geom_vline(xintercept=(0.016433-0.073830), linetype="dashed", color = "red", size=2)
# 
# e2 <- ggplot(example2_null, aes(x=V1)) + 
#   geom_histogram(aes(y=..density..), colour="black", fill="white")+
#   geom_density(alpha=.2, fill="#FF6666") +
#   labs(title = "Example SNP 2")+
#   theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
#   xlab("Delta_p") +
#   ylab(NULL) + 
#   theme_classic()+
#   theme(legend.title = element_text(size = 12),
#         legend.text=element_text(size=12)) + 
#   theme(text=element_text(family="Times New Roman", size=12, colour="black"))+
#   geom_vline(xintercept=(0.189444-0.097018 ), linetype="dashed", color = "red", size=2)
# 
# obs1 <- 0.016433-0.073830
# s_index_1 <- (sum(example1_null > obs1)+1)/101
# s_index_2 <- (sum(example1_null < obs1)+1)/101
# s_index = max(s_index_1, s_index_2)
# 
# obs2 <- 0.189444-0.097018 
# s_index_1 <- (sum(example2_null > obs2)+1)/101
# s_index_2 <- (sum(example2_null < obs2)+1)/101
# s_index = max(s_index_1, s_index_2)
# 
# e1 + e2 + theme(legend.position="right") +
#   theme(legend.title = element_text(size = 12),
#         legend.text=element_text(size=12)) 

############################################
# outlier detections among replicates CH20 #
############################################
# a reference allele at each locus increased (Pr[∆pi > ∆pdrifti]) or decreased (Pr[∆pi < ∆pdrifti]) in frequency more than expected under null model
library(export)
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Challenge_20/")
df = read.delim("out_selection_index.txt", header = TRUE, sep='\t')
df$SNP <- paste0(df$chromo, "_", df$position)
hist(df$pvalue)
outlier <- df[which(df$pvalue < 0.01),]
# 1-3 is generated from 100 permutation
df1 = read.delim("out_selection_index_1.txt", header = TRUE, sep='\t')
df2 = read.delim("out_selection_index_2.txt", header = TRUE, sep='\t')
df3 = read.delim("out_selection_index_3.txt", header = TRUE, sep='\t')
cor(df1$pvalue, df2$pvalue)

head(df1)
hist(df1$pvalue)
min(df1$pvalue)
df1$SNP <- paste0(df1$chromo, "_", df1$position)
df2$SNP <- paste0(df2$chromo, "_", df2$position)
df3$SNP <- paste0(df3$chromo, "_", df3$position)
outlier1 <- df1[which(df1$pvalue < 0.01),]
outlier2 <- df2[which(df2$pvalue < 0.01),]
outlier3 <- df3[which(df3$pvalue < 0.01),]
dim(outlier1)
dim(outlier2)
dim(outlier3)


df1$adj = p.adjust(df1$pvalue, method = 'BH')
df2$adj = p.adjust(df2$pvalue, method = 'BH')
df3$adj = p.adjust(df3$pvalue, method = 'BH')

length(intersect(outlier1$SNP, outlier2$SNP))
length(intersect(outlier2$SNP, outlier3$SNP))
length(intersect(intersect(outlier1$SNP, outlier2$SNP), outlier3$SNP))

x <- list(
  A = outlier1$SNP, 
  B = outlier2$SNP
)

library("ggvenn")
names(x) <- c("CH20 100 permutation trial1","CH20 100 permutation trial2")
fill_color = c("#0073C2FF", "#EFC000FF")

ggvenn(
  x, columns = c("CH20 100 permutation trial1", "CH20 100 permutation trial2"),
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

# 4-5 is generated from 200 permutation
df4 = read.delim("out_selection_index_4.txt", header = TRUE, sep='\t')
df5 = read.delim("out_selection_index_5.txt", header = TRUE, sep='\t')
df4$SNP <- paste0(df4$chromo, "_", df4$position)
df5$SNP <- paste0(df5$chromo, "_", df5$position)
cor(df4$pvalue, df5$pvalue)

outlier4 <- df4[which(df4$pvalue < 0.01),]
outlier5 <- df5[which(df5$pvalue < 0.01),]
dim(outlier4)
dim(outlier5)
length(intersect(outlier4$SNP, outlier5$SNP))

x <- list(
  A = outlier4$SNP, 
  B = outlier5$SNP
)

library("ggvenn")
names(x) <- c("CH20 200 permutation trial1","CH20 200 permutation trial2")
fill_color = c("#0073C2FF", "#EFC000FF")

ggvenn(
  x, columns = c("CH20 200 permutation trial1", "CH20 200 permutation trial2"),
  fill_color = c("#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

# share with SGS outliers (FDR < 0.1)
SGS_outlier <- read.delim("20_SGS_Sur_Ref_FDR_outlier.list",  header = TRUE, sep='\t')

share1 <- intersect(outlier1$SNP,  SGS_outlier$id)
share2 <- intersect(outlier2$SNP,  SGS_outlier$id)

share4 <- intersect(outlier4$SNP,  SGS_outlier$id)
share5 <- intersect(outlier5$SNP,  SGS_outlier$id)

x <- list(
  A = share1,
  B = share2,
  C = share4,
  D = share5
)

library("ggvenn")
names(x) <- c("CH20 common SNP SGS & permutation (100 trial 1)","CH20 common SNP SGS & permutation (100 trial 2)", "CH20 common SNP SGS & permutation (200 trial 1)","CH20 common SNP SGS & permutation (200 trial 2)")
fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")

ggvenn(
  x, columns = c("CH20 common SNP SGS & permutation (100 trial 1)","CH20 common SNP SGS & permutation (100 trial 2)", "CH20 common SNP SGS & permutation (200 trial 1)","CH20 common SNP SGS & permutation (200 trial 2)"),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

ggvenn(
  x, columns = c("CH20 common SNP SGS & permutation (100 trial 1)","CH20 common SNP SGS & permutation (100 trial 2)"),
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)
graph2ppt(file="share1",width=10,height=6)

ggvenn(
  x, columns = c("CH20 common SNP SGS & permutation (200 trial 1)","CH20 common SNP SGS & permutation (200 trial 2)"),
  fill_color = c("#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
graph2ppt(file="share2",width=10,height=6)

write.table(target,"target.bed", row.names = FALSE,col.names = FALSE, sep="\t", quote = FALSE)

GAE_outlier <- read.delim("CD5_shared_lfmm_rda.bed",  header = FALSE, sep='\t')
GAE_outlier$SNP <- paste0(GAE_outlier$V1, "_", GAE_outlier$V2)

####################################
# outlier output for each contrast #
####################################

library(export)

output_outlier <- function(name){
  # SGS outliers (ps < 1e-4)
  #name = "20_SGS_Sur_Ref"
  SGS_outlier <- read.delim(paste0(name, "_ps_outlier.list"),  header = TRUE, sep='\t')
  print(paste0("SGS outlier count is: ", dim(SGS_outlier)[1]))
  df = read.delim("out_selection_index.txt", header = TRUE, sep='\t')
  df$SNP <- paste0(df$chromo, "_", df$position)
  hist(df$pvalue)
  outlier <- df[which(df$pvalue < 0.005),]
  print(paste0("Permutation outlier count is: ", dim(outlier)[1]))
  share_id <- intersect(outlier$SNP,  SGS_outlier$id)
  print(paste0("Shared outlier count is: ", length(share_id)))
  bad_regex <- paste("SGS_")
  new_name  = trimws( sub(bad_regex, "", name) )
  write.table(outlier,paste0("./", new_name, "_permutation_outlier.txt"), row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
  
  SGS_shared <- SGS_outlier[which(SGS_outlier$id %in% share_id),]
  SGS_shared = SGS_shared[with(SGS_shared, order(chr, pos)),]
  
  pemt_shared <- df[which(df$SNP %in% share_id),]
  pemt_shared = pemt_shared[with(pemt_shared, order(chromo, position)),]
  
  shared_bed <- data.frame(SGS_shared$chr, SGS_shared$pos, SGS_shared$pos)
  write.table(shared_bed,paste0("./Shared/", new_name, "_shared_outlier.bed"), row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
  
  shared_bed_10K <- data.frame(SGS_shared$chr, SGS_shared$pos-5000, SGS_shared$pos+5000)
  write.table(shared_bed_10K,paste0("./Shared/", new_name, "_shared_10K.bed"), row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
  
  write.table(SGS_shared,paste0("./Shared/", new_name, "_shared_outlier.txt"), row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
}
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Challenge_19/")
output_outlier("19_SGS_Sur_Ref")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Challenge_20/")
output_outlier("20_SGS_Sur_Ref")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Wild_18/")
output_outlier("18_SGS_HC_NB")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Wild_19/")
output_outlier("19_SGS_HC_NB")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Wild_21/")
output_outlier("21_SGS_HC_NB")


#####################
# combine p-values ##
#####################
# 
# library("optparse")
# 
# # function to combine the p-value, multi options are available (this is the new method with correct df estimate)
# combinePValues <- function(..., 
#                            method=c("fisher", "z", "simes", "berger", "holm-middle"), 
#                            weights=NULL, log.p=FALSE, min.prop=0.5)
# {
#   input <- list(...)
#   if (length(input)==1L) {
#     return(unname(input[[1]])) # returning directly.
#   }
#   Np <- unique(lengths(input))
#   if (length(Np) != 1) {
#     stop("all p-value vectors must have the same length")
#   }
#   
#   method <- match.arg(method)
#   switch(method,
#          fisher={
#            if (log.p) {
#              all.logp <- input
#            } else {
#              all.logp <- lapply(input, FUN=log)
#            }
#            
#            n <- integer(Np)
#            X <- numeric(Np)
#            for (i in seq_along(all.logp)) {
#              current <- all.logp[[i]]
#              keep <- !is.na(current)
#              X[keep] <- X[keep] + current[keep]
#              n <- n + keep
#            }
#            
#            n[n==0] <- NA_real_ # ensure that we get NA outputs.
#            pchisq(-2*X, df=2*n, lower.tail=FALSE, log.p=log.p)
#          },
#          simes=combine_simes(input, log.p),
#          `holm-middle`=combine_holm_middle(input, log.p, min.prop),
#          z={
#            if (is.null(weights)) {
#              weights <- rep(1, length(input))
#            } else if (length(weights)!=length(input)) {
#              stop("'length(weights)' must be equal to number of vectors in '...'")
#            } else {
#              check <- unlist(lapply(weights, range))
#              if (any(is.na(check) | check<=0)) {
#                stop("weights must be positive")
#              }
#            }
#            
#            Z <- W2 <- numeric(Np)
#            for (i in seq_along(input)) {
#              current <- input[[i]]
#              keep <- !is.na(current)
#              Z[keep] <- Z[keep] + qnorm(current[keep], log.p=log.p) * weights[[i]]
#              W2 <- W2 + ifelse(keep, weights[[i]]^2, 0)
#            }
#            
#            # Combining p-values of 0 and 1 will yield zscores of -Inf + Inf => NaN.
#            # Here, we set them to 0 to get p-values of 0.5, as the Z-method doesn't
#            # give coherent answers when you have p-values at contradicting extremes.
#            Z[is.nan(Z)] <- 0
#            
#            W2[W2==0] <- NA_real_ # ensure we get NA outputs.
#            pnorm(Z/sqrt(W2), log.p=log.p) 
#          },
#          berger={
#            do.call(pmax, c(input, list(na.rm=TRUE)))
#          }
#   )
# }
# 
# #FUNCTION get.combined():
# #runs fisher.method() on a dataframe of p values
# #ARGUMENTS: DF = a two column data frame with p values in it
# get.combined = function(DF){
#   combined.p.values = c()
#   alt = c()
#   for (i in seq(1, length(DF[,1]))){ #(apparently did not know of nrow() yet)
#     x = DF[i,1]
#     y = DF[i,2]
#     #combined.p = fisher.method(c(x,y))
#     #combined.p = combinePValues(x,y, method='fisher')
#     combined.p = combinePValues(x,y, method='z')
#     combined.p.values = append(combined.p.values, combined.p)
#   }
#   return(combined.p.values)
# }
# 
# # Functio to combine two p-values @HG
# # ARGS: p1_name/p2_name - files with exact test p-values
# combine_p <- function(p1_name, p2_name){
#   dat_1 <- read.delim(p1_name, header = FALSE, sep='\t')
#   dat_2 <- read.delim(p2_name, header = FALSE, sep='\t')
#   p1 <- dat_1$V7
#   p2 <- dat_2$V7
#   cat('file1 length: ')
#   cat(length(p1))
#   cat('file2 length: ')
#   cat(length(p2))
#   cat('\n')
#   ps <- data.frame(p1, p2)
#   return(ps)
# }
# 
# SGS_ps <- read.delim("ps_19_SGS_Sur_Ref.txt",  header = FALSE, sep='\t')
# 
# ps = data.frame(df$pvalue, SGS_ps$V6)
# test_comb_ps <- get.combined(ps)
#######################
# outlier detections ##
#######################
# a reference allele at each locus increased (Pr[∆pi > ∆pdrifti]) or decreased (Pr[∆pi < ∆pdrifti]) in frequency more than expected under null model
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/Challenge_20/")
df = read.delim("out_selection_index.txt", header = TRUE, sep='\t')
head(df)
hist(df$pvalue)
min(df$pvalue)
df$SNP <- paste0(df$chromo, "_", df$position)
df$adj = p.adjust(df$pvalue, method = 'BH')

outlier <- df[which(df$pvalue < 0.01),]

outlier <- df[which(df$pvalue > 0.995),] # this is equivalent to a two-tailed probability of 99% or more that the allele frequency change was greater than expected under null model 
dim(outlier)
hist(outlier$deltap)
write.table(outlier,"permutation_outliers_n_45968.txt", row.names = FALSE,col.names = TRUE, sep="\t", quote = FALSE)

jpeg("delta_p_selection.jpg", width = 8, height = 6, units = 'in', res = 300)
p1 <- ggplot(df, aes(x=deltap, y=pvalue)) + 
  geom_point(alpha = 0.1)+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  xlab("Delta_p") +
  ylab("Selection Index") + 
  labs(caption = "S_index_i = max(Pr[∆pi > ∆p_drift_i], Pr[∆pi < ∆p_drift_i])")+
  theme_classic()

outlier <- df[which(df$pvalue > 0.995),]
dim(outlier)
hist(outlier$deltap)

p2<-ggplot(outlier1, aes(x=deltap, y=pvalue)) + 
  geom_point(alpha = 0.1)+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  xlab("Delta_p") +
  ylab(NULL)+
  labs(caption = "SNPs with selection index of 99.5% or greater (this is equivalent 
       to a two-tailed probability of 99% or more)")+
  theme_classic()

fig <-  p1 + p2  +
  theme(legend.title = element_text(size = 12),
        legend.text=element_text(size=12)) + 
  theme(text=element_text(family="Times New Roman", size=12, colour="black"))
fig
dev.off()

########################################
# compare SGS and permutation outliers #
########################################

SGS_outlier = read.delim("20_SGS_Sur_Ref_FDR_outlier.list", header = TRUE, sep='\t') 
SGS_Sur_Ref_20 <- SGS_outlier$id

head(outlier)

# Check shared SNPs
x <- list(
  A = outlier$SNP, 
  B = SGS_Sur_Ref_20
)

library("ggvenn")
names(x) <- c("2020_permutation_outliers", "2020_SGS_outliers")

ggvenn(
  x, columns = c("2020_permutation_outliers", "2020_SGS_outliers"),
  fill_color = c("#EFC000FF", "#bc5090"),
  stroke_size = 0.5, set_name_size = 4
)                 

length(intersect(outlier$SNP, SGS_Sur_Ref_20))                

#####################################
# compare the deltap and starting p #
#####################################
library(ggplot2)
library(hrbrthemes)
library(ggthemr)
library(dplyr)
library(stringr)
source("manhattan.R")
ggthemr("fresh")

outlier = read.delim("permutation_outliers_n_45968.txt", header = TRUE, sep='\t') 
permute_outliers <- data.frame(outlier$chromo, outlier$position, outlier$Sur20, outlier$Ref20, outlier$deltap, outlier$SNP)
colnames(permute_outliers) <- c("chromo","position","Sur20", "Ref20", "deltap", "SNP")
permute_outliers$type = "Permutation"

SGS_outlier = read.delim("20_SGS_Sur_Ref_FDR_outlier.list", header = TRUE, sep='\t') 
SGS_outliers <- data.frame(SGS_outlier$chr, SGS_outlier$pos, SGS_outlier$p1, SGS_outlier$p2, SGS_outlier$deltap, SGS_outlier$id)
colnames(SGS_outliers) <- c("chromo","position","Sur20", "Ref20", "deltap", "SNP")
SGS_outliers$type = "SGS"

shared_outliers <- permute_outliers[permute_outliers$SNP %in% SGS_outliers$SNP,]
permute_unique <- subset(permute_outliers, !(permute_outliers$SNP %in% shared_outliers$SNP))
SGS_unique <- subset(SGS_outliers, !(SGS_outliers$SNP %in% shared_outliers$SNP))

all_outliers <- rbind(shared_outliers, permute_unique, SGS_unique)
all_outliers %>% 
  mutate(chromo = str_replace(chromo, "NC_035780.1", "1")) %>% 
  mutate(chromo = str_replace(chromo, "NC_035781.1", "2")) %>% 
  mutate(chromo = str_replace(chromo, "NC_035782.1", "3")) %>% 
  mutate(chromo = str_replace(chromo, "NC_035783.1", "4")) %>% 
  mutate(chromo = str_replace(chromo, "NC_035784.1", "5")) %>% 
  mutate(chromo = str_replace(chromo, "NC_035785.1", "6")) %>% 
  mutate(chromo = str_replace(chromo, "NC_035786.1", "7")) %>% 
  mutate(chromo = str_replace(chromo, "NC_035787.1", "8")) %>%
  mutate(chromo = str_replace(chromo, "NC_035788.1", "9")) %>% 
  mutate(chromo = str_replace(chromo, "NC_035789.1", "10"))  -> all_outliers
all_outliers$chromo = as.numeric(all_outliers$chromo)
par(mar=c(5,5,4,2))
manhattan(chr="chromo",bp="position",p="deltap", snp = "SNP", all_outliers, highlight1 = shared_outliers$SNP,  logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(Delta~p), cex.lab=1.5) 

manhattan(chr="chromo",bp="position",p="deltap", snp = "SNP", all_outliers, highlight2 = SGS_outliers$SNP,  logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(Delta~p), cex.lab=1.5) 

plot(outlier$deltap)
plot(SGS_outlier$deltap)
