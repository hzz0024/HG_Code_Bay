library(ggplot2)
library(export)
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/03_realigned")
df <- read.table("Raw_data_summary.txt",header=TRUE,sep="\t")

p1 <- ggplot(df, aes(x=Group, y=Raw_bases, color=Group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(alpha = 0.1, shape=16, position=position_jitter(0.2)) +
  theme_classic()+
  labs(x=NULL, y = "Raw bases") + 
  theme(axis.text.x = element_text(angle=25,hjust=1)) +
  theme(text = element_text(size=8),
        legend.position = "none")
p1

p2 <- ggplot(df, aes(x=Group, y=Dedup_mapped_bases, color=Group)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(alpha = 0.1, shape=16, position=position_jitter(0.2)) +
  theme_classic()+
  labs(x=NULL, y = "Bases after mapping and filtering") + 
  theme(axis.text.x = element_text(angle=25,hjust=1)) +
  theme(text = element_text(size=8),
        legend.position = "none")
p2


graph2ppt(file="raw_read_diff",width=10,height=6)

###########################
### 2 sd above the mean ###
###########################
library(gridExtra)
library(dplyr)
meta_data <- df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_read = mean(Raw_bases),
    sd = sd(Raw_bases),
    up = mean(Raw_bases) + 2* sd(Raw_bases),
    down = mean(Raw_bases) - 2* sd(Raw_bases),
    up_cnt = sum(Raw_bases > (mean(Raw_bases) + 2* sd(Raw_bases))),
    down_cnt = sum(Raw_bases < (mean(Raw_bases) - 2* sd(Raw_bases))),
    mean_dup = mean(Dedup_mapped_bases),
    sd_dup = sd(Dedup_mapped_bases),
    up_dup = mean(Dedup_mapped_bases) + 2* sd(Dedup_mapped_bases),
    down_dup = mean(Dedup_mapped_bases) - 2* sd(Dedup_mapped_bases),
    up_dup_cnt = sum(Dedup_mapped_bases > (mean(Dedup_mapped_bases) + 2* sd(Dedup_mapped_bases))),
    down_dup_cnt = sum(Dedup_mapped_bases < (mean(Dedup_mapped_bases) - 2* sd(Dedup_mapped_bases))),
  )

t <- tableGrob(meta_data)    
grid.arrange(p1, p2, t)
graph2ppt(file="summary",width=14,height=10)

#############################################
### extract samples based on mean and 2sd ###
#############################################

extract <- function(group){
  db <- df[which(df$Group == group),]
  cut1 <- mean(db$Raw_bases) - 2*sd(db$Raw_bases)
  sample1 <- db$Bam_file[which(db$Raw_bases < cut1)]
  
  cut2 <- mean(db$Dedup_mapped_bases) - 2*sd(db$Dedup_mapped_bases)
  sample2 <- db$Bam_file[which(db$Dedup_mapped_bases < cut2)]
  
  sample = unique(c(sample1, sample2))
  write.table(sample, file = paste0(group, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  print(sample)
} 

extract("DB18_adult_wild")
extract("DB19_adult_wild_challenge")
extract("DB20_adult_challenge1")
extract("DB20_adult_challenge2")
extract("DB21_adult_wild")
extract("DB21_spat_wild")

###############################################
### Split bam files into Rows of Fixed Size ###
###############################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/04_summary")
bam <- read.table("bam_list_final.txt",header=FALSE,sep="\t")
bam_split <- split(bam$V1, ceiling(seq_along(bam$V1)/2))

for(i in seq_along(bam_split)){
  write.table(bam_split[[i]], paste0("./final_bam_list/bam_list_", i,".list"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}





      