# Patterns of ROH analysis
# This script contains the analyses for the first part of the paper,
# exploring ROH variation across the genome, ROH islands and deserts,
# and the association between ROH prevalence and recombination rate variation.

library(data.table)
library(tidyverse)
source("theme_simple.R")
library(windowscanr)
library(cowplot)
library(gt)
library(grid)
library(ggplotify)
library(patchwork)
library(viridis)
library(scales)
library(ggridges)
library(ggpubfigs)
library(gt)
library(lme4)
library(ggpubfigs)
library(broom)
library(dplyr)
library(export)
library(zoo)
options(scipen=9999)

#################
# plot for FROH #
#################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/Formal_plot")
library(cowplot)
library(reshape2)

df <- read.delim("Summary_ROH_per_breed_oyster_k200_individual_FROH.csv", header = TRUE, sep=',')
df$Source <- factor(df$Source, levels=c("Wild", "Selected"))
#   MEW1        MEW2      LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
col_gradient <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                   #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
                   "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#e1bb94", "#fbd0a5", "#b58383")
order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
df$FID <-factor(df$FID, levels=order1)
df$FROH1_=df$FROH1_/100

FROH <- df %>%
  ggplot(aes(FID, FROH1_, fill=FID)) +
  #geom_violin(trim = FALSE, alpha=0.4, show.legend = FALSE, width=6) +
  geom_boxplot(width = 0.4, alpha=0.8, aes(fill=FID), outlier.shape = NA)  +
  geom_jitter(alpha=0.2, width = 0.1)+
  scale_fill_manual(values = col_gradient) +
  scale_color_manual(values = col_gradient) +
  #scale_fill_manual(values = c("#049DD9", "#F25C05"))+
  xlab(NULL) + ylab(expression("Genomic"~"Inbreeding"~"Coefficient"~"("~italic(F)[ROH]~")"))+
  theme_classic()+
  theme(#=strip.placement = 'outside',
    panel.spacing = unit(0, "lines"),
    legend.position="none",
    axis.title.y = element_text(margin=margin(r=-50))) + # control the label position y is r 
  theme(text=element_text(family="Times New Roman", size=12, colour="black"),
        axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))
 
FROH <- FROH + ggtitle("(a)")
FROH

# DT = read.delim("ROH_analyse_population_x.hom.txt", header = TRUE, sep='\t')
# plotdat <- data.frame(ROH=DT$NSEG, Population=DT$FID)
# 
# wilcox_dt = read.delim("ROH_analyse_population_x.hom.wildsel.txt", header = TRUE, sep='\t')
# wilcox.test(wilcox_dt$NSEG[wilcox_dt$FID == "WILD"],wilcox_dt$NSEG[wilcox_dt$FID == "SEL"])


#################
# plot for SROH #
#################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/Formal_plot")
library(cowplot)
library(scales)
library(dplyr)
library(magrittr)

df <- read.delim("Summary_ROH_per_breed_oyster_k200_individual_SROH.csv", header = TRUE, sep=',')
df %>%
  group_by(FID) %>%
  summarise_at(vars(F_ROHall, SROH02_05, SROH05_1, SROH1_2, SROH2_4, SROH4_8, SROH8_16, SROH16_, SROH1_, NSEG, NSEG02_05, NSEG05_1, NSEG1_2, NSEG2_4, NSEG4_8, NSEG8_16, NSEG16_, NSEG1_, KB, KBAVG), sum, ra.rm=FALSE) ->  df_summary
df_summary %<>% mutate_at(10:18, as.integer) # change the NSEG to integers
df_summary <- as.data.frame(df_summary)

write_delim(df_summary, path = "./SROH_summary.txt")
SROH_df <- read.delim("SROH_summary_by_pop.txt", header = TRUE, sep='\t')
wilcox.test(SROH_df$SROH1_[SROH_df$Source == "Selected"], SROH_df$SROH1_[SROH_df$Source == "Wild"])
wilcox.test(SROH_df$NSEG1_[SROH_df$Source == "Selected"], SROH_df$NSEG1_[SROH_df$Source == "Wild"])

df_summary %>% 
  group_by(FID) %>%
  summarise(
    FID=FID,
    p02_05 = (SROH02_05/F_ROHall)*100,
    p05_1 = (SROH05_1/F_ROHall)*100,
    p1_2 = (SROH1_2/F_ROHall)*100,
    p2_4 = (SROH2_4/F_ROHall)*100,
    p4_8 = (SROH4_8/F_ROHall)*100,
    p8_16 = (SROH8_16/F_ROHall)*100,
    p16_ = (SROH16_/F_ROHall)*100,
    np02_05 = (NSEG02_05/NSEG)*100,
    np05_1 = (NSEG05_1/NSEG)*100,
    np1_2 = (NSEG1_2/NSEG)*100,
    np2_4 = (NSEG2_4/NSEG)*100,
    np4_8 = (NSEG4_8/NSEG)*100,
    np8_16 = (NSEG8_16/NSEG)*100,
    np16_ = (NSEG16_/NSEG)*100) -> df_percentage
df_percentage
write_delim(df_percentage, path = "./SROH_perc_summary.txt")
# _______________________________________________________________________________
# plot for sum of ROH across individuals in each population SROH
extract_df <- data.frame(df$FID, df$SROH02_05, df$SROH05_1, df$SROH1_2, df$SROH2_4, df$SROH4_8)
extract_df <- data.frame(df_summary$FID, df_summary$SROH02_05, df_summary$SROH05_1, df_summary$SROH1_2, df_summary$SROH2_4, df_summary$SROH4_8)



colnames(extract_df) = c("Pop", "SROH02_05", "SROH05_1", "SROH1_2", "SROH2_4", "SROH4_8")
plot_df <- melt(extract_df, id.vars=c("Pop"),
                variable.name = "Interval", 
                value.name = "Values")
plot_df$Values <- plot_df$Values/1000
order1 = rev(c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2"))
plot_df$Pop <-factor(plot_df$Pop, levels=order1)
#   MEW1        MEW2      LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
cbPalette <- rev(c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
                "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#e1bb94", "#fbd0a5", "#b58383"))

SROH_by_pop <- plot_df %>%
  ggplot(aes(Interval, Values, fill = factor(Pop))) +
  geom_bar(stat = "identity",  width = 0.8, position=position_dodge(width = 0.8), alpha=0.8) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) +
  xlab(NULL) + ylab(expression("Total"~"Sum"~"of"~"ROH"~"(Mb)"))+
  scale_x_discrete(labels=c("SROH02_05" = "0.2<ROH<0.5Mb", 
                            "SROH05_1" = "0.5<ROH<1Mb",
                            "SROH1_2" = "1<ROH<2Mb",
                            "SROH2_4" = "2<ROH<4Mb",
                            "SROH4_8" = "4<ROH<8Mb"))+
  theme_classic()+
  theme(legend.title = element_text(size = 10),
        legend.text=element_text(size=10), 
        legend.position = c("none"))+
  theme(text=element_text(family="Times New Roman", face="bold", size=12, colour="black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  coord_flip()

SROH_by_pop <- SROH_by_pop + ggtitle("(b)")
SROH_by_pop
# _______________________________________________________________________________
# plot for the number of ROH segments NROH
extract_df <- data.frame(df_summary$FID, df_summary$NSEG02_05, df_summary$NSEG05_1, df_summary$NSEG1_2, df_summary$NSEG2_4, df_summary$NSEG4_8)
colnames(extract_df) = c("Pop", "NROH02_05", "NROH05_1", "NROH1_2", "NROH2_4", "NROH4_8")
plot_df <- melt(extract_df, id.vars=c("Pop"),
                variable.name = "Interval", 
                value.name = "Values")
order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
plot_df$Pop <-factor(plot_df$Pop, levels=order1)
#   MEW1        MEW2      LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
cbPalette <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
                "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#e1bb94", "#fbd0a5", "#b58383")

NROH_by_pop <- plot_df %>%
  ggplot(aes(Interval, Values, fill = factor(Pop))) +
  geom_bar(stat = "identity",  width = 0.8, position=position_dodge(width = 0.8), alpha=0.8) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) +
  xlab(NULL) + ylab(expression("Total"~"Sum"~"of"~"ROH"~"(Mb)"))+
  scale_x_discrete(labels=c("NROH02_05" = "0.2<ROH<0.5Mb", 
                            "NROH05_1" = "0.5<ROH<1Mb",
                            "NROH1_2" = "1<ROH<2Mb",
                            "NROH2_4" = "2<ROH<4Mb",
                            "NROH4_8" = "4<ROH<8Mb"))+
  theme_classic()+
  theme(legend.title = element_text(size = 10),
        legend.text=element_text(size=10), 
        legend.position = c("right"))+
  theme(text=element_text(family="Times New Roman", face="bold", size=12, colour="black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  coord_flip()

NROH_by_pop <- NROH_by_pop + guides(fill=guide_legend(title="Population/Line"))
NROH_by_pop
# _______________________________________________________________________________
# Start to plot for selected vs wild groups
df$Source <-factor(df$Source, levels=c("Selected", "Wild"))
extract_df <- data.frame(df$FID, df$SROH02_05, df$SROH05_1, df$SROH1_2, df$SROH2_4, df$SROH4_8, df$SROH8_16, df$Source)
colnames(extract_df) = c("Pop", "SROH02_05", "SROH05_1", "SROH1_2", "SROH2_4", "SROH4_8", "SROH8_16", "Source")
plot_df <- melt(extract_df, id.vars=c("Pop","Source"),
                 variable.name = "Interval", 
                 value.name = "Values")
plot_df$Values <- plot_df$Values/1000
                  # Selected     Wild      
col_gradient <- c( "orange", "skyblue") 
SROH <- plot_df %>%
  ggplot(aes(Interval, Values)) +
  #geom_violin(trim = FALSE, alpha=0.4, show.legend = FALSE, width=2) +
  geom_boxplot(width = 0.4, alpha=0.8, aes(fill=Source), outlier.shape = NA)  +
  #geom_jitter(alpha=0.2, width = 0.1)+
  #scale_fill_manual(values = col_gradient) +
  scale_color_manual(values = col_gradient) +
  xlab(NULL) + ylab(expression("Total"~"Sum"~"of"~"ROH"~"(Mb)"))+
  scale_x_discrete(labels=c("SROH02_05" = "0.2<ROH<0.5Mb", 
                            "SROH05_1" = "0.5<ROH<1Mb",
                            "SROH1_2" = "1<ROH<2Mb",
                            "SROH2_4" = "2<ROH<4Mb",
                            "SROH4_8" = "4<ROH<8Mb",
                            "SROH8_16" = "8<ROH<16Mb"))+
  ylim(0,2500)+
  theme_classic()+
  theme(legend.title = element_text(size = 10),
        legend.text=element_text(size=10), 
        legend.position = c(0.7, 0.9)) +
  theme(text=element_text(family="Times New Roman", face="bold", size=12, colour="black"),
       axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  coord_flip()

SROH <- SROH+ggtitle("(b)") + guides(fill=guide_legend(title="Origin"))
SROH


#########################
# plot for SROH vs NROH #
#########################
library(ggrepel)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/Formal_plot")

ROH_plot <- read_delim("TableS3_SROH_NROH_summary.txt", delim = "\t")

order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
ROH_plot$Pop <-factor(ROH_plot$Pop, levels=order1)
col <- c( "#F62A00", "#325A98")

SROH_NSEG <- ggplot(ROH_plot, aes(x = SROH1_, y = NSEG1_, color = Origin, shape=Origin)) +
  geom_point(size=2.5)+
  scale_color_manual(values = col) +
  geom_label_repel(aes(label = Pop),size = 3.2, fill = "white", box.padding = 0.2, max.overlaps = 20, show.legend = FALSE)+
  xlab("Sum total length of ROH (Mb)") + ylab("Total number of ROH (NROH)") + 
  theme_classic()+
  theme(legend.title =element_text(size = 12),
        legend.text=element_text(size=12), 
        legend.position = c("top"))+
  theme(legend.position=c(.2,.85), 
        text=element_text(family="Times New Roman", size=12, colour="black"))

SROH_NSEG <- SROH_NSEG + ggtitle("(b)")
SROH_NSEG

#############################
# plot for ROH distribution #
#############################
#setwd("~/Dropbox/Mac/Documents/HG/Github/References/ROH/sheep_ID-master")
setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/All_n_509_formal_run_plot/")
# Chr lengths
chr_data <- read_delim("data/chromosome_info_Cv30.txt", delim = "\t") %>% 
        dplyr::rename(size_BP = Length, CHR = Part) %>% 
        mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
        .[2:11, ] %>% 
        summarise(sum_KB = sum(size_KB)) %>% 
        as.numeric()

#~~ ROH for survival data subset
file_path <- "output/ROH/All_ROH_distribution.txt"
roh_lengths <- fread(file_path)

# max ROH length
roh_lengths[which.max(roh_lengths$KB), ]

# ROH overview -----------------------------------------------------------------
# descriptive ROH statistics
num_roh_per_ind <- roh_lengths %>% group_by(IID) %>% tally() 
summary(num_roh_per_ind$n)
sd(num_roh_per_ind$n)

# inbreeding coefficients
froh <- roh_lengths %>%
        dplyr::group_by(IID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/autosomal_genome_size) %>% 
        mutate(FROH_cent = FROH - mean(FROH))
mean(froh$FROH)
range(froh$FROH)

# longest ROH
roh_lengths %>% arrange(-KB)

# longest ROH proportional to chr size
chr_sizes <- chr_data %>% .[-1, ] %>% 
        mutate(CHR = str_replace(CHR, "Chromosome ", "")) %>% 
        mutate(CHR = as.integer(CHR))

roh_lengths %>% 
        arrange(desc(KB)) %>% 
        left_join(chr_sizes, by = "CHR") %>% 
        mutate(prop_chr = KB / size_KB) %>% 
        arrange(desc(prop_chr))
#plot(froh$KBSUM, num_roh_per_ind$n)

# ROH length and abundance in the 1% least and most inbred individuals 
num_roh_per_ind %>% 
        left_join(froh) %>% 
       # top_frac(-0.01, FROH) %>%    # top 1% least inbred individuals
        top_n(-7, FROH) %>%      # top 1% most inbred individuals
        summarise(mean(n), mean(KBAVG))

# Supplementary Figure FROH / ROH across individuals ---------------------------
p_froh <- ggplot(froh, aes(FROH)) +
        geom_histogram(bins = 100,  color = "white", size = 0.1, position = "identity",
                       alpha = 1, fill = "#4c566a") +
        ylab("individuals") +
        xlab(expression(inbreeding~coefficient~F[ROH])) +
        scale_y_continuous(expand = c(0, 0)) 
        #theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12) 
p_froh

p_roh <- ggplot(num_roh_per_ind, aes(n)) +
        #geom_histogram(binwidth = 1,  fill = "#E5E9F0", color = "black",  size = 0.1) +
        geom_histogram(binwidth = 1, color = "white", size = 0.1, position = "identity",
                       alpha = 1, fill = "#4c566a") +
        ylab("individuals") +
        xlab("ROH per genome") +
        scale_y_continuous(expand = c(0, 0)) 
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12, 
                     base_family = "Lato") 
p_roh

p_roh_dist <- p_froh + p_roh + plot_annotation(tag_levels = 'a') &
                theme(plot.tag = element_text(face = "bold"))
p_roh_dist
ggsave("figs/Sup_ROH_dist.jpg", p_roh_dist, width = 7, height = 2.5)

#############################
# plot for ROH heatmap      #
#############################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/All_n_509_formal_run_plot")
all_roh <- roh_lengths %>% 
         dplyr::group_by(IID) %>% 
         dplyr::summarise(sum_roh = sum(KB)) %>% 
         ungroup() %>% 
         arrange(desc(sum_roh))

# edit sample_id_to_modify.txt file to have a vector following your target list
write_delim(all_roh, path = "./sample_id_to_modify.txt")
# after editing, now reload the sample list for next step
sample_list <- read.delim("sample_id_to_modify_sorted.csv", header = TRUE, sep=',')

# for 2 clusters

df <- roh_lengths %>%
        mutate(POS1 = POS1 / 1e+6,
               POS2 = POS2 / 1e+6,
               MB = KB / 1000)

df <- df %>% 
        mutate(IID = factor(IID, levels = sample_list$IID))

yax <- data.frame(IID = fct_inorder(levels(df$IID))) %>%
        mutate(yax = seq(from = 2,
                         to = 2*length(unique(df$IID)),
                         by = 2)) 
num_ind <- dim(yax)[1]

df <- left_join(df, yax, by = "IID")

order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
df$FID <-factor(df$FID, levels=order1)

shade <- df %>%
        dplyr::group_by(CHR) %>%
        dplyr::summarise(min = min(POS1), max = max(POS2)) %>%
        mutate(min = case_when(CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10  ~ 0,
                               TRUE ~ min)) %>%
        mutate(max = case_when(CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10  ~ 0,
                               TRUE ~ max))

#col <- c("#2b2c2e", "#FF0000")
#col <- c("#1E3231", "#9aadbf")
col <- c( "#f67321", "#000000")
chr_names <- as.character(1:10)
names(chr_names) <- as.character(1:10)
#chr_names[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)] <- ""

# df is the freqment plot without labels
df %>% 
        filter(MB > 0) %>% 
        filter(CHR %in% 1:10) %>% 
        ggplot() +
        geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
                  alpha=0.5, fill = "#eceff4") + # "#f7f7f7" "#eceff4"
        geom_hline(data = yax, aes(yintercept = yax), color = "#ffffff", size = 0.4) + 
        #geom_hline(data = yax, aes(yintercept = yax), color = "black", size = 0.4) +
        geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 1.6, ymax = yax + 1.8, 
                      fill = as.factor(CHR)),  col = "grey", size = 0, alpha = 1, inherit.aes=FALSE) + 
        scale_fill_manual(values = rep(col, 10)) + 
        scale_color_manual(values = rep(col, 10)) +
        scale_y_reverse(expand = c(0, 0)) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13, base_family = "Helvetica") +
        facet_grid(~CHR,scales = 'free_x', space = 'free_x', switch = 'x',
                   labeller = as_labeller(chr_names)) +
        theme(#=strip.placement = 'outside',
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                panel.spacing = unit(0, "lines"),
                plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
                axis.line.x = element_blank(),
                legend.position="right",
                axis.title.x = element_text(margin=margin(t=-25)), # control the label position x is t
                axis.title.y = element_text(margin=margin(r=-30)), # control the label position y is r
                axis.text.y = element_text(colour = "white"),
                axis.line.y = element_blank()) +
        coord_cartesian(clip = 'off') +
        xlab("Chromosome") +
        ylab("Individuals") -> ROH_per_ind
ROH_per_ind
# start to add the popultion labels
# first is the number of samples per populations
FID_cnt <- sample_list %>%
  dplyr::group_by(FID) %>% 
  dplyr::summarise(n=n())
target <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
FID_cnt <- FID_cnt[match(target, FID_cnt$FID),]
FID_cnt

gp.label <- data.frame(FID_cnt$FID, FID_cnt$n, c(rep("Wild", 9), rep("Selected", 8)))
colnames(gp.label) <- c("POP", "N", "Type")
gp.label$Y <- rollmean(c(0, cumsum(gp.label$N*2)),2)
gp.label$X <- x_axis <- min(ggplot_build(ROH_per_ind)$data[[1]]$xmin)-10
gp.label$POP <- factor(gp.label$POP, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2"))
gp.label$CHR <- 1 # this step is important to ensure the points are added to CHR 1 panel!
                #   MEW1        MEW2      LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
cbPalette <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                   #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
                   "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#e1bb94", "#fbd0a5", "#b58383")
gp.line <- data.frame(seq(1:length(c(0, cumsum(gp.label$N*2)))) , c(0, cumsum(gp.label$N*2)))
colnames(gp.line) <- c("Number", "Pos")

DROH <-  ROH_per_ind +
  geom_point(data = gp.label,aes(x=X,y=Y,color=POP,shape=Type),size=4, alpha=0.95) +
  geom_hline(data = gp.line, aes(yintercept = Pos), color = "grey", size=0.8, alpha=0.5) + 
  scale_x_continuous(expand = expand_scale(mult = c(0.1, 0)))+
  #scale_x_continuous(limits = c(min(gp.label$X) - ddelta/2 , max(ggplot_build(outlier.gp$genotypes)$data[[1]]$x))) + 
  scale_shape_manual(values=c(15, 17), name="Origin") + 
  scale_color_manual(values=cbPalette, name="Population/Line") +
  theme(legend.position="right") +
  theme(legend.title = element_text(size = 10),
        legend.text=element_text(size=10)) + 
  guides(fill = "none")+ 
  # guides(color="none")+
  # guides(shape = "none") +
  theme(text=element_text(family="Times New Roman", size=12, colour="black"))

DROH <- DROH+ ggtitle("(c)")
DROH
layout <-" 
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
BBBBBCCCCCCCCC
BBBBBCCCCCCCCC
BBBBBCCCCCCCCC
BBBBBCCCCCCCCC
"
final_ROH <- FROH + SROH_NSEG + DROH + plot_layout(design=layout, guides="keep")
final_ROH
ggsave("figs/final_ROH.jpg", final_ROH, width = 12, height = 10, dpi = 300)
graph2ppt(file="figs/final_ROH", width=12, height=10) 

###### end of ROH plot ######


































ggsave("figs/roh_per_ind.jpg", fig, width = 7, height = 3, dpi = 300)
graph2ppt(file="figs/roh_per_ind", width=8, height=4) 

# needs to be a grob to display axis labels correctly
pg <- ggplotGrob(ROH_per_ind)

for(i in which(grepl("strip-b", pg$layout$name))){
        pg$grobs[[i]]$layout$clip <- "off"
}

ROH_per_ind_grob <- as.ggplot(pg)
ggsave("figs/Fig1A.pdf", ROH_per_ind_grob, width = 6, height = 2.5)
#ggsave("figs/fig1b_roh_per_ind_5Mb.jpg", ROH_per_ind_grob, width = 6, height = 3.5)

#p1 / p_roh_classes + plot_layout(heights = c(1, 0.8))

# ROH classes ------------------------------------------------------------------

# # expected ROH length using cM/Mb from Johnston et al (2016)
# length_dist <- data.frame(g = c(1, 2,2^2, 2^3, 2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13)) %>%
#         mutate(ROH_length_cM = 100 / (2*g)) %>% 
#         mutate(ROH_length_Mb = ROH_length_cM * 0.7816743)
# 
# prop_IBD_df <- roh_lengths %>%
#         mutate(length_Mb = KB/1000) %>%
#         mutate(class = case_when(#length_Mb >= 39.083715000 ~ 1,
#                                 # length_Mb < 39.083715000 & length_Mb >= 19.541857500 ~ 2,
#                                  length_Mb >= 19.541857500 ~ 2,
#                                  length_Mb < 19.541857500 & length_Mb >= 9.770928750 ~ 4,
#                                  #length_Mb < 9.770928750 & length_Mb >= 6.513952500 ~ 6,
#                                  length_Mb < 9.770928750& length_Mb >= 4.885464375 ~ 8,
#                                  # length_Mb < 4.885464375 & length_Mb >= 3.908371500 ~ 10,
#                                  length_Mb < 4.885464375 & length_Mb >= 2.442732188 ~ 16,
#                                  length_Mb < 2.442732188 & length_Mb >= 1.2 ~ 32)) %>% # 0.610683047 1.221366094
#         mutate(length_class = case_when(
#                 #class == 1 ~ ">39 (1G)",
#                 #class == 2 ~ ">19.5-39 (2g)",
#                 class == 2 ~ ">19.5 (2g)",
#                 class == 4 ~ "9.8-19.5 (2-4g)",
#                 #class == 6 ~ "6.5-9.7 (6G)",
#                 class == 8 ~ "4.9-9.8 (4-8g)",
#                 # class == 10 ~ "3.9-4.9 (10G",
#                 class == 16 ~ "2.4-4.9 (8-16g)",
#                 class == 32 ~ "1.2-2.4 (16-32g)"
#                 # class == 128 ~ "0.6-0.3 (128G)"
#         )) %>% 
#         mutate(length_class = fct_reorder(length_class, class)) %>% 
#         mutate(IID = as.character(IID)) %>% 
#         group_by(IID, class, length_class) %>%
#         dplyr::summarise(prop_IBD = sum(length_Mb / (autosomal_genome_size/1000))) #%>% 
# 
# # add IBD of non-ROH snps if wanted
# #  bind_rows(homs) 
# 
# prop_IBD_df_with_0 <- prop_IBD_df %>% 
#         # add missing length classes as 0
#         ungroup() %>% 
#         tidyr::complete(length_class, nesting(IID)) %>% 
#         mutate(class = ifelse(is.na(class), length_class, class)) %>% 
#         mutate(prop_IBD = ifelse(is.na(prop_IBD), 0, prop_IBD))
# 
# prop_IBD_df_with_0 %>% 
#         group_by(length_class) %>% 
#         summarise(mean(prop_IBD))
# 
# prop_IBD_df_with_0 %>% 
#         group_by(length_class) %>% 
#        #summarise(sum(prop_IBD > 0)/ 5925)
#         filter(prop_IBD > 0) %>% 
#         summarise(mean(prop_IBD))
# 
# library(viridis)
# library(gghalves)
# col_pal <- plasma(6)
# col_pal <- paste0("#", (c("21295c","204683","1763a1","0a96d6","65bee2","BAE2F2")))
# col_pal <- paste0("#", (c("01161e","124559","598392","84a3a1","aec3b0","eff6e0")))
# #col_pal <- paste0("#", c("432371","714674","9f6976","cc8b79","e39d7a","faae7b"))
# p_roh_length <- prop_IBD_df_with_0 %>% 
#         mutate(prop_IBD = prop_IBD * 100) %>% 
#         ggplot(aes(length_class, prop_IBD, fill = length_class)) +
#         geom_half_point(side = "l", shape = 21, alpha = 0.3, stroke = 0.1, size =2, color = "#4c566a",
#                         transformation_params = list(height = 0, width = 1.3, seed = 1)) +
#         geom_half_boxplot(side = "r", outlier.color = NA,
#                           width = 0.6, lwd = 0.3, color = "black",
#                           alpha = 0.8) +
#         theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13,
#                         base_family = "Helvetica") +
#         ylab("% genome") +
#         scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
#         scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
#         theme(legend.position = "none",
#               plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
#               #axis.ticks.x = element_blank(),
#               axis.title=element_text(size = rel(1.1)), 
#               axis.text = element_text(color = "black")) + 
#         xlab("ROH length class in Mb (~ generations to MRCA)") 
# p_roh_length
# 
# ggsave("pics/Fig1B.jpg", p_roh_length, width = 6.3, height = 2.65)
# ggsave("pics/Fig1B.pdf", p_roh_length, width = 6.3, height = 2.65)
#~~~ ROH density ---------------------------------------------------------------
hom_sum <- fread("output/ROH/ROH_analyse_population_x.hom.summary")
cnt_ind <- 262

hom_sum <- hom_sum %>%
        mutate(MB = BP / 1000000,
               KB = BP / 1000,
               index = 1:nrow(.))

running_roh <- winScan(x = hom_sum,
                       groups = "CHR",
                       position = "KB",
                       values = "UNAFF",
                       win_size = 500,
                       win_step = 500,
                       funs = c("mean"))

# ROH sharing 1
# remove windows without snps
running_roh %>% 
        mutate(UNAFF_mean = UNAFF_mean/cnt_ind) %>% 
        filter(UNAFF_n > 0) -> running_roh_p

#cor(running_roh$UNAFF_mean, running_roh$UNAFF_n, use = "complete.obs")
#plot(running_roh$UNAFF_mean, running_roh$UNAFF_n)

# check distribution of snps
UNFAFF_dis <- ggplot(running_roh, aes(win_mid, UNAFF_n)) + 
  geom_point() + geom_smooth(span = 0.1) +
  facet_wrap(~CHR, scales = "free_x")
UNFAFF_dis
ggsave("figs/UNAFF_n_distribution.jpg", UNFAFF_dis, width = 18, height = 12, dpi = 300)

# ROH across the genome plot
fill_cols <- viridis(20, option = "D")
qn <- scales::rescale(quantile(running_roh_p$UNAFF_mean,
                               probs=seq(0, 1, length.out=length(fill_cols))))

p1 <- ggplot(running_roh_p, aes(x = win_start, y = 0.5, fill = UNAFF_mean)) + 
        #geom_tile(color = "grey", size = 0) +
        geom_tile() +
        theme_simple(base_size = 13, grid_lines = FALSE, base_family = "Times") + 
        scale_y_continuous(expand = c(0,0))+
        scale_x_continuous(expand = expansion(mult = c(0, 0.5), # @HG https://stackoverflow.com/questions/36561030/how-expand-ggplot-bar-scale-on-one-side-but-not-the-other-without-manual-limits
                                                 add = c(0, 0)), 
                           breaks = seq(0, 100000, by = 20000),
                           labels = as.character(seq(0, 100, 20)))+
        ylab("Chromosome") +
        scale_fill_gradientn("% of oyster with ROH",
                             colors = rev(fill_cols), 
                             values = qn,
                             breaks = c(0.1, 0.2, 0.3),
                             labels = c(10, 20, 30)) +
        facet_grid(CHR~., switch="both") +
        xlab("Position in Mb") +
        theme(panel.spacing.y=unit(0.1, "lines"),
              axis.title.x = element_text(margin=margin(t=5)),
              axis.title.y = element_text(margin=margin(r=5)),
              axis.text.y = element_blank(),
              axis.text.x = element_text(color = "black"),
              axis.ticks.x = element_line(size = 0.3),
              #axis.title.y = element_blank(),
              plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = c(0.75,0.73),
              legend.direction = "horizontal",
              strip.text.y.left = element_text(size = 10, angle = 0),
              axis.line.x = element_line(size = 0.1)
              ) +
        guides(fill = guide_colourbar(title.position = "bottom" ,
                                      barwidth = 10.95, barheight = 0.5))
p1
ggsave("figs/ROH_dis_main.jpg", p1, width = 18, height = 12, dpi = 300)
graph2ppt(file="figs/ROH_dis_main", width=12, height=10) 

x <- running_roh_p$UNAFF_mean
y <- density(x, n = 2^12)
p2 <- ggplot(data.frame(x = y$x, y = y$y), aes(x, y)) + 
        geom_line() + 
        geom_segment(aes(xend = x, yend = 0, color = x)) +
        scale_x_continuous(expand = c(0, 0)) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_family = "Times") +
        scale_color_gradientn("Proportion of Sheep with ROH",
                             colors = rev(fill_cols), values = qn,
                             breaks = c(0.1, 0.2, 0.3),
                             labels = c(10, 20, 30)) +
        theme(legend.position = "none",
              plot.margin = margin(1, 1, 1, 1, unit = "cm"),
              axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank()) + 
        scale_y_continuous(expand = c(0, 0), breaks = c(1, 3, 5)) +
        ylab("Density")
        
p2
ggsave("pics/Fig1C_legend.jpg", p2, width = 4, height = 2)
# ggplot(running_roh_p, aes(UNAFF_mean, "test", fill = ..x..)) +
#         geom_density_ridges_gradient(scale = 2.5, lwd = 0.1) +
#         scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) +
#        # theme_simple() +
#         scale_fill_gradientn("Proportion of Sheep with ROH",
#                              colors = rev(fill_cols), values = qn,
#                              breaks = c(0.1,0.3, 0.5, 0.7, 0.9)) +
#         theme(legend.position = "none",
#               plot.margin = margin(1, 1, 1, 1, unit = "cm"))

p2
ggsave("figs/roh_genome_legend.jpg", p2, width = 5, height = 2.5)

# try simple combined plot
p_roh_comb_simple <- plot_grid(ROH_per_ind_grob, p1, nrow = 2, 
                               rel_heights = c(0.658, 1), label_size = 15, 
                               labels = c("A", "B"), align = "v")
p_roh_comb_simple
ggsave("figs/roh_patterns_simple.jpg", p_roh_comb_simple, width = 7, height = 6.5)

# combine legend density and plot in keynote or somewhere else for full plot.

# save all
p_roh_length
grid1 <- plot_grid(ROH_per_ind_grob, p_roh_length, nrow = 2, 
          rel_heights = c(1,0.658), label_size = 15, 
          labels = c("A", "B"), align = "v")
grid2 <- plot_grid(grid1, p1, nrow = 1)
grid2
ROH_per_ind_grob

# ROH ISLANDS AND DESERTS ------------------------------------------------------

#~~~ ROH density
hom_sum <- fread("output/ROH/ROH_analyse_population_x.hom.summary") # ROH_surv_subset/
head(hom_sum)

hom_sum <- hom_sum %>%
        mutate(MB = BP / 1000000,
               KB = BP / 1000,
               index = 1:nrow(.))# %>% 
head(hom_sum)
# overview over SNPs with most and least ROH
hom_sum %>% 
        filter(UNAFF > 0) %>% 
        mutate(prop_roh = UNAFF/cnt_ind) %>% 
        arrange(prop_roh) %>% 
        .[1:50, ]

hom_sum %>% 
        filter(UNAFF > 0) %>% 
        mutate(prop_roh = UNAFF/cnt_ind) %>% 
        arrange(desc(prop_roh)) %>% 
        .[1:50, ]

# count ROH in running windows of 500 Kb
# UNAFF_n is number of SNPs in a window
# UNAFF_mean is mean ROH prevalence in a window
running_roh <- winScan(x = hom_sum,
                       groups = "CHR",
                       position = "KB",
                       values = "UNAFF",
                       win_size = 500,
                       win_step = 500,
                       funs = c("mean"),
                       cores = 8)
head(running_roh)
plot(running_roh$UNAFF_mean, running_roh$UNAFF_n)
cor(running_roh$UNAFF_mean, running_roh$UNAFF_n, use = "complete")

# remove NAs
running_roh <- na.omit(running_roh)
max(running_roh$UNAFF_mean)
min(running_roh$UNAFF_mean)

# check dist
hist(running_roh$UNAFF_mean, xlim = c(0,100), breaks = 10)

# roh desert
# filter windows with too few snps (lowest percentile)
cut_window <- quantile(running_roh$UNAFF_n, probs = c(0.01))
# calculate 0.5% and 99.5% quantile 
quants <- quantile(running_roh$UNAFF_mean, na.rm = TRUE, probs= c(0.005, 0.995))
# obtain the index of deserts 
# here quants[1] = 0, we are looking for regions without ROH
desert_index <- length(running_roh[running_roh$UNAFF_mean == quants[1], ][,1])

roh_deserts <- running_roh %>% 
        filter(UNAFF_n > cut_window) %>% 
        as_tibble() %>% 
        mutate(prop_roh = UNAFF_mean/cnt_ind) %>%  # 7691 cnt_ind
        arrange(prop_roh) %>% 
        # top 0.5% of windows
        .[1:desert_index, ]
mean(roh_deserts$prop_roh)

# make table for supplementary
roh_deserts %>% 
        mutate(win_start = round(win_start/1000, 2), win_end = round(win_end/1000, 2), 
               prop_roh = round(prop_roh * 100, 2)) %>% 
        dplyr::select(CHR, win_start, win_end, prop_roh, UNAFF_n) %>%
        setNames(c("Chromosome", "WinStart", "WinEnd", "% of individuals with ROH", "N (SNPs)")) %>% 
        gt() %>% 
        tab_header(
                title = "Top 0.5% ROH deserts",
                subtitle = "ROH density measured in 500Kb running windows"
        ) %>% 
        gtsave(filename = "figs/tables/roh_desert.png")

# obtain the index of islands 
# here quants[1] = 0, we are looking for regions with a lot of ROH
islands_index <- length(running_roh[running_roh$UNAFF_mean > quants[2], ][,1])

roh_islands <- running_roh %>% 
        filter(UNAFF_n > cut_window) %>% 
        mutate(prop_roh = UNAFF_mean/cnt_ind) %>%  # 7691 cnt_ind
        arrange(desc(prop_roh)) %>% 
        .[1:islands_index, ]

mean(roh_islands$prop_roh)

# make table for supplementary
roh_islands %>% 
        mutate(win_start = round(win_start/1000, 2), win_end = round(win_end/1000, 2), 
               prop_roh = round(prop_roh * 100, 2)) %>% 
        dplyr::select(CHR, win_start, win_end, prop_roh, UNAFF_n) %>%
        setNames(c("Chromosome", "WinStart", "WinEnd", "% of individuals with ROH", "N (SNPs)")) %>% 
        gt() %>% 
        tab_header(
                title = "Top 0.5% ROH islands",
                subtitle = "ROH density measured in 500Kb running windows"
        ) %>% 
        gtsave(filename = "figs/tables/roh_islands.png")


quants <- quantile(running_roh$UNAFF_mean, na.rm = TRUE, probs= c(0.01, 0.99))

# # check that genome assembly is ok where deserts and islands are ---------------
# # load linkage map with interpolated SNP positions
# lmap <- read_delim("data/Oar3.1_Interpolated.txt", "\t") %>% 
#         setNames(c("chr", "snp_name", "bp", "cM")) %>% 
#         mutate(mb_pos = bp/1000000,
#                kb_pos = bp/1000)
# 
# win_area <- function(win_mid, CHR, ...) {
#         diffs <- lmap %>% 
#                 filter(chr == CHR) %>% 
#                 mutate(diff = abs(win_mid - kb_pos)) %>% 
#                 arrange(diff) %>% 
#                 top_n(-500)
# }
# 
# # deserts
# all_des <- pmap(roh_deserts, win_area) %>% 
#                 map(as_tibble) %>% 
#                 bind_rows(.id = "desert_num") %>% 
#                 mutate(desert = paste0(desert_num, " | Chr. ", chr),
#                        desert = fct_inorder(desert))
#  
# all_des$win_start <- rep(roh_deserts$win_start, each = 500)
# all_des$win_end <- rep(roh_deserts$win_end, each = 500)
# 
# # ROH deserts and islands tend to be very good at coinciding with genome
# # assembly errors. Let's check that this is not the case by 
# # by plotting genetic vs. physical SNP positions around deserts and islands.]
# 
# # deserts
# p_des_rec <- ggplot(all_des, aes(mb_pos, cM)) +
#         geom_point(shape = 21, fill = "#eceff4", stroke = 0.05, size = 2) +
#         facet_wrap(~desert, scales = "free") +
#         scale_x_continuous( breaks = scales::pretty_breaks(3)) +
#         scale_y_continuous( breaks = scales::pretty_breaks(3)) +
#         geom_vline(aes(xintercept = win_start/1000)) +
#         geom_vline(aes(xintercept = win_end/1000)) +
#         #theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
#         xlab("Mb") +
#         ggtitle("Genetic vs. physical SNP positions in regions with ROH deserts")
# p_des_rec 
# # ggsave("figs/pot_sup_ROH_des_rec.jpg", width = 8, height = 7)
# 
# # islands
# all_isl <- pmap(roh_islands, win_area) %>% 
#         map(as_tibble) %>% 
#         bind_rows(.id = "island_num") %>% 
#         mutate(island = paste0(island_num, " | Chr. ", chr),
#                island = fct_inorder(island))
# 
# all_isl$win_start <- rep(roh_islands$win_start, each = 500)
# all_isl$win_end <- rep(roh_islands$win_end, each = 500)
# 
# p_isl <- ggplot(all_isl, aes(mb_pos, cM)) +
#         geom_point(shape = 21, fill = "#eceff4", stroke = 0.05, size = 2) +
#         facet_wrap(~island, scales = "free") +
#         scale_x_continuous( breaks = scales::pretty_breaks(3)) +
#         scale_y_continuous( breaks = scales::pretty_breaks(3)) +
#         geom_vline(aes(xintercept = win_start/1000)) +
#         geom_vline(aes(xintercept = win_end/1000)) +
#         #theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
#         xlab("Mb") +
#         ggtitle("Genetic vs. physical SNP positions in regions with ROH islands")
# p_isl
# # ggsave("figs/pot_sup_ROH_isl_rec.jpg", width = 8, height = 7)


# save islands and deserts
roh_extremes <- bind_rows(roh_islands, roh_deserts, .id = "extreme") %>% 
                mutate(extreme = ifelse(extreme == 1, "island", "desert"))
write_delim(roh_extremes, path = "output/roh_islands_deserts.txt")

col <- c("#FF0000")

df %>% 
        filter(MB > 0) %>% 
        filter(CHR %in% 1:10) %>% 
        ggplot() +
        geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
                  alpha=0.5, fill = "slategray") + # "#f7f7f7" "#eceff4"
        geom_hline(data = yax, aes(yintercept = yax), color = "#f8f9fb", size = 0.4) +
        geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 1.6, ymax = yax + 1.8, 
                      fill = as.factor(CHR)),  col = "#FF0000", size = 0.2, alpha = 1) + 
        scale_fill_manual(values = rep(col, 10)) + 
        scale_color_manual(values = rep(col, 10)) +
        scale_y_reverse(expand = c(0, 0)) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13, base_family = "Helvetica") +
        facet_grid(~CHR,scales = 'free_x', space = 'free_x', switch = 'x',
                   labeller = as_labeller(chr_names)) +
        theme(#=strip.placement = 'outside',
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                panel.spacing = unit(0, "lines"),
                plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
                axis.line.x = element_blank(),
                legend.position="none",
                axis.title.x = element_text(margin=margin(t=0)),
                axis.title.y = element_text(margin=margin(r=0)),
                axis.text.y = element_text(colour = "white"),
                axis.line.y = element_blank()) +
        coord_cartesian(clip = 'off') +
        xlab("Chromosome") +
        ylab("Individuals") -> ROH_per_ind

ROH_per_ind +
geom_vline(data = roh_islands, aes(xintercept = win_mid/1000), color = "red", size = 1, alpha=0.2)  -> highlight

highlight 
graph2ppt(file="figs/roh_per_ind_island", width=14, height=6) 

roh_islands

####################
#  ROH counts plot #
####################

# the input is located at /Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/All_n_507_for_NSA_plot/plot/
# use 14_ROH_plink.R to generated the input from All_n_507_for_NSA_plot folder.
library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/All_n_507_for_NSA_plot/plot")

data_process <- function(headname, pop){
  #headname = "DBW1"
  name = paste0("Empirical.relatedness.",headname, ".txt")
  DT = read.delim(name, header = TRUE, sep=' ')
  dat <- data.frame(relat=DT$ritland, Population=pop) # @change
  dat$relat <- as.numeric(dat$relat)
  return(dat)
}

DT = read.delim("ROH_analyse_population_x.hom.txt", header = TRUE, sep='\t')
plotdat <- data.frame(ROH=DT$NSEG, Population=DT$FID)

wilcox_dt = read.delim("ROH_analyse_population_x.hom.wildsel.txt", header = TRUE, sep='\t')
wilcox.test(wilcox_dt$NSEG[wilcox_dt$FID == "WILD"],wilcox_dt$NSEG[wilcox_dt$FID == "SEL"])

#   MEW1        MEW2      LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
col <- c( "#0A2C86", "#325A98",  "#1D92BD", "#3DB9C1", "#C4E9B3", "#7BD2BF", "#ECF6B9", "#EEE8AA", 
          #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
          "#F9476B", "#FC709F","#E376B7", "#CF7FBC",  "#A36DC1", "#FEB22B", "#F36616", "#D83B1C", "#FF9117")
plotdat$Population = factor(plotdat$Population, levels=c("MEW1", "MEW2","LIW1","LIW2","DBW1","DBW2","NCW1","NCW2","DBX1","DBX2","DBX3","UNC1","UNC2","UMFS","NEH1","NEH2", "MEH2"))
#jpeg("Relatedness_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
p <- ggplot(plotdat, aes(x=Population, y=ROH, fill=Population)) +
  geom_boxplot(width = 0.8, outlier.size = 0.1) + 
  #geom_jitter(alpha = 0.03, width = 0.2)+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = mean, geom='point', shape=20, size=2, color='red', fill='red') +
  labs(x=NULL, y = "ROH segments (500KB/window)") 
p +  scale_fill_manual(values=col) + 
  #scale_y_continuous(limits=c(0, 0.25)) +
  theme(text = element_text(size=20),
        legend.position = "none")+
  coord_cartesian(ylim = c(0, 20))

graph2ppt(file="ROH",width=10,height=6)

