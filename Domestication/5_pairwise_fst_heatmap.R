library(ggplot2)
library(ComplexHeatmap)
library(circlize)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/05_pairwise_fst")
dt = read.delim("05_pairwise_fst.txt", header = TRUE, sep='\t')
head(dt)

plot_dt <- as.matrix(dt)
rownames(plot_dt) <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
colnames(plot_dt) <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
plot_dt <- as.matrix(plot_dt)
class(plot_dt)

fa = rep(c("WILD", "DOM"), times = c(9, 8))
fa_col = c("WILD" = "#aecdc2", "DOM" = "#f0b8b8")
dend1 = cluster_between_groups(plot_dt, fa)

Heatmap(plot_dt, name = "Pairwise Fst",
        col = colorRamp2(c(0, 0.1, 0.2), c("white", "skyblue", "orange")),
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        width = unit(14, "cm"), height = unit(8, "cm"),
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        top_annotation = HeatmapAnnotation(Group = fa, col = list(Group = fa_col)))

graph2ppt(file="Pairwise_fst",width=10,height=6)
