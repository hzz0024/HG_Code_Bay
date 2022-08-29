# Usage: Rscript -i infile.covar -c component1-component2 -a annotation.file -o outfile.pdf --x_min minimum_x_value --x_max maximum_x_vlaue --y_min minimum_y_value --y_max maximum_y_value
# Rscript --verbose plotPCAngsd_mod_label.R -i wild_227_D100maxD450_minQ20_minMAF05_SNPe6_no227inv.cov.npy -c 1-2 -a wild_227.txt -o wild_227_D100maxD450_minQ20_minMAF05_SNPe6_no227inv.PCAngsd.WGS.pc1-2.pdf --x_min 0 --x_max 1 --y_min 0 --y_max 1
# change .txt filenames in red below
library(methods)
library(optparse)
library(ggplot2)
library(RcppCNPy)
library(ggrepel)
library(wesanderson)
setwd('.')

# colorblind friendly: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- wes_palette("Zissou1", 18, type = "continuous")
#cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for fills, add: scale_fill_manual(values=cbPalette)
# To use for line and point colors, add: scale_colour_manual(values=cbPalette)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file (output from ngsCovar)'),
                    make_option(c('-c','--comp'), action='store', type='character', default=1-2, help='Components to plot'),
                    make_option(c('-a','--annot_file'), action='store', type='character', default=NULL, help='Annotation file with individual classification (2 column TSV with ID and ANNOTATION)'),
                    make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file'),
                    make_option(c('--x_min'), action='store', type='double', default=0.0, help='X-axis minimum'),
                    make_option(c('--x_max'), action='store', type='double', default=0.0, help='X-axis maximum'),
                    make_option(c('--y_min'), action='store', type='double', default=0.0, help='Y-axis minimum'),
                    make_option(c('--y_max'), action='store', type='double', default=0.0, help='Y-axis maximum')
)
opt <- parse_args(OptionParser(option_list = option_list))

# Annotation file is in plink cluster format

#################################################################################

# Read input file
covar <- npyLoad(opt$in_file) #%>%
#covar <- read.table(c, stringsAsFactors=F, sep="\t", fill=TRUE);

# Read annot file
annot <- read.table(opt$annot_file, sep="\t", header=T); # note that plink cluster files are usually tab-separated

# Parse components to analyze
comp <- as.numeric(strsplit(opt$comp, "-", fixed=TRUE)[[1]])

# Eigenvalues
eig <- eigen(covar, symm=TRUE);

#write.table(eig$val, file = "eigen_values_pcangsd_wild_D100maxD450_minQ20_minMAF05_SNPe6_no227invSNPs_pc1-2.txt", quote = FALSE)
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Pop <- factor(annot$CLUSTER)

#sink('eigenvectors_eigen_scores_pcangsd_wild_D100maxD450_minQ20_minMAF05_SNPe6_no227invSNPs_pc1-2.txt')
print(PC)
sink()

title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")

x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")


idx_in_range = PC[x_axis]>opt$x_min & PC[x_axis]<opt$x_max & PC[y_axis]>opt$y_min & PC[y_axis]<opt$y_max
print(idx_in_range)

ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop")) + 
  ggtitle(title) + 
  scale_colour_manual(values=cbPalette, breaks=c("HC_18", "ARN_18", "COH_18", "SR_18", "NB_18", "HC_19", "ARN_19", "COH_19", "SR_19", "NB_19", "HC_21", "ARN_21", "COH_21", "SR_21", "BS_21", "BEN_21", "NAN_21", "NB_21")) + 
  # add label for each sample, edit the nudge_x and nudge_y to adjust the lable position
  geom_text_repel(data=PC[idx_in_range,], aes_string(x=x_axis, y=y_axis, color="Pop"), label=annot$IID[idx_in_range], size=2, nudge_x = 0.01, nudge_y = 0.01)

ggsave(opt$out_file)
unlink("Rplots.pdf", force=TRUE)
