SAMPLE = TRUE
GLOBAL_NAME = 'SR_REF_COH_ARN_'
file1 = "fish_REF_CH_REF-CH-SR-HC.txt"
file2 = "fish_SR_HC_REF-CH-SR-HC.txt"
#================================
#check if odd or even number
is.odd <- function(x) x %% 2 != 0

#performs fisher's method to combine a set of p values
#ARGUMENTS: ps = a vector of p values you want combined
fisher.method=function(ps) {
  ftest=-2*sum(log(ps))
  df=length(2*ps)
  pv=1-pchisq(q=ftest,df=df)
}

#FUNCTION get.combined():
#runs fisher.method() on a dataframe of p values
#ARGUMENTS: DF = a two column data frame with p values in it
get.combined = function(DF){
  combined.p.values = c()
  alt = c()
  for (i in seq(1, length(DF[,1]))){ #(apparently did not know of nrow() yet)
    x = DF[i,1]
    y = DF[i,2]
    combined.p = fisher.method(c(x,y))
    combined.p.values = append(combined.p.values, combined.p)
  }
  return(combined.p.values)
}

# Functio to combine two p-values @HG
# ARGS: p1_name/p2_name - files with exact test p-values
combine_p <- function(p1_name, p2_name){
  dat_1 <- read.delim(p1_name, header = FALSE, sep='\t')
  dat_2 <- read.delim(p2_name, header = FALSE, sep='\t')
  p1 <- dat_1$V7
  p2 <- dat_2$V7
  ps <- data.frame(p1, p2)
  return(ps)
}

# Combined p-value adjustment and output the SNPs with fdr < 0.01 and < 0.05 @HG
# ARGS: ps - a vector of combined p-values; dat - file with exact test p-values; out1/out5 - output names for SNPs with fdr < 0.01 and fdr < 0.05
fisher_method <- function(ps, dat, out1, out5){
  
  ps = get.combined(ps)
  adj = p.adjust(ps, method = 'BH')
  idx = adj < 0.05
  dat5 = paste0(dat$V1[idx],'\t',dat$V2[idx],'\t',adj[idx])
  message(paste0('0.05: ',length(dat$V1[idx])))
  write.table(dat5, out5, quote = FALSE, row.names = FALSE)
  idx = adj < 0.01
  dat1 = paste0(dat$V1[idx],'\t',dat$V2[idx],'\t',adj[idx])
  message(paste0('0.01: ',length(dat$V1[idx])))
  write.table(dat1, out1, quote = FALSE, row.names = FALSE)
}



#======================================================================
#FUNCTION plot.scan()
#Adjusts a set of p values for FDR and plots them based on map locations and adds bars below to indicate boostrap significances
#ARGUMENTS:
# DF = the dataframe
# MAIN = the title you want for the plot
# LEFT = vector of the left sides of bootstrap bars
# Y = the Y coordinate for where to plot the bootstrap bars
# LEN = vector of lengths for the bootstrap bars
do.plot = function(pval.df, MAIN, LEFT, Y, LEN, threshold=0.05){##NOTE THIS FUNCTION IS REPEATED BELOW
  p = get.combined(pval.df)
  DF = data.frame(VARIANT, LG, CM, p)
  
  CEX = 1.2
  adj = p.adjust(DF$p, method = 'BH')
  #set up grey and black points alternating by LG
  colors = is.odd(DF$LG)
  colors[colors == TRUE] <- 'black'
  colors[colors == FALSE] <- 'grey'
  #make significant points red
  x = adj < threshold
  for (i in 1:length(x)){
    if (x[i] == TRUE){
      colors[i] <- 'red'
    }
  }
  # comment out if all points need to be plotted
  #DF = DF[x==TRUE,] #@HG
  #colors = colors[x==TRUE] #@HG
  #CM = CM[x==TRUE] #@HG
  plot(-log(p, 10)~CM, data = DF, main = NULL, pch = 19, col = colors, cex = point.size, axes = F, xlab = "Chromosome", ylab = expression('-log'['10']*'(p)'), ylim=setYlim)
  axis(1, at = endpoints, labels = F)
  axis(2, labels = NULL, las = 1, cex.axis = CEX)
  mtext(lg.names, side = 1, at = mids, line = .75, cex = CEX)
  par(lend = 2)
  if (PLOT.SCALE == TRUE){
    segments(LEFT, Y, (LEFT + LEN), Y, lwd = 2)
    TEXT = paste(LEN, "cM")
    labPos = LEFT + .5 * LEN
    #text(labPos, y = (Y + .4), labels = TEXT)
  }
  rownames(DF) = DF$VARIANT
  #return(DF)
}

## NOW WE CAN PLOT THE MANHATTAN PLOTS FOR EACH REPLICATE

dat = read.delim(file1, header=FALSE, sep='\t')
ps <- combine_p(file1, file2)

# decide whether use test sample or not
if(SAMPLE==TRUE){
  idx = sample(seq(1,length(dat[,1])), 50000)
  ps = ps[idx,]
  rownames(ps) <- seq(length=nrow(ps))
  dat = dat[idx,]
  rownames(dat) <- seq(length=nrow(dat))
}


idxs = order(dat$V1, dat$V2)
dat = dat[idxs,]
ps = ps[idxs,]

fisher_method(ps, dat, out1=paste0(GLOBAL_NAME, 'out_0.01.txt'), out5=paste0(GLOBAL_NAME, 'out_0.05.txt')) #write outliers to output files, no return

pos = dat$V2
chr = as.numeric(factor(dat$V1))
posmin <- tapply(pos,chr, min);
posmax <- tapply(pos,chr, max);
posshift <- head(c(0,cumsum(posmax)),-1);
names(posshift) <- levels(chr)
genpos <- pos + posshift[chr];

LG = as.numeric(factor(dat$V1))
CM = genpos
VARIANT = paste0(dat$V1,'_',dat$V2)
lg.names = 1:10
point.size = .3
setYlim = c(-.2, 8)
#setXlim = c(0, 10)

#SET UP LABELS FOR THE 14 LINKAGE GROUPS
lgs = 1:10

#GATHER THE ENDPOINTS AND MIDPOINTS OF EACH LINKAGE GROUP FOR PLOTTING PURPOSES
endpoints = c(0)
for (i in lgs){
  sub = CM[LG == i]
  end = max(sub)
  endpoints = append(endpoints, end)
}
mids = c()
for (i in seq(1:(length(endpoints) - 1))){
  left = endpoints[i]
  right = endpoints[(i + 1)]
  point = (left + right) / 2
  mids = append(mids, point)
}

head(ps)
par(mfrow = c(1,1))
PLOT.SCALE = FALSE
png(paste0(GLOBAL_NAME, 'plot_0.05.jpg'))
do.plot(ps, paste0(GLOBAL_NAME, 'PLOT'), -1 , -1 , -1, 0.05)
dev.off()
png(paste0(GLOBAL_NAME, 'plot_0.01.jpg'))
do.plot(ps, paste0(GLOBAL_NAME, 'PLOT'), 1, 5, 1, 0.01)
dev.off()
