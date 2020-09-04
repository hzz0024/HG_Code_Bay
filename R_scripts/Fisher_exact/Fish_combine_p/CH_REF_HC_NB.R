SAMPLE = FALSE
GLOBAL_NAME = 'CH_REF_HC_NB_'
file1 = "fish_CH_REF_HCNB.txt"
file2 = "fish_HC_NB.txt"
#================================

is.odd <- function(x) x %% 2 != 0
#FUNCTION fisher.method():
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

########### Combine two p-values ##############
combine_p <- function(p1_name, p2_name){
  dat_1 <- read.delim(p1_name, header = FALSE, sep='\t')
  dat_2 <- read.delim(p2_name, header = FALSE, sep='\t')
  p1 <- dat_1$V7
  p2 <- dat_2$V7
  ps <- data.frame(p1, p2)
  return(ps)
}

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
plot.scan = function(DF, MAIN, LEFT, Y, LEN, threshold){##NOTE THIS FUNCTION IS REPEATED BELOW
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
  
  DF = DF[x==TRUE,] #@ww
  colors = colors[x==TRUE] #@ww
  CM = CM[x==TRUE]
  
  plot(-log(p, 10)~CM, data = DF, main = NULL, pch = 19, col = colors, cex = point.size, axes = F, xlab = "Chromosome", ylab = expression('-log'['10']*'(p)'), ylim=setYlim)
  axis(1, at = endpoints, labels = F)
  axis(2, labels = NULL, las = 1, cex.axis = CEX)
  mtext(lg.names, side = 1, at = mids, line = .75, cex = CEX)
  par(lend = 2)
  if (PLOT.SCALE == TRUE){
    segments(LEFT, Y, (LEFT + LEN), Y, lwd = 2)
    TEXT = paste(LEN, "cM")
    labPos = LEFT + .5 * LEN
    text(labPos, y = (Y + .4), labels = TEXT)
  }
}

#FUNCTION do.plot()
#Runs get.combined and plot.scan to build a Manhattan plot and return the modified plot data
do.plot = function(pval.df, MAIN, LEFT, Y, LEN, threshold){
  p = get.combined(pval.df) 
  plot.dat = data.frame(VARIANT, LG, CM, p)
  plot.scan(plot.dat, MAIN, LEFT, Y, LEN, threshold)
  rownames(plot.dat) = plot.dat$VARIANT
  return(plot.dat)
}

## NOW WE CAN PLOT THE MANHATTAN PLOTS FOR EACH REPLICATE
####### DO AC ONE-TAILED PLOT ##############

dat = read.delim(file1, header=FALSE, sep='\t')
ps <- combine_p(file1, file2)

if(SAMPLE==TRUE){
  idx = sample(seq(1,length(dat[,1])), 50000)
  ps = ps[idx,]
  dat = dat[idx,]
}


fisher_method(ps, dat, out1=paste0(GLOBAL_NAME, 'out_0.01.txt'), out5=paste0(GLOBAL_NAME, 'out_0.05.txt')) #write outliers to output files, no return

LG = as.numeric(factor(dat$V1))
CM = dat$V2 
CM = LG*100000000+CM
#LG = dat$lg
VARIANT = paste0(dat$V1,'_',dat$V2)
#VARIANT = dat$variant
lg.names = 1:10
point.size = .3
setYlim = c(-.2, 8)

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
PLOT.SCALE = TRUE
png(paste0(GLOBAL_NAME, 'plot_0.05.jpg'))
do.plot(ps, paste0(GLOBAL_NAME, 'PLOT'), 100, 5, 100, 0.05)
dev.off()
png(paste0(GLOBAL_NAME, 'plot_0.01.jpg'))
do.plot(ps, paste0(GLOBAL_NAME, 'PLOT'), 100, 5, 100, 0.01)
dev.off()


