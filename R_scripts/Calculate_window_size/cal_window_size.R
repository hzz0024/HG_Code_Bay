GLOBAL_NAME = 'REF-CH-SR-HC_'
file1 = "fish_REF_CH_REF-CH-SR-HC.txt"
file2 = "fish_SR_HC_REF-CH-SR-HC.txt"
#================================

## Produce datasets for manhattan plots 
dat = read.delim(file1, header=FALSE, sep='\t')
idxs = order(dat$V1, dat$V2)
dat = dat[idxs,]
#ps = ps[idxs,]

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
point.size = .3
setYlim = c(-.2, 8)

DF = data.frame(VARIANT, LG, CM)
#SET UP LABELS FOR THE 10 LINKAGE GROUPS
lgs = 1:10
lg.names = lgs

#FUNCTION setup.windows
#Divides the linkage groups into windows each with number of markers equal to grab.size
setup.windows = function(grab.size, lg.tot){
  first.whole = floor(lg.tot/grab.size)##sets up the first windows up until the last
  left.sides = c()
  right.sides = c()
  left = 0
  for (i in 1:first.whole){
    right = grab.size*i
    left.sides = append(left.sides, left)
    right.sides = append(right.sides, right)
    left = right  ##now make the new left the previous right side, this way all adjacent windows will overlap
  }
  last.right = lg.tot
  last.left = lg.tot - (grab.size-1)
  left.sides = append(left.sides, last.left)
  right.sides = append(right.sides, last.right)
  window.sides = data.frame(left.sides, right.sides)
  return(window.sides)
}

#FUNCTION plot.boot.bars
#This combines the do.plot() function used above with the bootstrapping functions to plot
#a finalized figure that shows the combined p values and color coded bars below regions
#with significant bootstrapping signal. Also outputs the bootstrap coordiantes for later use.
plot.boot.bars = function(DF, DF.ps, grab.size, cut, null.dist, MAIN, bar.Ycoord, bar.width, easy.cut, hard.cut){
  #ARGUMENTS: 
  # DF = the dataframe with the p values output for a particular cross (output from function do.plot()) 
  # DF.ps = data frame with the one-sided p values for each individual replicate
  # grab.size = the window size (in number of loci) used to build the null distribution beging used
  # cut = the cutoff p value cutoff for the null distribution
  # null.dist = the null distribution output from get.bootstrap()
  # MAIN = title
  # barYcoord = the vertical coordinate for placing the bootstrap bars
  # bar.width = the vertical width of the bootstrap bar
  # easy.cut = the unconservative p value cutoff for bootstrap bars
  # hard.cut = the conservative cutoff
  x = c()
  boot.ps = c()
  sig.counts = c()
  
  obs.count.c = c()
  for (i in lgs){
    window.left.bounds = c()
    window.right.bounds = c()
    lg.sub = DF[DF$LG == i,] ##pull out the linkage groups one at a time
    ##set up the windows. This is more complicated than just seq(by = grab.size, because we want to include every marker)
    window.sides = setup.windows(grab.size, length(lg.sub[,1]))##sets up the windows as a dataframe of left and right sides
    for (win in seq(1:(length(window.sides[,1])))){
      left = window.sides[win,1]
      right = window.sides[win,2]
      sub = lg.sub[left:right,]
      win.left.cM = sub[1,'CM']##grabs the left side of this window in cM
      win.right.cM = sub[length(sub[,1]),'CM']##grabs the right side of the window
      window.left.bounds = append(window.left.bounds, win.left.cM)
      window.right.bounds = append(window.right.bounds, win.right.cM)
    }
    w = sum(window.right.bounds - window.left.bounds)/length(window.left.bounds) 
    #print(paste0("Chr_",i))
    print(floor(w))
  }
}

for(grab.size in c(15, 25, 50, 150, 1500, 10000)){
  cut = 0.05
  ps.boot = c()
  ps = c()
  
  EASY.BAR.COLOR = 'darkgreen'
  EASY.CUT = 0.05
  HARD.BAR.COLOR = 'green'
  HARD.CUT = 0.01
  
  print(grab.size)
  plot.boot.bars(DF, ps, grab.size, cut, ps.boot, 'Test', .825*setYlim[1], -.75*setYlim[1], 0.1, 0.05)
}

