#install_github('tavareshugo/windowscanr')
library(windowscanr)
source("manhattan.R")
####################################
##########  read input     ########
####################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/11_SGS/11_SGS_boot_window")
pname = "ps_Del19_challenge.txt"
dat = read.delim(pname, header = FALSE, sep='\t',col.names=c('chromo', 'position', 'p1', 'p0', 'delta_p', 'ps', 'raw_candidates'))
dat = dat[order(dat[,1], dat[,2]),]
dat$ps[dat$ps == 0] = 0.00001
dat$adj = p.adjust(dat$ps, method = 'BH')

#check if odd or even number
is.odd <- function(x) x %% 2 != 0
# function to combine the p-value, multi options are available (this is the new method with correct df estimate)

SAMPLE = FALSE
if(SAMPLE==TRUE){
  #idx = seq(1,50000) # sampling first 50000
  idx = sample(seq(1,length(dat[,1])), 5000)
  dat = dat[idx,]
  rownames(dat) <- seq(length=nrow(dat))
}

dat = dat[order(dat[,1], dat[,2]),]

pos = dat$position
chr = as.numeric(factor(dat$chromo))
posmin <- tapply(pos,chr, min);
posmax <- tapply(pos,chr, max);
posshift <- head(c(0,cumsum(posmax)),-1);
names(posshift) <- levels(chr)
genpos <- pos + posshift[chr];

reverse_pos <- function(genpos){
  chr = sum(genpos > posshift)
  pos = genpos - posshift[chr]
  return(paste0(chr,'_',pos))
}

LG = as.numeric(factor(dat$chromo))
CM = genpos
VARIANT = paste0(dat$chromo,'_',dat$position)
lg.names = 1:10
lgs = lg.names
point.size = .3
setYlim = c(-.2, 2)
p = dat$ps
#@change: adjust p
#p = p.adjust(p, 'BH')

DF = data.frame(VARIANT, LG, CM, p)

#SET UP LABELS FOR THE 10 LINKAGE GROUPS
lgs = lg.names

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

do.plot = function(DF, threshold=0.05){##NOTE THIS FUNCTION IS REPEATED BELOW
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
  plot(-log(adj, 10)~CM, data = DF, main = NULL, pch = 19, col = colors, cex = point.size, axes = F, xlab = "Chromosome", ylab = expression('-log'['10']*'(fdr)'), ylim=setYlim)
  axis(1, at = endpoints, labels = F)
  axis(2, labels = NULL, las = 1, cex.axis = CEX)
  mtext(lg.names, side = 1, at = mids, line = .75, cex = CEX)
  par(lend = 2)
  rownames(DF) = DF$VARIANT
  #return(DF)
}

do.plot(DF, threshold = 0.05)




reps = 10^5

get.bootstrap = function(DF, id, grab.size, cut){
  counts = c()
  ps = c()
  for (i in 1:reps){
    sub = sample(DF$p, grab.size)
    sub2 = sub[sub < cut]
    new.count = length(sub2)
    new.p = mean(sub)
    counts = append(counts, new.count)
    ps = append(ps, new.p)
  }
  results = data.frame(counts, ps)
  write.table(results, paste(GLOBAL_NAME, id, grab.size, cut, sep = "_"), row.names = F, quote = F)
  return(results)
}##Arugments: DF = the dataframe; id = a string to id the file from when exported

#SET UP THE FUNCTIONS FOR PERFORMING THE BOOTSTRAPPING

# FUNCTION get.boot.p
# A putatively selected variant is one with an unadjustd p value 0.05
# This function compares a count of putatively selected variants
# To the null distribution and returns the proportion of null entires 
# with equal more greater number of putatively selected variants
get.boot.p = function(null, obs){
  tot = length(null)
  #more.extreme.count = length(null[null <= obs]) # @HG change
  more.extreme.count = length(null[null >= obs])
  p.value = more.extreme.count/tot
  print(p.value)
  return(p.value)
}

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
  window.left.bounds = c()
  window.right.bounds = c()
  for (i in lgs){
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
      #get p value for counts of significant markers (turn this on and off)
      sig.sub = sub$p[sub$p < cut]
      obs.count = length(sig.sub)
      p = get.boot.p(null.dist, obs.count)
      boot.ps = append(boot.ps, p)
      cM = median(sub$CM)
      x = append(x, cM)
    }
  }
  #x=x[1:51]
  #window.left.bounds = window.left.bounds[1:51]
  #window.right.bounds = window.right.bounds[1:51]
  #boot.ps = boot.ps[1:51]
  result1 = data.frame(x, boot.ps, window.left.bounds, window.right.bounds)
  #result1$adj.ps = p.adjust(result1$boot.ps, method = "BH")
  result1$adj.ps = result1$boot.ps
  #write.table(result1, paste(GLOBAL_NAME, grab.size, cut, sep = "_"), row.names = F, quote = F)
  easy = result1[result1$adj.ps > hard.cut,]
  easy = easy[easy$adj.ps < easy.cut,]
  hard = result1[result1$adj.ps < hard.cut,]
  print(head(easy))
  print(nrow(easy))
  do.plot(DF.ps, MAIN, 100, 5, 100, threshold = 0.05)
  par(xpd=TRUE)
  top = bar.Ycoord + 0.5*bar.width
  bottom = bar.Ycoord - 0.5*bar.width
  print(seq(bottom, top, by = 0.01))
  for (y in seq(bottom, top, by = 0.01)){
    if (nrow(easy) >= 1){
      for (i in 1:nrow(easy)){
        segments(easy$window.left.bounds[i], y, easy$window.right.bounds[i], y, col = EASY.BAR.COLOR, lwd = 2)
      }
    }
    if (nrow(hard) >= 1){
      for (i in 1:nrow(hard)){
        segments(hard$window.left.bounds[i], y, hard$window.right.bounds[i], y, col = HARD.BAR.COLOR, lwd = 2)
      }
    }
  }
  par(xpd=FALSE)
  if (nrow(hard) > 0 || nrow(easy) > 0){
    print("easy:")
    print(easy)
    print("hard:")
    print(hard)
    COLNAMES = c('lefts', 'rights', 'cut')
    
    if(nrow(easy) > 0){
      e.lefts = easy$window.left.bounds
      e.rights = easy$window.right.bounds
      output.e = data.frame(e.lefts, e.rights)
      output.e$cut = "easy"
      if(nrow(hard) == 0){
        colnames(output.e) = COLNAMES
        output = output.e
      }
    }
    
    if(nrow(hard) > 0){
      h.lefts = hard$window.left.bounds
      h.rights = hard$window.right.bounds
      output.h = data.frame(h.lefts, h.rights)
      output.h$cut = "hard"
      if(nrow(easy) == 0){
        colnames(output.h) = COLNAMES
        output = output.h
      }
    }
    
    if(nrow(hard) > 0 && nrow(easy) > 0){
      colnames(output.e) = COLNAMES
      colnames(output.h) = COLNAMES
      output = rbind(output.e, output.h)
    }
    
    #return(output)
  }
  #return(result1)
  return(list(result1, output))
}





























mean_na_rm <- function(x, ...){
  mean(x, na.rm = TRUE, ...)
}


pos_win <- winScan(x = dat, 
                   groups = "chromo", 
                   position = "position", 
                   values = c("ps"), 
                   win_size = 500,
                   win_step = 250,
                   funs = c("mean_na_rm"))

pos_win_noNA <- subset(pos_win, !is.na(ps_mean_na_rm))
pos_win_noNA
plot(-log(pos_win_noNA$ps_mean_na_rm, 10))
set.seed(1)
##step1. get the null distributions for counts of markers below a p theshold
##note that the p-values in the dataframe are not transformed (i.e. not fdr corrected), and these are used to get the window p values
get.bootstrap = function(DF, id, reps, cut){
  chromo = c()
  win_start = c()
  win_end = c()
  win_mid = c()
  ps_n = c()
  ps_mean_na_rm = c()
  ps = c()
  for (i in 1:dim(DF)[1]){
    chr_ = DF$chromo[i]
    ws = DF$win_start[i]
    we = DF$win_end[i]
    wm = DF$win_mid[i]
    snp_cnt = DF$ps_n[i]
    ps_mean = DF$ps_mean_na_rm[i]
    sub = sample(DF$ps_mean_na_rm, reps)
    obs = DF$ps_mean_na_rm[i]
    extreme.count = length(sub[sub < obs])
    p.value = extreme.count/reps
    chromo = c(chromo, chr_)
    win_start = c(win_start, ws)
    win_end = c(win_end, we)
    win_mid = c(win_mid, wm)
    ps_n = c(ps_n, snp_cnt)
    ps_mean_na_rm = c(ps_mean_na_rm, ps_mean)
    ps = c(ps, p.value)
    # if(p.value < cut){
    #   chromo = c(chromo, chr_)
    #   win_start = c(win_start, ws)
    #   win_end = c(win_end, we)
    #   win_mid = c(win_mid, wm)
    #   ps_n = c(ps_n, snp_cnt)
    #   ps_mean_na_rm = c(ps_mean_na_rm, ps_mean)
    #   ps = c(ps, p.value)
    # }
  }
  results = data.frame(chromo,win_start, win_end, win_mid, ps_n, ps_mean_na_rm, ps)
  #write.table(results, paste(GLOBAL_NAME, id, reps, cut, sep = "_"), row.names = F, quote = F)
  return(results)
}##Arugments: DF = the dataframe; id = a string to id the file from when exported

##step1. get the null distributions for counts of markers below a p theshold
##note that the p-values in the dataframe are not transformed (i.e. not fdr corrected), and these are used to get the window p values

set.seed(1)
dat <- get.bootstrap(pos_win_noNA, "test", 50, 0.05)
dat1 <- get.bootstrap(pos_win_noNA, "test", 50, 0.05)

dat$SNP = paste0(dat$chromo,'_',dat$win_mid)
dat1$SNP = paste0(dat1$chromo,'_',dat1$win_mid)
length(intersect(dat$SNP, dat1$SNP))

# replace chromosome if it is numerical
chr_num_list = seq(10)
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(1:10)) 
  pos_win_noNA$chromo[pos_win_noNA$chromo==chr_str_list[i]] = chr_num_list[i]
pos_win_noNA$chromo = as.numeric(pos_win_noNA$chromo)
pos_win_noNA$SNP = paste0(pos_win_noNA$chromo,'_',pos_win_noNA$position)
is.odd <- function(x) x %% 2 != 0
colors = is.odd(pos_win_noNA$chromo)

colors[colors == TRUE] <- 'black'
colors[colors == FALSE] <- 'grey'

#make significant points red
threshold = 0.05
x = dat$ps < threshold
for (i in 1:length(x)){
  if (x[i] == TRUE){
    colors[i] <- 'red'
  }
}


manhattan(chr="chromo",bp="win_mid",p="ps_mean_na_rm", snp = "SNP", pos_win_noNA, logp=TRUE, cex.axis = 1.2,
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('Surv - Ref '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 1, col = "black", cex.lab=1.5)


plot(-log(ps_mean_na_rm, 10)~win_mid, data = pos_win_noNA, main = NULL, pch = 19, col = colors, cex = 0.2, axes = F, xlab = "Chromosome", ylab = expression('-log'['10']*'(fdr)'))
axis(1, at = 1000000, labels = F)
axis(2, labels = NULL, las = 1, cex.axis = 1)

axis(2, labels = NULL, las = 1, cex.axis = 1)

mtext(lg.names, side = 1, at = mids, line = .75, cex = 1)

segments(0, -0.1, 100000, -0.1, lty = 1, col = "red", lwd = 8)
