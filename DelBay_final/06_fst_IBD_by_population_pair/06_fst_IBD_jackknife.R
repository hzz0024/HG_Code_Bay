setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/06_pairwise_fst_IBD/WILD21")
rm(list=ls())

fst = read.csv("WILD21.csv", row.names=1)

geo = read.csv("Pop_distance.csv", row.names=1)
# linearize fst, remove diagonals
fst[upper.tri(fst, diag=T)] = NA
fstlin = fst/(1-fst)

geo[upper.tri(geo, diag=T)] = NA
i = is.na(geo)
fst[i] = NA

WILD21_fst = as.matrix(fstlin)
WILD21_geo = as.matrix(geo)

library(vegan)
library(smatr)

sma(Fst ~ Distance, data=neutral, slope.test=0)
l = line.cis(y=WILD21_fst[lower.tri(WILD21_fst)], x= WILD21_geo[lower.tri(WILD21_geo)])
##### WILD21 and WILD21 separately with RMA regression lines from smatr
library(smatr)
quartz(width=9, height=5)
#par(mfrow=c(1,2)) # for paper?
par(mfrow=c(1,1), cex=1.6, mar = c(4,3,3,0.5), mgp=c(2,0.5,0), omi=c(0,0,0,0), cex.axis=0.9) # for ppt
pchs = 20
lwds = 1
xlims = c(0,32)
ylim = c(min(WILD21_fst, WILD21_fst, na.rm=T), max(WILD21_fst, WILD21_fst, na.rm=T))
plot(WILD21_geo[1:100], WILD21_fst[1:100], pch=pchs, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="WILD21", ylim=ylim, xlim = xlims, bty="n")
l = line.cis(y=WILD21_fst[lower.tri(WILD21_fst)], x= WILD21_geo[lower.tri(WILD21_geo)])
abline(a=l$coef[1], b=l$coef[2], col="black", lwd=lwds) # RMA linear model estimate

mantel(WILD21_geo, WILD21_fst, method="pearson", permutations=10000)

#########################
## homegrown mantel test to work with missing values
## After Piepho, HP 2005 Permutation tests for the correlation among genetic distances and measures of heterosis. Theor Appl Gen 111:95-99

mantel.me = function(fst, geo, nperm = 1000, use = "pairwise.complete.obs"){
	fstdist = as.dist(fst)
	fstsym = as.matrix(fstdist) # convert to symmetic for permuting
	r = cor(fst[lower.tri(fst)], geo[lower.tri(geo)], use="pairwise.complete.obs") # observed r
	
	n = dim(fst)[1]
	rp = numeric(0)
	for(i in 1:nperm){
		bootn = sample(1:n, n, replace=FALSE)
		bootfst = fstsym[bootn,bootn] # permut rows and columns
		rp = c(rp, cor(bootfst[lower.tri(fst)], geo[lower.tri(geo)], use=use))
	}
	p = (sum(abs(rp)>=abs(r))+1)/(nperm+1)
	out = list(r = r, p=p)
	return(out)
}

####################################
#### Jackknife over populations
####################################
# See "Arlequin jackknife/" for jackknife over loci

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/06_pairwise_fst_IBD/WILD21")
fst = read.csv("WILD21.csv", row.names=1)
geo = read.csv("Pop_distance.csv", row.names=1)
# linearize fst, remove diagonals
fst[upper.tri(fst, diag=T)] = NA
fstlin = fst/(1-fst)
geo[upper.tri(geo, diag=T)] = NA
i = is.na(geo)
fst[i] = NA

WILD21_fst = as.matrix(fstlin)
WILD21_geo = as.matrix(geo)

library(vegan)
library(smatr) # for Reduced Major Axis regression

maxpops = 8
mant = data.frame(WILD21_b = numeric(maxpops), WILD21_rmab = numeric(maxpops), WILD21_r = numeric(maxpops), WILD21_r2 = numeric(maxpops), 
WILD21_p = numeric(maxpops), WILD21_b = numeric(maxpops), WILD21_rmab = numeric(maxpops), WILD21_r = numeric(maxpops), WILD21_r2 = numeric(maxpops), WILD21_p = numeric(maxpops)) # OLS slope, RMA slope, r, r2, and Mantel p-values for each jacknife

## Plots and regression jackknife (mantel is later)
# WILD21
quartz(width=8, height=6)
par(mfrow=c(4,2))
popnames = c("HC", "ARN", "COH", "SR", "BS", "BEN", "NAN", "NB")
pops = 1:length(popnames)
for(i in pops){
	jackpops = pops[-i]
	jackgen = WILD21_fst[jackpops, jackpops]
	jackgeo = WILD21_geo[jackpops, jackpops]

	plot(jackgeo[1:length(jackgeo)], jackgen[1:length(jackgen)], pch=16, main=paste("W/out", popnames[i]), xlab="Distance (km)", ylab="Fst/(1-Fst)")
	l = lm(jackgen[lower.tri(jackgen)] ~ jackgeo[lower.tri(jackgeo)]) # OLS regression
	#lines(jackgeo[lower.tri(jackgen)], l$fitted.values, col="orange")

	mant$WILD21_b[i] = l$coefficients[2]
	mant$WILD21_r2[i] = summary(l)$r.squared	
	
	l = line.cis(y = jackgen[lower.tri(jackgen)], x= jackgeo[lower.tri(jackgeo)], alpha=0.05, method="SMA", intercept=TRUE)
	mant$WILD21_rmab[i] = l$coef[2]
	abline(a=l$coef[1], b = l$coef[2], lty=2, col="red")
}


## Mantels
# WILD21
pops = 1:8
for(i in pops){
	jackpops = pops[-i]
	jackgen = WILD21_fst[jackpops, jackpops]
	jackgeo = WILD21_geo[jackpops, jackpops]

	l = mantel(jackgeo, jackgen, method="pearson", permutations = 10000)
	mant$WILD21_p[i] = l$signif
	mant$WILD21_r[i] = l$statistic
}

row.names(mant) = paste("No", 1:8)

# Jackknife analysis on WILD21: See Sokal & Rolf

	# RMA
St = line.cis(y=WILD21_fst[lower.tri(WILD21_fst)], x= WILD21_geo[lower.tri(WILD21_geo)])$coef[2] # the observed value RMA
St
pops = 1:8
len = length(pops)
ps = numeric(len)
for(i in pops){
	ps[i] = len*St-(len-1)*mant$WILD21_rmab[i] # the pseudovalues
}
Sthat = mean(ps) # the jackknifed mean
Sthat
sst = sqrt(sum((St-ps)^2)/(len*(len-1))) # the jackknifed standard error
sst
ts = (St-0)/sst # the t-value
dt(ts, df=(len-1)) # the p-value: a mantel test is more appropriate than this, right?

write.csv(mant, paste("MantelPops",Sys.Date(),".csv", sep=""))

########################
## Can we estimate De and sigma from slope and intercept?
########################

# Mean values across WILD21 and WILD21
a = -0.003370571 # intercept = A1/(4Nsigma) Rousset 1997 Eq. 7
b = 3.009e-05 # slope = 1/(4Nsigma^2) Rousset 1997 Eq. 7

# WILD21
a =  -2.976235e-03
b = 1.910388e-05

# WILD21
a = -3.764907e-03
b =  1.022595e-04

## combine and solve the two parts of Eq. 7: sigma = a/(b*A1) and N = A1/(4asigma) = 1/(4bsigma^2)

A1 = -0.8238 # for Gaussian dispersal, from Sawyer 1977 Adv Appl Prob Eq. 2.4, as referenced by Rousset 1997 

s = a/(b*A1)
s # 135 km (mean), 189 km (WILD21), 47 km (WILD21)
N = A1/(4*a*s)
N # 0.45 (mean), 0.37 (WILD21), 1.22 (WILD21)



#################################################################
## Resampling approach to estimate 95% CI on dispersal distance
#################################################################
library(locfit)

# From Rousset method: slope = 1/4Dsigma2
len = 10000 # how many iterations?

# RMA slopes and 95% CIs for calculation of log-normal SE
slopes = data.frame(region=c("WILD21", "WILD21"), m = c(0.847e-4, 1.89e-4), l95 = c(0.63018e-4, 1.6007e-4), u95 = c(1.1375e-4, 2.2364e-4), sdlog = c(NA, NA))
for(i in 1:2){
	j = log(slopes$m[i])-log(slopes$l95[i])
	k = log(slopes$u95[i])-log(slopes$m[i])
	print(paste(j,k))
	slopes$sdlog[i] = mean(c(j,k)/1.96)
}

# Bootstrap distributions of MNe medians for use in resampling
library(boot)
samplemedian <- function(x, d) { # have to set up a function that boot can use (uses indices d)
	return(median(x[d]))
}

mne = read.csv("~/Documents/Stanford/Philippines/2008/Genotyping/Analysis/APCL/090507/MNE/Aclarkii_2009-05-13_MLNEfocalTopTwo/MLNE Summary 090513.csv")
mne$region[mne$Pop<=17] = "WILD21"
mne$region[mne$Pop>=18] = "WILD21"
cbootall = boot(as.numeric(mne$ML_Ne[mne$region=="WILD21"])/25, statistic=samplemedian, R=len, sim="ordinary")
lbootall = boot(as.numeric(mne$ML_Ne[mne$region=="WILD21"])/25, statistic=samplemedian, R=len, sim="ordinary")

mne = read.csv("~/Documents/Stanford/Philippines/2008/Genotyping/Analysis/APCL/090507/MNE/Aclarkii_2009-06-29_MLNEfocalTopTwoSurrTwo/MLNE-Flanking 090629.csv") # MNe-Flanking method
mne$region[mne$Pop<=17] = "WILD21"
mne$region[mne$Pop>=18] = "WILD21"
cbootflanking = boot(as.numeric(mne$ML_Ne[mne$region=="WILD21"])/25, statistic=samplemedian, R=len, sim="ordinary")
lbootflanking = boot(as.numeric(mne$ML_Ne[mne$region=="WILD21"]/25), statistic=samplemedian, R=len, sim="ordinary")



# Dataframe of means and sd for m, D, Ne/N, and De
	# Using RMA linear model regressions from R with approximate SE from 31.8% CI using line.cis in smatr
#priors = list(names=c("WILD21 Demographic", "WILD21 Demographic", "WILD21 MNE-All", "WILD21 MNE-All", "WILD21 MNE-Flanking", "WILD21 MNe-Flanking"), m = c(1.89e-4, 0.847e-4, 1.89e-4, 0.847e-4, 1.89e-4, 0.847e-4), mse = c(0.318e-4, 0.127e-4, 0.318e-4, 0.127e-4, 0.318e-4, 0.127e-4), D = c(144, 317, NA, NA, NA, NA), Dse = c(96, 107, NA, NA, NA, NA), nen = c(0.57, 0.57, NA, NA, NA, NA), nense = c(0.065, 0.065, NA, NA, NA, NA), De = c(NA, NA, 4.0, 3.7, 13, 21), sigma = rep(NA, 6), sigmas = vector("list", 6)) 

	# Using RMA linear model regressions from R with approximate SE from log-normal CI using line.cis in smatr
priors = list(names=c("WILD21 Demographic", "WILD21 Demographic", "WILD21 MNE-All", "WILD21 MNE-All", "WILD21 MNE-Flanking", "WILD21 MNe-Flanking", "WILD21 De=1000", "WILD21 De=1000", "WILD21 nen=10e-5", "WILD21 nen=10e-5", "WILD21 MNe-All boot", "WILD21 MNe-All boot", "WILD21 MNe-Flanking boot", "WILD21 MNe-Flanking boot"), 
m = slopes$m[c(2,1,2,1,2,1,2,1,2,1,2,1,2,1)], # IBD slopes
msdlog = slopes$sdlog[c(2,1,2,1,2,1,2,1,2,1,2,1,2,1)], # IBD slope sdlog (for use in rlnorm())
D = c(144, 317, NA, NA, NA, NA, 1000, 1000, 144, 317, NA, NA, NA, NA), # census density
Dse = c(96, 107, NA, NA, NA, NA, 0, 0, 96, 107, NA, NA, NA, NA), # SE of D
nen = c(0.57, 0.57, NA, NA, NA, NA, 1, 1, 10e-5, 10e-5, NA, NA, NA, NA), # Ne/N ratio
nense = c(0.065, 0.065, NA, NA, NA, NA, 0, 0, 0, 0, NA, NA, NA, NA), # SE of nen
De = c(NA, NA, 4.0, 3.7, 13, 21, NA, NA, NA, NA, 4, 3.7, 13, 21), # effective density
Dese = list(NA, NA, 0, 0, 0, 0, NA, NA, NA, NA, list(lbootall$t), list(cbootall$t), list(lbootflanking$t), list(cbootflanking$t)), # SE of De
sigma = rep(NA, 14), # the point estimate
sigmas = vector("list", 14), # the sampled estimates
Des = vector("list", 14)) 

for(i in 1:length(priors$names)){
	if(is.na(priors$De[i])){
		D = rnorm(len, mean = priors$D[i], sd = priors$Dse[i])
		nen = rnorm(len, mean=priors$nen[i], sd = priors$nense[i])
		De = D*nen
		priors$Des[i] = list(De)
		priors$sigma[i] = sqrt(1/(4*priors$D[i]*priors$nen[i]*priors$m[i])) # the point estimate
	} else {
		if(is.numeric(priors$Dese[[i]])){ # if we have a number for SE of De
			De = rnorm(len, mean=priors$De[i], sd = priors$Dese[[i]])
		}
		if(is.list(priors$Dese[[i]])){ # if we have a list of bootstrapped De values to approximate SE of De
			De = unlist(priors$Dese[[i]])
		}
		priors$Des[i] = list(De)
		priors$sigma[i] = sqrt(1/(4*priors$De[i]*priors$m[i])) # the point estimate
	}
	m = rlnorm(len, mean = log(priors$m[i]), sd = priors$msdlog[i])
	sigma = sqrt(1/(4*De*m))
	priors$sigmas[i] = list(sigma) 
}

# Print results
for(i in 1:length(priors$names)){
	temp = priors$sigmas[[i]]
	print(paste(priors$names[i], "  Point:", round(priors$sigma[i], digits=1), "  Mean:", round(mean(temp, na.rm=T), digits=1), "  Median:", round(median(temp, na.rm=T), digits=1), " Sd:", round(sd(temp, na.rm=T), digits=1), "  L95:", round(quantile(temp, 0.025, na.rm=T), digits=1), "  U95:", round(quantile(temp, 0.975, na.rm=T), digits=1)))
}

# Plot
par(mfrow=c(3,5))
numbreaks = 100
breaks = seq(-1, 5, length.out=numbreaks)
for(i in 1:14){
	temp = priors$sigmas[[i]]
	hist(log10(temp), breaks=breaks, xaxt="n", main=priors$names[i], col="grey", border="grey")
	labinds = seq(1,numbreaks, length.out=7)
	axis(1, at=breaks[labinds], labels=round(10^breaks[labinds]))
}








###################################################################
### Plot IBD sigma vs. De with best estimates as an illustration
###################################################################
	rousset = function(D,m){
		s = sqrt(1/(4*D*m))
		return(s)	
	}

	# Island Mean
	#m = 3.009e-05 # slope of IBD using WILD21 & WILD21 together (neg Fsts remain), removed APCL285,286 and fixed APCL250, using 090514 GEarth distances
	3mse = 1.441e-5 # standard error on the mean slope
	#mWILD21 = 1.910388e-05 # ordinary least squares
	#mWILD21 = 1.022595e-04
	mWILD21 = 0.847e-04 # RMA regression
	mWILD21 = 1.89e-04

	# Effective density estimates
	#Dlow = 0.16 # spp density from He of Neighborhoods (mu = 5e-4, 24000km, Theta=2.87)
	#Dlow = 3.17 # spp density from Ne = theta/(4mu) divided by study area (452.7 km)
	Dhigh = 252 # mean WILD21-WILD21 adult density
	#Dhigh = 182 # mean WILD21-WILD21 adult density * Ne/N from Vk = 3.532 (Jones et al. 2005)

	#Dmlnelow = 10 # ave across both regions
	#Dmlnehigh = 50

	#DmlnelowWILD21 = 8.1 # MNe-All (mean across sites)
	#DmlnelowWILD21 = 11.8
	#DmlnehighWILD21 = 38.8 # MNe-Flanking (mean across sites)
	#DmlnehighWILD21 = 73.3
	
	DmlnelowWILD21 = 3.7 # MNe-All (median across sites)
	DmlnelowWILD21 = 4.0
	DmlnehighWILD21 = 21.0 # MNe-Flanking (median across sites)
	DmlnehighWILD21 = 13.1

	DreproWILD21 =  317*0.57 # census * Ne/N from Vk from Jones et al 2005 and Planes et al 2009
	DreproWILD21 = 144*0.57

	# Plot params
	Drange = 10^seq(-3, 4, by=0.1) # total De range to plot
	sigmarange = rousset(Drange, m) # the sigma to match Derange
	sigmarangelow = rousset(Drange, m+mse)
	sigmarangehigh = rousset(Drange, m-mse)
	xticksm = c(1:10*rep(c(.001, .01, .1, 1, 10, 100, 1000, 10000),c(10, 10, 10, 10, 10, 10, 10, 10)))
	xticklg = c(0.001, 0.01, 0.1,1,10,100,1000, 10000)
	yticksm = c(1:10*rep(c(0.1, 1, 10, 100, 1000),rep(10, 5)))
	yticklg = c(0.1, 1,10,100,1000)
	ylims = c(0.1, 3000)
	lwd1 = 2
	lwd2 = 2
	col1 = "black" # for paper
	col2 = "gray60"
	#col1 = "blue" # for ppt
	#col2 = "red"
	lty1 = 2
	lty2 = 3
	lty3 = 4

	# Dispersal estimates
	rousset(DmlnelowWILD21, mWILD21)
	rousset(DmlnelowWILD21, mWILD21)
	rousset(DmlnehighWILD21, mWILD21)
	rousset(DmlnehighWILD21, mWILD21)
	rousset(DreproWILD21, mWILD21)
	rousset(DreproWILD21, mWILD21)
	rousset(3.7, mWILD21) # median of MLNE using all pops
	rousset(4.0, mWILD21) # median of MLNE using all pops
	rousset(21.0, mWILD21) # median of MLNE using two surrounding pops
	rousset(13.9, mWILD21) # median of MLNE using two surrounding pops

	median(c(rousset(DmlnelowWILD21, mWILD21), rousset(DmlnelowWILD21, mWILD21), rousset(DmlnehighWILD21, mWILD21), rousset(DmlnehighWILD21, mWILD21), rousset(DreproWILD21, mWILD21), rousset(DreproWILD21, mWILD21)))
	mean(c(rousset(DmlnelowWILD21, mWILD21), rousset(DmlnelowWILD21, mWILD21), rousset(DmlnehighWILD21, mWILD21), rousset(DmlnehighWILD21, mWILD21), rousset(DreproWILD21, mWILD21), rousset(DreproWILD21, mWILD21)))

	
	# Plot using mean slope
	plot(Drange, sigmarange, bty="l", log="xy", lwd=2, type="l", xaxt="n", yaxt="n", xlab="Effective density (De)", ylab="Dispersal spread (sigma)", ylim=ylims)
	polygon(x=c(Dlow, Dlow, min(Drange), min(Drange), Dlow, Dhigh, Dhigh), y=c(min(ylims), rousset(Dhigh,m), rousset(Dhigh,m), rousset(Dlow,m), rousset(Dlow,m), rousset(Dhigh,m), min(ylims)), col="grey", border=NA) # polygon for Dhigh to Dlow range
	polygon(c(Drange,rev(Drange)), c(sigmarangelow, rev(sigmarangehigh)),col="dark grey", border=NA) #plot a range for the sigma line from +/- mse
	lines(x=c(Dmlnelow, Dmlnelow, min(Drange), min(Drange)), y=c(min(ylims), rousset(Dmlnelow,m), rousset(Dmlnelow,m), rousset(Dmlnelow,m)), lty=2, lwd=1) # plot Dmlnelow
	lines(x=c(Dmlnehigh, Dmlnehigh, min(Drange), min(Drange)), y=c(min(ylims), rousset(Dmlnehigh,m), rousset(Dmlnehigh,m), rousset(Dmlnehigh,m)), lty=2, lwd=1) # plot Dmlnehigh
	lines(Drange, sigmarange, lwd=2) # replot the sigma line
	axis(1, labels=F, at=xticksm, tcl=-0.3)
	axis(1, labels=xticklg, at=xticklg)
	axis(2, labels=F, at=yticksm, tcl=-0.3)
	axis(2, labels=yticklg, at=yticklg)

	# Plot WILD21 and WILD21 separately on one graph
	plot(Drange, rousset(Drange, mWILD21), bty="l", log="xy", lwd=lwd1, type="l", xaxt="n", yaxt="n", xlab="Effective density (adults/km)", ylab="Dispersal spread (km)", ylim=ylims, col=col1) # line for WILD21 slope
	lines(Drange, rousset(Drange, mWILD21), lwd= lwd2, col=col2) # WILD21 slope
	lines(x=c(DmlnelowWILD21, DmlnelowWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DmlnelowWILD21,mWILD21), rousset(DmlnelowWILD21,mWILD21), rousset(DmlnelowWILD21,mWILD21)), lty=lty1, lwd= lwd1, col=col1) # plot DmlnelowWILD21
	lines(x=c(DmlnelowWILD21, DmlnelowWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DmlnelowWILD21,mWILD21), rousset(DmlnelowWILD21,mWILD21), rousset(DmlnelowWILD21,mWILD21)), lty=lty1, lwd= lwd2, col=col2) # plot DmlnelowWILD21
	lines(x=c(DmlnehighWILD21, DmlnehighWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DmlnehighWILD21,mWILD21), rousset(DmlnehighWILD21,mWILD21), rousset(DmlnehighWILD21,mWILD21)), lty=lty2, lwd= lwd1, col=col1) # plot DmlnehighWILD21
	lines(x=c(DmlnehighWILD21, DmlnehighWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DmlnehighWILD21,mWILD21), rousset(DmlnehighWILD21,mWILD21), rousset(DmlnehighWILD21,mWILD21)), lty=lty2, lwd= lwd2, col=col2) # plot DmlnehighWILD21
	lines(x=c(DreproWILD21, DreproWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DreproWILD21,mWILD21), rousset(DreproWILD21,mWILD21), rousset(DreproWILD21,mWILD21)), lty=lty3, lwd= lwd1, col=col1) # plot DreproWILD21
	lines(x=c(DreproWILD21, DreproWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DreproWILD21,mWILD21), rousset(DreproWILD21,mWILD21), rousset(DreproWILD21,mWILD21)), lty=lty3, lwd= lwd2, col=col2) # plot DmlnehighWILD21
	axis(1, labels=F, at=xticksm, tcl=-0.3)
	axis(1, labels=xticklg, at=xticklg)
	axis(2, labels=F, at=yticksm, tcl=-0.3)
	axis(2, labels=yticklg, at=yticklg)

	# Plot WILD21 and WILD21 separately on one graph for PPT
	lwd1 = 4
	lwd2 = 4

	par(mai = c(1,1,0.5,0.5))
	plot(Drange, rousset(Drange, mWILD21), bty="l", log="xy", lwd=lwd1, type="l", xaxt="n", yaxt="n", xlab="Effective density (De)", ylab="Dispersal spread (sigma)", ylim=ylims, col=col1, cex.lab=1.5) # line for WILD21 slope
	lines(Drange, rousset(Drange, mWILD21), lwd= lwd2, col=col2) # WILD21 slope
	axis(1, labels=F, at=xticksm, tcl=-0.3)
	axis(1, labels=xticklg, at=xticklg)
	axis(2, labels=F, at=yticksm, tcl=-0.3)
	axis(2, labels=yticklg, at=yticklg)


	lines(x=c(DmlnelowWILD21, DmlnelowWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DmlnelowWILD21,mWILD21), rousset(DmlnelowWILD21,mWILD21), rousset(DmlnelowWILD21,mWILD21)), lty=lty1, lwd= lwd1, col=col1) # plot DmlnelowWILD21
	lines(x=c(DmlnelowWILD21, DmlnelowWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DmlnelowWILD21,mWILD21), rousset(DmlnelowWILD21,mWILD21), rousset(DmlnelowWILD21,mWILD21)), lty=lty1, lwd= lwd2, col=col2) # plot DmlnelowWILD21
	lines(x=c(DmlnehighWILD21, DmlnehighWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DmlnehighWILD21,mWILD21), rousset(DmlnehighWILD21,mWILD21), rousset(DmlnehighWILD21,mWILD21)), lty=lty2, lwd= lwd1, col=col1) # plot DmlnehighWILD21
	lines(x=c(DmlnehighWILD21, DmlnehighWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DmlnehighWILD21,mWILD21), rousset(DmlnehighWILD21,mWILD21), rousset(DmlnehighWILD21,mWILD21)), lty=lty2, lwd= lwd2, col=col2) # plot DmlnehighWILD21
	lines(x=c(DreproWILD21, DreproWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DreproWILD21,mWILD21), rousset(DreproWILD21,mWILD21), rousset(DreproWILD21,mWILD21)), lty=lty3, lwd= lwd1, col=col1) # plot DreproWILD21
	lines(x=c(DreproWILD21, DreproWILD21, min(Drange), min(Drange)), y=c(min(ylims), rousset(DreproWILD21,mWILD21), rousset(DreproWILD21,mWILD21), rousset(DreproWILD21,mWILD21)), lty=lty3, lwd= lwd2, col=col2) # plot DmlnehighWILD21



########### Compare mean and standard deviation for various distributions 
# Lockwood defined mean as mean of half the distribution (e.g., the positive half)

# sd = 10

# normal
u = rnorm(10000, 0, 7)
sd = sd(u)
mean(u[u>0])
mean(u[u>0])/sd # 80%

# lognormal
u = rlnorm(10000, 0, 1.526)
v = -rlnorm(10000, 0, 1.526)
sd(c(u,v))
mean(u[u>0])
mean(u[u>0])/sd(c(u,v)) # 33%

# laplace
sd = 10
u = runif(len, -1/2, 1/2)
u = -sd/sqrt(2)*sign(u)*log(1-2*abs(u))
#hist(u)
sd(u)
mean(u[u>0])
mean(u[u>0])/sd(u) # 70%

# Weibull
u = rweibull(10000, .408)
v = -rweibull(10000, .408)
sd(c(u,v))
mean(u[u>0])
mean(u[u>0])/sd(c(u,v)) # 32%

# Gamma
u = rgamma(10000, 9.6)
v = -rgamma(10000, 9.6)
sd(c(u,v))
mean(u[u>0])
mean(u[u>0])/sd(c(u,v)) # 95%
