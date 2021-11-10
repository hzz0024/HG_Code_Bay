library(gtools)
library(reshape2)
library(LDheatmap)
library(ggplot2)

pname1 = "./outlier/HC.chr5.ngsld.output"
dat1 = read.delim(pname1, header = FALSE, sep='\t')
my.results = dat1[with(dat1, order(V1, V2)),]
my.results <- my.results[complete.cases(my.results), ]
colnames(my.results)=c('Snp1', 'Snp2', 'Distance', "A", "B", "C", "Rsquared")
### Plot the Decay of LD with Distance
plot(my.results$Distance, y=my.results$Rsquared, xlab="Distance between SNPs (bp)", ylab=expression(R^2), pch=20, col=rgb(0,0,0,alpha=0.2), main="HC SGS Outlier Decay of Linkage Disequilibrium on Chr5")

### Add a smoothed moving average line
bins=seq(1,max(my.results$Distance), by=100) # Get all numbers between 1 and the maximum Distance value, in increments of 500
my.means=rep(0, length(bins)) # Set up an empty vector to hold the mean r-squared for each "bin"
LD.averages=data.frame(bins, my.means) # Create the results table (to hold calculations)

### Go through all of the "bin" values,
### For each one, find the subset of data that
### falls in that bin, and get the mean r-squared
for (i in 1:length(bins)) {
  data.interval=subset(my.results, (my.results$Distance >= bins[i] & my.results$Distance < (bins[i]+100))) # Get the subset of data that falls in the ith bin
  LD.averages$my.means[i]=mean(data.interval$Rsquared) # Calculate the average R-squared value for this data subset, and save in my "LD.averages" results table
}
LD.averages

### Add points and a line for my moving average on top of my exisiting plot
points(x=LD.averages$bins, y=LD.averages$my.means, col="red", pch=20) 

### Find the point when LD drops below 0.1
### The which statement tells me which rows of the table have means
### greater than 0.1,
### the brackets are giving me the subset of the "bins" with averages
### that correspond to these rows, and the max statement tells me which bin has
### the greatest distance value
LD.drop = max(LD.averages$bins[which(LD.averages$my.means>0.1)])
abline(v=LD.drop, col="blue", lty=2) # Add a vertical line corresponding to the drop off point
text(900, 0.8, "LD drops below 0.1", col="blue",cex=2) # label that line

### Find the LD half life
LD.half = (max(LD.averages$my.means))/2 # Calculate half of the maximum average LD value
LD.half.point = max(LD.averages$bins[which(LD.averages$my.means>LD.half)]) # Same strategy as when I found the 0.1 drop off point
abline(v=LD.half.point, col="darkgreen", lty=2) # add a vertical line at the half life point
text(LD.half.point, 0.9, "LD half life", col="darkgreen", cex=2) # label this line

### Calculate rho based on Hill Weir equation:
n=50 # get the number of samples by finding the number of genotypes
hill.weir.eq = (Rsquared~(((10+(rho*Distance))/((2+(rho*Distance))*(11+(rho*Distance))))*(1+(((3+(rho*Distance))*(12+(12*(rho*Distance))+((rho*Distance)**2)))/(n*(2+(rho*Distance))*(11+(rho*Distance))))))) # Define the formula for the Hill-Weir equation

rho.start=0.1 # Set an arbitrary starting value for rho
m=nls(formula=hill.weir.eq, data=my.results, start=list(rho=rho.start)) # Run the function to test different rho values and find the best fit to the data
results.m=summary(m) # Get a summary of results from the regression model

### Plot Expected R-squared
rho.estimate=results.m$parameters[1] # extract the rho estimate that produced the best fit
Distance=sort(my.results$Distance) # sort the distance results from lowest to highest

exp.rsquared=(((10+(rho.estimate*Distance))/((2+(rho.estimate*Distance))*(11+(rho.estimate*Distance))))*(1+(((3+(rho.estimate*Distance))*(12+(12*(rho.estimate*Distance))+((rho.estimate*Distance)**2)))/(n*(2+(rho.estimate*Distance))*(11+(rho.estimate*Distance)))))) # Use the best rho estimate inside of the Hill-Weir formula to calculate the expected r-squared at each distance value

lines(Distance, exp.rsquared, col="purple", lwd=2) # plot a line for the expected values (the same as a best fit line)
legend(600,0.94, c("Means", "Expected R2"), lty=c(1,1), col=c("red", "purple")) # add a legend


####################################
##########  plot r2 values  ########
####################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/outliter_block_plot/500bp/NB_HC/")
format <- function (pname1, pname2, i){
  #pname1 = "./HC.chr1.ngsld.output"
  dat1 = read.delim(pname1, header = FALSE, sep='\t')
  dat1 = dat1[with(dat1, order(V1, V2)),]
  #pname2 = "./NB.chr1.ngsld.output"
  dat2 = read.delim(pname2, header = FALSE, sep='\t')
  dat2 = dat2[with(dat2, order(V1, V2)),]
  df = data.frame(paste0("Chr", i), dat1$V7, dat2$V7)
  colnames(df)=c('Chr', 'HC', 'NB')
  df1 <- df[!is.na(df$HC), ]
  df2 <- df1[!is.na(df1$NB), ]
  return(df2)
}

chr1 <- format("HC.chr1.ngsld.output", "NB.chr1.ngsld.output", 1)
chr2 <- format("HC.chr2.ngsld.output", "NB.chr2.ngsld.output", 2)
chr3 <- format("HC.chr3.ngsld.output", "NB.chr3.ngsld.output", 3)
chr4 <- format("HC.chr4.ngsld.output", "NB.chr4.ngsld.output", 4)
chr5 <- format("HC.chr5.ngsld.output", "NB.chr5.ngsld.output", 5)
chr6 <- format("HC.chr6.ngsld.output", "NB.chr6.ngsld.output", 6)
chr7 <- format("HC.chr7.ngsld.output", "NB.chr7.ngsld.output", 7)
chr8 <- format("HC.chr8.ngsld.output", "NB.chr8.ngsld.output", 8)
chr9 <- format("HC.chr9.ngsld.output", "NB.chr9.ngsld.output", 9)
chr10 <- format("HC.chr10.ngsld.output", "NB.chr10.ngsld.output", 10)

df <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)

df_m <- melt(df, id.vars=c("Chr"))
df_f <- df_m[!is.na(df_m$value), ]
#then plot
means <- aggregate(value ~  variable, df_f, mean)

tiff("./HC_NB_LD_R2_500bp.jpg", units="in", width=16, height=12, res=300)
p2 <- ggplot(df_f, aes(x=factor(variable),y=value,fill=factor(variable)))+
  geom_boxplot(alpha = .7) + 
  labs(title="Wild HC-NB LD R2 comparison around SGS outlier candidates (window = 500bp)") +facet_wrap(~Chr) +
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1)) + 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6,vjust = -0.5)+
  guides(fill = "none") +
  theme_bw() +
  ylim(0,1) +
  labs( x = NULL, y = "R2")
p2
dev.off()

t.test(chr1$HC, chr1$NB)
t.test(chr2$HC, chr2$NB)
t.test(chr3$HC, chr3$NB)
t.test(chr4$HC, chr4$NB)
t.test(chr5$HC, chr5$NB)
t.test(chr6$HC, chr6$NB)
t.test(chr7$HC, chr7$NB)
t.test(chr8$HC, chr8$NB)
t.test(chr9$HC, chr9$NB)
t.test(chr10$HC, chr10$NB)

############## 2019 challenge ###########
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/outliter_block_plot/500bp/CHR19_REF19//")
format <- function (pname1, pname2, i){
  #pname1 = "./CHR19.chr1.ngsld.output"
  dat1 = read.delim(pname1, header = FALSE, sep='\t')
  dat1 = dat1[with(dat1, order(V1, V2)),]
  #pname2 = "./REF19.chr1.ngsld.output"
  dat2 = read.delim(pname2, header = FALSE, sep='\t')
  dat2 = dat2[with(dat2, order(V1, V2)),]
  df = data.frame(paste0("Chr", i), dat1$V7, dat2$V7)
  colnames(df)=c('Chr', 'CHR19', 'REF19')
  df1 <- df[!is.na(df$CHR19), ]
  df2 <- df1[!is.na(df1$REF19), ]
  return(df2)
}

chr1 <- format("CHR19.chr1.ngsld.output", "REF19.chr1.ngsld.output", 1)
chr2 <- format("CHR19.chr2.ngsld.output", "REF19.chr2.ngsld.output", 2)
chr3 <- format("CHR19.chr3.ngsld.output", "REF19.chr3.ngsld.output", 3)
chr4 <- format("CHR19.chr4.ngsld.output", "REF19.chr4.ngsld.output", 4)
chr5 <- format("CHR19.chr5.ngsld.output", "REF19.chr5.ngsld.output", 5)
chr6 <- format("CHR19.chr6.ngsld.output", "REF19.chr6.ngsld.output", 6)
chr7 <- format("CHR19.chr7.ngsld.output", "REF19.chr7.ngsld.output", 7)
chr8 <- format("CHR19.chr8.ngsld.output", "REF19.chr8.ngsld.output", 8)
chr9 <- format("CHR19.chr9.ngsld.output", "REF19.chr9.ngsld.output", 9)
chr10 <- format("CHR19.chr10.ngsld.output", "REF19.chr10.ngsld.output", 10)

df <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)

df_m <- melt(df, id.vars=c("Chr"))
df_f <- df_m[!is.na(df_m$value), ]
#then plot
means <- aggregate(value ~  variable, df_f, mean)

tiff("./CHR19_REF19_LD_R2_500bp.jpg", units="in", width=16, height=12, res=300)
p2 <- ggplot(df_f, aes(x=factor(variable),y=value,fill=factor(variable)))+
  geom_boxplot(alpha = .7) + 
  labs(title="2019 CHR19-REF19 LD R2 comparison around SGS outlier candidates (window = 500bp)") +facet_wrap(~Chr) +
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1)) + 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6,vjust = -0.5)+
  guides(fill = "none") +
  theme_bw() +
  ylim(0,1) +
  labs( x = NULL, y = "R2")
p2
dev.off()

t.test(chr1$CHR19, chr1$REF19)
t.test(chr2$CHR19, chr2$REF19)
t.test(chr3$CHR19, chr3$REF19)
t.test(chr4$CHR19, chr4$REF19)
t.test(chr5$CHR19, chr5$REF19)
t.test(chr6$CHR19, chr6$REF19)
t.test(chr7$CHR19, chr7$REF19)
t.test(chr8$CHR19, chr8$REF19)
t.test(chr9$CHR19, chr9$REF19)
t.test(chr10$CHR19, chr10$REF19)

############## SR-HC ###########
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/outliter_block_plot/500bp/SR_HC/")
format <- function (pname1, pname2, i){
  #pname1 = "./CHR19.chr1.ngsld.output"
  dat1 = read.delim(pname1, header = FALSE, sep='\t')
  dat1 = dat1[with(dat1, order(V1, V2)),]
  #pname2 = "./SR.chr1.ngsld.output"
  dat2 = read.delim(pname2, header = FALSE, sep='\t')
  dat2 = dat2[with(dat2, order(V1, V2)),]
  df = data.frame(paste0("Chr", i), dat1$V7, dat2$V7)
  colnames(df)=c('Chr', 'HC', 'SR')
  df1 <- df[!is.na(df$HC), ]
  df2 <- df1[!is.na(df1$SR), ]
  return(df2)
}

chr1 <- format("HC.chr1.ngsld.output", "SR.chr1.ngsld.output", 1)
chr2 <- format("HC.chr2.ngsld.output", "SR.chr2.ngsld.output", 2)
chr3 <- format("HC.chr3.ngsld.output", "SR.chr3.ngsld.output", 3)
chr4 <- format("HC.chr4.ngsld.output", "SR.chr4.ngsld.output", 4)
chr5 <- format("HC.chr5.ngsld.output", "SR.chr5.ngsld.output", 5)
chr6 <- format("HC.chr6.ngsld.output", "SR.chr6.ngsld.output", 6)
chr7 <- format("HC.chr7.ngsld.output", "SR.chr7.ngsld.output", 7)
chr8 <- format("HC.chr8.ngsld.output", "SR.chr8.ngsld.output", 8)
chr9 <- format("HC.chr9.ngsld.output", "SR.chr9.ngsld.output", 9)
chr10 <- format("HC.chr10.ngsld.output", "SR.chr10.ngsld.output", 10)

df <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)

df_m <- melt(df, id.vars=c("Chr"))
df_f <- df_m[!is.na(df_m$value), ]
#then plot
means <- aggregate(value ~  variable, df_f, mean)

tiff("./HC_SR_LD_R2_500bp.jpg", units="in", width=16, height=12, res=300)
p2 <- ggplot(df_f, aes(x=factor(variable),y=value,fill=factor(variable)))+
  geom_boxplot(alpha = .7) + 
  labs(title="Wild HC-SR LD R2 comparison around SGS outlier candidates (window = 500bp)") +facet_wrap(~Chr) +
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1)) + 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6,vjust = -0.5)+
  guides(fill = "none") +
  theme_bw() +
  ylim(0,1) +
  labs( x = NULL, y = "R2")
p2
dev.off()


t.test(chr1$HC, chr1$SR)
t.test(chr2$HC, chr2$SR)
t.test(chr3$HC, chr3$SR)
t.test(chr4$HC, chr4$SR)
t.test(chr5$HC, chr5$SR)
t.test(chr6$HC, chr6$SR)
t.test(chr7$HC, chr7$SR)
t.test(chr8$HC, chr8$SR)
t.test(chr9$HC, chr9$SR)
t.test(chr10$HC, chr10$SR)


####################################
##########  plot LD block   ########
####################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD")
TMP_FILE = "NB.ngsld.output"

df <- read.table(TMP_FILE, header=FALSE, stringsAsFactors=FALSE)
colnames(df)=c('snp1', 'snp2', 'dis', 'r1', 'D1', 'D2', 'r')
r <- df[complete.cases(df), ]
id <- unique(mixedsort(c(r[,"snp1"],r[,"snp2"])))
posStart <- head(id,1)
posEnd <- tail(id,1)
r <- rbind(r, c(posStart,posStart,0,NA,NA,NA,NA), c(posEnd,posEnd,0,NA,NA,NA,NA))

SNPs = read.delim("chr5_ADF1.txt", header = FALSE, sep='\t')

for (ld in c("r")) {
  m <- apply(acast(r, snp1 ~ snp2, value.var=ld, drop=FALSE),2,as.numeric)
  rownames(m) <- colnames(m)
  m <- m[mixedorder(rownames(m)),mixedorder(colnames(m))]
  id <- rownames(m)
  dist <- as.numeric(sub(".*:","",id))

  # Save plot
  tiff(paste("LD_blocks", ld,"jpg", sep="."), units="in", width=6, height=6, res=300)
  #pdf(paste("LD_blocks", ld,"pdf", sep="."), width=10, height=10)
  LDheatmap(m, genetic.distances=dist, geneMapLabelX=0.75, geneMapLabelY=0.25, title = NULL, color="blueToRed", LDmeasure=ld, SNP.name = SNPs$V1)
  require(grid)
  grid.edit("symbols", pch = 20, gp = gpar(cex = 1, col = "red")) # change the symbol size
  grid.edit(gPath("ldheatmap", "geneMap","SNPnames"), gp = gpar(cex=0.2, col = "red"))  # change the SNP label size
  dev.off()
}
