####################################
##########  Delta_p plot ###########
####################################

# customize the color
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  ## Save the color
  invisible(t.col)
}
mycol <- t_col("black", perc = 50, name = "lt.grey")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/power")
# load the dataset
fname = "starting_af_power.txt"  
df_maf = read.delim(fname, header = TRUE, sep='\t')
head(df_maf)
dim(df_maf)
df_maf$id = paste0(df_maf$chromo,'_',df_maf$position)
df_maf$order = seq(1:nrow(df_maf))
head(df_maf)

plot(df_maf$order, df_maf$p0,col="red",pch=20,cex=1,ylim=c(0,1),xlab="SNPs",ylab="Allele frequency",xlim=c(0,50)) # ylab=expression(italic("p"))
for (l in 1:nrow(df_maf))
{
  segments(l,df_maf$p1[l],l,df_maf$p0[l],col=mycol,lwd=2.0, lty = "dotted")
  #text(df_maf$order[l], df_maf$p1[l], labels=df_maf$delta_p[l],data=df_maf, cex=0.7, font=2, pos=3)
  text(df_maf$order[l], df_maf$p1[l], labels=df_maf$ps[l],data=df_maf, cex=0.7, font=2, pos=3)
  if (df_maf$delta_p[l]<=0)
  {
    points(l,df_maf$p1[l],col=mycol,pch=6,cex=1)
  }
  if (df_maf$delta_p[l]>0)
  {
    points(l,df_maf$p1[l],col=mycol,pch=17,cex=1)
  }
}

graph2ppt(file="deltap_test_plot.pptx", width=10, height=8)
