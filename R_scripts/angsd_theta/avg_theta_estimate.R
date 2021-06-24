file1 = 'Del20_challenge.thetas.window.idx.pestPG' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
theta1 = dat1$tW/dat1$nSites
avg1 = mean(theta1, na.rm=TRUE)
avg1

file1 = 'WILD_all.thetas.window.idx.pestPG' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
theta1 = dat1$tW/dat1$nSites
avg1 = mean(theta1, na.rm=TRUE)
avg1
