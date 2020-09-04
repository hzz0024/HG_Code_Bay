a <- rbind(c(29,26),c(7,10))
fisher.test(a)
fisher.test(a, alternative = "greater")$p.value*2
###############################
file1 = 'REF_maf0.05_pctind0.7_cv30_anc.mafs'
file2 = 'REF_maf0.05_pctind0.7_cv30.mafs'

dat1 = read.delim(file1, header = TRUE, sep='\t')
dat2 = read.delim(file2, header = TRUE, sep='\t')

a = sum(dat1$chromo == dat2$chromo) == length(dat1$chromo)
b = sum(dat1$position == dat2$position) == length(dat1$position)
c = sum(dat1$major == dat2$major) == length(dat1$major)
d = sum(dat1$minor == dat2$minor) == length(dat1$minor)
e = sum(dat1$anc == dat2$anc) == length(dat1$anc)


index = dat1$major != dat2$major | dat1$minor != dat2$minor
diff_list = paste(dat1$chromo[index], dat1$position[index], dat1$major[index], dat1$minor[index], dat2$major[index], dat2$minor[index], sep = ' ')

datt1 = dat1[index,]
datt2 = dat2[index,]
index =  ! (datt1$major == datt2$minor & datt1$minor == datt2$major)
diff_list = paste(datt1$chromo[index], datt1$position[index], datt1$major[index], datt1$minor[index], datt2$major[index], datt2$minor[index], sep = ' ')





