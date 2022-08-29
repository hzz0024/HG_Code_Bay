library(SoyNAM)
setwd("~/Dropbox/HG/Deep_learning_for_domestication/SoyNAM")
pop1 <- BLUP(trait="yield",family="all",env=3:4,dereg=FALSE,
     MAF=0.05,use.check=TRUE,impute="FM",rm.rep=TRUE)

pop2 <- BLUP(trait="yield",family="all",env="all",dereg=FALSE,
             MAF=0.05,use.check=TRUE,impute="FM",rm.rep=TRUE)

pop$Gen[,1]

test1 <- as.data.frame(pop$Gen)
test1$Year <- data.line$year[which(rownames(test1) %in% data.line$strain)]
test1$Pop <- pop$Fam
test1$Yield <- pop$Phen
test1[1:2, 4401:4403]
write.table(test1, "yield_n5487_snp4401.csv", row.names=F, col.names = T, quote=F, sep=",")

test2 <- data.line[which(data.line$strain) %in% rownames(test1),]
factor(test2$year)
factor(test2$location)
factor(test2$lodging)
factor(test2$size)
# ############ gen in #############
# The dataset including yield components collected at Purdue University (West Lafayette, Indiana)
# was used to investigate genomic prediction (Xavier et al. 2016) and interaction among traits
# (Xavier et al. 2017). This dataset contains the genotypic information in the matrix "gen.in",
# with missing values imputed using the software MaCH (Li et al. 2010). Similar to the datasets
# previously described, phenotypes are allocated into two objects, lines ("data.line.in") and checks
# ("data.checks.in"). Information on these data objects include year, location, environment (combination of year and location), strain, family, set (set in each environment), spot (combination of set                                                                                                                                                           and environment), the spatial coordinates of the field plots (BLOCK, ROW and COLUMN), plant
# height (in centimeters), R1 (number of days to flowering), R8 (number of days to maturity), lodging
# (score from 1 to 5), yield (in bu/ac), leaf shape (ratio length:width), number of nodes in the main
# stem, number of pods in the main stem, number of pods per node, average canopy coverage, rate
# of canopy coverage, growing degree day to flowering (GDD_R1), growing degree day to maturity
# (GDD_R1), and length of reproductive period in terms of growing degree day (GDD_REP).

raw_df <- as.data.frame(gen.in)
rownames(raw_df)
df_filter <- raw_df[which(rownames(raw_df) %in% data.line$strain),]
df_filter$row.names <- rownames(df_filter)
df_sort <- df_filter[order(df_filter$row.names),]

trait_df <- data.line[which(data.line$strain %in% df_sort$row.names),]
trait_df_sort <- trait_df[order(trait_df$strain),]

sum(trait_df_sort$strain == df_sort$row.names)

test5 <- data.line[which(rownames(test3) %in% data.line$strain),]
gen.in[1:2,1:2]
dim(gen.in)

data(soynam) 
geno = gen.raw[grep('DS1',rownames(gen.raw)),]
geno[,1]
fam=rownames(geno)
fam=gsub('DS1.-','',fam)
fam=gsub('...$','',fam,perl = T)
fam=as.numeric(fam)

# CHR
chr=rep(NA,20)
for(i in 1:20) chr[i]=length(grep(
  paste("Gm",sprintf("%02d",i),sep=''),
  colnames(geno)));rm(i)

trait = "yield"
test=dcast(data.check,environ+spot~strain,value.var=trait,mean)

# in data.check 46 strain in total, then rowmean across the strains, 
# in data.check 1753 environ + spot combination in total
# environ IA_2012 IA_2013 IL_2011 IL_2012 IL_2013 IN_2012 IN_2013 KS_2012 KS_2013 MI_2012 MO_2012 MO_2013 NE_2011 NE_2012 OHmc_2012 OHmc_2013 OHmi_2012 OHmi_2013  (18)
# spot (1753) 




