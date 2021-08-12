library("ggplot2")

plotclass <- function(headname, titlename){
  #jpeg(paste0(headname,".jpg"), width = 6, height = 6, units = 'in', res = 300)
  DT = read.delim(headname, header = FALSE, sep='\t')
  SNP_cnt <- length(DT$V1)
  print(paste0("number of SNPs is ", SNP_cnt))
  cnt <- data.frame(table(DT$V1))
  colnames(cnt) <- c("Class", "Count")
  write.table(cnt, file = paste0(headname,".csv"), sep = ",", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  # plot as a bar chart
  p <- ggplot(cnt, aes(x=Class, y=Count, fill=Class)) + geom_bar(stat="identity")
  p + ggtitle(titlename) + theme(plot.title=element_text(face="bold")) +
    xlab("") + ylab("Count")+ labs(x = "", fill = "Class") +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 
    #scale_fill_manual(values=c("#4c4c4c", "#86BB8D", "#68a4bd", "#ff9900"), name="Class")
  #dev.off()
}

headname = "Genome-wide.noinvers.variant_function"

jpeg("CS_HC-HCVA_CLP.outliers.variant_function.jpg", width = 6, height = 6, units = 'in', res = 300)
plotclass("CS_HC-HCVA_CLP.outliers.variant_function", "SNPs respond to salinity adaptation (n=463)") 
dev.off()
jpeg("Genome-wide.noinvers.chr2.variant_function.jpg", width = 6, height = 6, units = 'in', res = 300)
plotclass("Genome-wide.noinvers.chr2.variant_function","SNPs in chromosome 2 (outside inversions, n=668282)") 
dev.off()
jpeg("95.outlier.SNPs.inversion.variant_function.jpg", width = 6, height = 6, units = 'in', res = 300)
plotclass("95.outlier.SNPs.inversion.variant_function","Jon's outlier SNP list within the inversion (n=2269)")
dev.off()
jpeg("95.outlier.SNPs.no_inversion.variant_function.jpg", width = 6, height = 6, units = 'in', res = 300)
plotclass("95.outlier.SNPs.no_inversion.variant_function","Jon's outlier SNP list outside the inversion (n=4397)")
dev.off()
jpeg("95.outlier.SNPs.wildae.no_inversion.variant_function.jpg", width = 6, height = 6, units = 'in', res = 300)
plotclass("95.outlier.SNPs.wildae.no_inversion.variant_function","Jon's outlier SNP list from wild Atlantic populations (n=315)")
dev.off()
jpeg("Genome-wide.invers.variant_function.jpg", width = 6, height = 6, units = 'in', res = 300)
plotclass("Genome-wide.invers.variant_function","Genome-wide SNPs (within inversions, n=1125493)") 
dev.off()
jpeg("Genome-wide.noinvers.variant_function.jpg", width = 6, height = 6, units = 'in', res = 300)
plotclass("Genome-wide.noinvers.variant_function","Genome-wide SNPs (outside inversions, n=4640654)") 
dev.off()
################### Stacked barplots for SNPs within coding regions ###################
require(reshape)
require(ggplot2)


headname <- "CS_HC-HCVA_CLP.outliers.exonic_variant_function"
jpeg(paste0(headname,".jpg"), width = 6, height = 6, units = 'in', res = 300)
df = read.delim(headname, header = FALSE, sep='\t')
length(df$V1)
df$V12=unlist(lapply(strsplit(df$V3, ":"), `[[`, 1))
df$V13=unlist(lapply(strsplit(df$V2, " "), `[[`, 1))
tmp = data.frame(df$V13, df$V12)
names(tmp) = c('Class', 'Gene')
tmp[tmp$Gene=='gene6173',]$Gene = 'CDP-diacylglycerol--glycerol-3\n-phosphate 3-phosphatidyltransferase'
tmp[tmp$Gene=='gene6174',]$Gene = 'dynein beta chain, ciliary-like'
tmp = as.data.frame(table(tmp))
tmp$Gene <- factor(tmp$Gene, levels =c('CDP-diacylglycerol--glycerol-3\n-phosphate 3-phosphatidyltransferase','dynein beta chain, ciliary-like'))
# write.table(tmp, file = "test.txt", sep = "\t",
#             row.names = TRUE, col.names = NA, quote = FALSE)
p <- qplot(Gene, data=tmp, geom="bar", weight=Freq, fill=Class)
p + ggtitle("Exonic SNP loci respond to salinity adaptation (n=216)") +
    scale_fill_manual(values=c("#4c4c4c", "#86BB8D", "#68a4bd", "#ff9900"), name="Variant\nclass")+ 
    xlab("") + ylab("Count") +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
          axis.title=element_text(size=14),
          plot.title = element_text(size = 14, face = "bold"))
dev.off()

headname <- "95.outlier.SNPs.inversion.exonic_variant_function"
jpeg(paste0(headname,".jpg"), width = 8, height = 6, units = 'in', res = 300)
df = read.delim(headname, header = FALSE, sep='\t')
length(df$V1)
df$V12=unlist(lapply(strsplit(df$V3, ":"), `[[`, 1))
df$V13=unlist(lapply(strsplit(df$V2, " "), `[[`, 1))
tmp = data.frame(df$V13, df$V12)
names(tmp) = c('Class', 'Gene')
tmp[tmp$Gene=='gene10032',]$Gene = 'cell wall protein DAN4-like'
tmp[tmp$Gene=='gene10035',]$Gene = 'dynein heavy chain 8, axonemal-like'
tmp[tmp$Gene=='gene10044',]$Gene = 'uncharacterized LOC111123444'  
tmp[tmp$Gene=='gene10062',]$Gene = 'sodium bicarbonate cotransporter 3-like' 
tmp[tmp$Gene=='gene10063',]$Gene = 'uncharacterized LOC111122613' 
tmp[tmp$Gene=='gene10064',]$Gene = 'uncharacterized LOC111127523' 
tmp[tmp$Gene=='gene22309',]$Gene = 'cyclic nucleotide-gated channel \nrod photoreceptor subunit alpha-like ' 
tmp[tmp$Gene=='gene32280',]$Gene = 'uncharacterized LOC111109272' 
tmp = as.data.frame(table(tmp))
#tmp$Gene <- factor(tmp$Gene, levels =c('PGS1', 'DYHC'))
p <- qplot(Gene, data=tmp, geom="bar", weight=Freq, fill=Class)
p + ggtitle("Exonic SNP within the inversion (n=32)") +
  scale_fill_manual(values=c("#4c4c4c", "#86BB8D", "#68a4bd", "#ff9900"), name="Variant\nclass")+ 
  xlab("") + ylab("Count") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.title=element_text(size=14),
        plot.title = element_text(size = 14, face = "bold"))
dev.off()


headname <- "95.outlier.SNPs.no_inversion.exonic_variant_function"
jpeg(paste0(headname,".jpg"), width = 10, height = 6, units = 'in', res = 300)
df = read.delim(headname, header = FALSE, sep='\t')
length(df$V1)
df$V12=unlist(lapply(strsplit(df$V3, ":"), `[[`, 1))
df$V13=unlist(lapply(strsplit(df$V2, " "), `[[`, 1))
tmp = data.frame(df$V13, df$V12)
names(tmp) = c('Class', 'Gene')
tmp[tmp$Gene=='gene10326',]$Gene = 'T-box transcription factor TBX1-A-like'
tmp[tmp$Gene=='gene10329',]$Gene = 'low-density lipoprotein receptor-like'
tmp[tmp$Gene=='gene10332',]$Gene = 'T-box transcription factor mls-1-like'
tmp[tmp$Gene=='gene10333',]$Gene = 'T-box-containing protein TBX6L-like '
tmp[tmp$Gene=='gene10334',]$Gene = 'T-box transcription factor TBX20-like '
tmp[tmp$Gene=='gene10336',]$Gene = 'von Willebrand factor C domain-containing protein 2-like'
tmp[tmp$Gene=='gene10338',]$Gene = 'N-acetylneuraminate lyase B-like'
tmp[tmp$Gene=='gene10340',]$Gene = 'uncharacterized LOC111126902 '
tmp[tmp$Gene=='gene13891',]$Gene = 'proline-rich protein PRCC-like'
tmp[tmp$Gene=='gene13894',]$Gene = 'uncharacterized LOC111131263'
tmp[tmp$Gene=='gene13901',]$Gene = 'arginase, hepatic-like'
tmp[tmp$Gene=='gene17785',]$Gene = 'uncharacterized LOC111134359' 
tmp[tmp$Gene=='gene17788',]$Gene = 'atrial natriuretic peptide receptor 1-like'
tmp[tmp$Gene=='gene19770',]$Gene = 'phosphatidylinositol glycan anchor biosynthesis class U protein-like' 
tmp[tmp$Gene=='gene19771',]$Gene = 'tRNA-dihydrouridine(20) synthase [NAD(P)+]-like'
tmp[tmp$Gene=='gene19773',]$Gene = 'nuclear pore complex protein Nup93-like '
tmp[tmp$Gene=='gene19774',]$Gene = 'uncharacterized LOC111137708'
tmp[tmp$Gene=='gene19778',]$Gene = 'H(+)/Cl(-) exchange transporter 7-like'
tmp[tmp$Gene=='gene19779',]$Gene = 'guanine nucleotide-binding protein G(s) subunit alpha-like'
tmp[tmp$Gene=='gene19780',]$Gene = 'protein fantom-like'
tmp[tmp$Gene=='gene19781',]$Gene = 'ceramide kinase-like'
tmp[tmp$Gene=='gene19837',]$Gene = 'poly [ADP-ribose] polymerase 2-like'
tmp[tmp$Gene=='gene19839',]$Gene = 'poly [ADP-ribose] polymerase 2 isoform X1'
tmp[tmp$Gene=='gene22326',]$Gene = 'RRP12-like protein'
tmp[tmp$Gene=='gene22327',]$Gene = 'threonine synthase-like 2'
tmp[tmp$Gene=='gene22328',]$Gene = 'uncharacterized LOC111133424'
tmp[tmp$Gene=='gene25575',]$Gene = 'cell division cycle 5-like protein'
tmp[tmp$Gene=='gene25576',]$Gene = 'dystrophin-like'
tmp[tmp$Gene=='gene25577',]$Gene = 'uveal autoantigen with coiled-coil domains and ankyrin repeats-like '
tmp[tmp$Gene=='gene25592',]$Gene = 'nacrein-like protein'
tmp[tmp$Gene=='gene25593',]$Gene = 'uncharacterized LOC111104062'
tmp[tmp$Gene=='gene26206',]$Gene = 'ashwin-like '
tmp[tmp$Gene=='gene26207',]$Gene = 'snRNA-activating protein complex subunit 1-like'
tmp[tmp$Gene=='gene26208',]$Gene = 'uncharacterized LOC111104812'
tmp[tmp$Gene=='gene5427',]$Gene = 'adhesion G protein-coupled receptor L3-like '
tmp[tmp$Gene=='gene5428',]$Gene = 'ferric-chelate reductase 1-like'
tmp[tmp$Gene=='gene5430',]$Gene = 'uncharacterized LOC111120629'
tmp = as.data.frame(table(tmp))
#tmp$Gene <- factor(tmp$Gene, levels =c('PGS1', 'DYHC'))
p <- qplot(Gene, data=tmp, geom="bar", weight=Freq, fill=Class)
p + ggtitle("Exonic SNP outside the inversion (n=365)") +
  scale_fill_manual(values=c("#4c4c4c", "#86BB8D", "#68a4bd", "#ff9900"), name="Variant\nclass")+ 
  xlab("") + ylab("Count") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=6),
        axis.title=element_text(size=14),
        plot.title = element_text(size = 14, face = "bold")) +
  theme(plot.margin = unit(c(0.1, 0.1 ,0.1, 2), "cm"))
dev.off()


headname <- "95.outlier.SNPs.wildae.no_inversion.exonic_variant_function"
jpeg(paste0(headname,".jpg"), width = 6, height = 6, units = 'in', res = 300)
df = read.delim(headname, header = FALSE, sep='\t')
length(df$V1)
df$V12=unlist(lapply(strsplit(df$V3, ":"), `[[`, 1))
df$V13=unlist(lapply(strsplit(df$V2, " "), `[[`, 1))
tmp = data.frame(df$V13, df$V12)
names(tmp) = c('Class', 'Gene')
tmp[tmp$Gene=='gene133',]$Gene = 'protocadherin Fat 4-like'
tmp[tmp$Gene=='gene3231',]$Gene = 'elongation factor 1-alpha'
tmp[tmp$Gene=='gene9070',]$Gene = 'uncharacterized LOC111122899'  
tmp[tmp$Gene=='gene10750',]$Gene = 'barH-like 1 homeobox protein' 
tmp[tmp$Gene=='gene10751',]$Gene = 'barH-like 1 homeobox protein' 
tmp[tmp$Gene=='gene10752',]$Gene = 'sulfiredoxin-1-like' 
tmp[tmp$Gene=='gene10753',]$Gene = 'von Hippel-Lindau disease tumor suppressor-like' 
tmp[tmp$Gene=='gene10754',]$Gene = 'cap-specific mRNA (nucleoside-2-O-)\n-methyltransferase 2-like' 
tmp[tmp$Gene=='gene10755',]$Gene = 'uncharacterized LOC111127120' 
tmp[tmp$Gene=='gene10756',]$Gene = 'uncharacterized LOC111123760' 
tmp[tmp$Gene=='gene13457',]$Gene = 'patatin-like phospholipase domain-containing protein 7' 
tmp[tmp$Gene=='gene13458',]$Gene = 'maleylacetoacetate isomerase-like'
tmp[tmp$Gene=='gene13460',]$Gene = 'protein O-mannosyl-transferase 2-like' 
tmp[tmp$Gene=='gene13461',]$Gene = 'uncharacterized LOC111131424' 
tmp[tmp$Gene=='gene13903',]$Gene = 'thyrostimulin alpha-2 subunit-like' 
tmp[tmp$Gene=='gene15830',]$Gene = 'zinc finger protein-like 1' 
tmp = as.data.frame(table(tmp))
#tmp$Gene <- factor(tmp$Gene, levels =c('PGS1', 'DYHC'))
p <- qplot(Gene, data=tmp, geom="bar", weight=Freq, fill=Class)
p + ggtitle("Exonic SNP in wild Atlantic populations (n=47)") +
  scale_fill_manual(values=c("#4c4c4c", "#86BB8D", "#68a4bd", "#ff9900"), name="Variant\nclass")+ 
  xlab("") + ylab("Count") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=6),
        axis.title=element_text(size=14),
        plot.title = element_text(size = 14, face = "bold")) +
  theme(plot.margin = unit(c(0.1, 0.1 ,0.1, 2), "cm"))
dev.off()

################## counting syn vs non-syn ##################

sum(tmp$Class == "synonymous")
sum(tmp$Class == "nonsynonymous")
sum(tmp$Class == "stopgain")
sum(tmp$Class == "stoploss")
sum(tmp$Class == "unknown")

headname <- "CS_HC-HCVA_CLP.outliers.exonic_variant_function"
headname <- "95.outlier.SNPs.inversion.exonic_variant_function"
headname <- "95.outlier.SNPs.no_inversion.exonic_variant_function"
headname <- "95.outlier.SNPs.wildae.no_inversion.exonic_variant_function"
headname <- "Genome-wide.invers.exonic_variant_function"
headname <- "Genome-wide.noinvers.exonic_variant_function"
headname <- "Genome-wide.noinvers.chr2.exonic_variant_function"

plotratio <- function(headname, titlename){
#headname <- "CS_HC-HCVA_CLP_ratio.txt"
  df = read.delim(headname, header = TRUE, sep='\t')
  p <- qplot(Group, data=df, geom="bar", weight=Ratio, fill=Class)
  p + ggtitle(titlename) +
    scale_fill_manual(values=c("#4c4c4c", "#86BB8D", "#68a4bd", "#ff9900"), name="Variant\nclass")+ 
    xlab("") + ylab("Frequency") +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12),
          axis.title=element_text(size=14),
          plot.title = element_text(size = 14, face = "bold")) +
    theme(plot.margin = unit(c(0.1, 0.1 ,0.1, 2), "cm"))
}
jpeg(paste0("CS_HC-HCVA_CLP_ratio.jpg"), width = 6, height = 6, units = 'in', res = 300)
plotratio("CS_HC-HCVA_CLP_ratio.txt","Exonic genome-wide SNPs vs. \noutlier SNPs respond to salinity adaptation")
dev.off()
jpeg(paste0("95.outlier.SNPs.inversion.jpg"), width = 6, height = 6, units = 'in', res = 300)
plotratio("95.outlier.SNPs.inversion.txt","Exonic genome-wide SNPs within inversions\n vs outlier SNPs within inversions")
dev.off()
jpeg(paste0("95.outlier.SNPs.noinversion.jpg"), width = 6, height = 6, units = 'in', res = 300)
plotratio("95.outlier.SNPs.noinversion.txt","Exonic genome-wide SNPs outside inversions\n vs outlier SNPs outside inversions")
dev.off()
jpeg(paste0("95.outlier.SNPs.wildae.no_inversion.jpg"), width = 6, height = 6, units = 'in', res = 300)
plotratio("95.outlier.SNPs.wildae.no_inversion.txt","Exonic outlier SNPs in wild Atlantic contrasts\n vs genome-wide SNPs outside inversions")
dev.off()

############### chi-square test ###############

chisq_test <- function(tag_name,ref_name){
  #tag_name <- "CS_HC-HCVA_CLP.outliers.exonic_variant_function"
  df1 = read.delim(tag_name, header = FALSE, sep='\t')
  cnt1 = length(df1$V1)
  print(paste0("SNP counts in ", tag_name, " is ", cnt1))
  df1$V13="Target"
  df1$V14=unlist(lapply(strsplit(df1$V2, " "), `[[`, 1))
  tmp1 = data.frame(df1$V13, df1$V14)
  names(tmp1) = c('Data', 'Class')
  tmp1<- tmp1[which(tmp1$Class %in% "synonymous" | tmp1$Class %in% "nonsynonymous"),]
  
  #ref_name <- "Genome-wide.noinvers.exonic_variant_function"
  df2 = read.delim(ref_name, header = FALSE, sep='\t')
  cnt2 = length(df2$V1)
  print(paste0("SNP counts in ", ref_name, " is ", cnt2))
  df2$V13="Reference"
  df2$V14=unlist(lapply(strsplit(df2$V2, " "), `[[`, 1))
  tmp2 = data.frame(df2$V13, df2$V14)
  names(tmp2) = c('Data', 'Class')
  tmp2<- tmp2[which(tmp2$Class %in% "synonymous" | tmp2$Class %in% "nonsynonymous"),]
  
  dat <- rbind(tmp1,tmp2) 
  tab <- table(dat$Data,dat$Class)
  chi <- chisq.test(tab)
  print(chi$observed)
  print(chi$expected)
  print(paste0("p-value:",chi$p.value, "; chi-square:",chi$statistic, "; df:",chi$parameter))
  return(tab)
}
setwd("/Users/ryan/Documents/Ryan_workplace/CVreseq_annotation/annovar_results")
tab1<-chisq_test("95.outlier.SNPs.no_inversion.exonic_variant_function","Genome-wide.noinvers.exonic_variant_function")
chisq_test("95.outlier.SNPs.wildae.no_inversion.exonic_variant_function","Genome-wide.noinvers.exonic_variant_function")
chisq_test("CS_HC-HCVA_CLP.outliers.exonic_variant_function","Genome-wide.noinvers.exonic_variant_function")


