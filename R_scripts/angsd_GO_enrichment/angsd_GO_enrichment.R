setwd("~/Downloads/physalia_adaptation_course-master/05_day5")

install.packages("splitstackshape")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("goseq")

library(goseq)
library(splitstackshape)
library(dplyr)

transcript_info = read.table("06_go/transcript_go_simplified.txt", header=T, stringsAsFactors=F, sep="\t")
head(transcript_info)
row.names(transcript_info)<-transcript_info$TranscriptName

#split the GO terms
go_split = cSplit(transcript_info, "GeneGo", sep = ";")
go_split$contig<- transcript_info$TranscriptName
head(go_split)

library(data.table)
#linearize the matrix
terms = colnames(select(go_split, contains("GeneGo")))
go_long = melt(go_split,measure.vars=terms, id.vars = "contig", na.rm=TRUE)
head(go_long)
go_ready = as.data.frame(go_long[,c(1, 3)])
head(go_ready)

#upload transcript intersecting with snps
all_transcripts<-read.table("06_go/all_snps.transcript", header=F)
colnames(all_transcripts)[1]<-"TranscriptName"
dim(all_transcripts) #how many?
head(all_transcripts)

library(dplyr)
#add size info
all_transcripts<-left_join(all_transcripts,transcript_info[,c(1,3)])
dim(all_transcripts)
head(all_transcripts)

#make unique
all_transcripts_unique<- all_transcripts %>% distinct(TranscriptName,.keep_all = TRUE)
dim(all_transcripts_unique) # see the matrix has reduced...
#and we need to name the rows b transcript names for goseq...
row.names(all_transcripts_unique)<-all_transcripts_unique$TranscriptName
head(all_transcripts_unique)

#transcripts in outliers
outliers_transcripts_temp_rda<-read.table("06_go/outlier_temp_rda.transcript", header=F)
colnames(outliers_transcripts_temp_rda)<-"TranscriptName"
head(outliers_transcripts_temp_rda)
dim(outliers_transcripts_temp_rda)

all_transcripts_unique$outliers_temp_rda<- as.numeric(all_transcripts_unique$TranscriptName %in% outliers_transcripts_temp_rda$TranscriptName)
head(all_transcripts_unique)

measured_genes = as.vector(all_transcripts_unique$TranscriptName)
outliers_genes = as.vector(all_transcripts_unique$outliers_temp_rda)
length = as.vector(all_transcripts_unique$length)
pwf_outliers = nullp(outliers_genes, bias.data=length)
row.names(pwf_outliers)<-row.names(all_transcripts_unique)
head(pwf_outliers)

enrich_outliers = goseq(pwf_outliers,gene2cat=go_ready, use_genes_without_cat=TRUE)
head(enrich_outliers)

enrich_outliers$over_represented_padjust<-p.adjust(enrich_outliers$over_represented_pvalue,method="BH")
head(enrich_outliers)
write.table(enrich_outliers, "06_go/GO_enrich_temp_RDA.txt", sep="\t")

enrich_outliers[which(enrich_outliers$over_represented_padjust<0.10),]



