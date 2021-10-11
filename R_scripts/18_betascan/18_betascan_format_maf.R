#####################################
########## format mafs file  ########
#####################################
format_mafs <- function(headname){
  dat <-read.table(paste0(headname,"_all_minq20_minmq30_CV30_masked.mafs"), header = T)
  # do formatting for each chromosome
  for(j in c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')){ 
    DT = dat[which(dat$chromo %in% j),]
    # filter out SNPs with maf < 0.05 or maf > 0.95 in each population
    chr_ = c()
    pos_ = c()
    g_a1_ = c()
    g_tot_ = c()
    #for(i in seq(10)){
    for(i in seq(length(DT[,1]))){
      chr = DT$chromo[i]
      pos = DT$position[i]
      g_a1 = round(DT$knownEM[i]*DT$nInd[i]*2)
      g_tot = DT$nInd[i]*2
      if (DT$knownEM[i] > 0.05 && DT$knownEM[i] < 0.95){
        pos_ = c(pos_, pos)
        g_a1_ = c(g_a1_, g_a1)
        g_tot_ = c(g_tot_, g_tot)
      }
    }
    output = data.frame(pos_, g_a1_, g_tot_)
    output = output[order(output$pos_),]
    print(length(output$pos_))
    write.table(output, file = paste0(headname, "_", j ,".txt"), sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
  }
}

format_mafs("REF19")
format_mafs("CHR19")
format_mafs("HC")
format_mafs("NB")




#####################################
########## plot betascores   ########
#####################################

name1 = 'highlight.txt' 
h1 = read.delim(name1, header = FALSE, sep=',')
h1 = as.list(h1)

source("manhattan.R")
setwd("D:/Dropbox/cornell/DelBay19/18_betascan/test")
format_manhattan <- function(headname, chromosome){
  dat <-read.table(headname, header = T)
  colnames(dat)=c('position', 'Beta1.')
  dat$chromo <- chromosome
  dat$SNP <- paste0(dat$chromo, "_", dat$position)
  colnames(dat)=c('position', 'Beta1.', 'chromo', 'SNP' )
  manhattan(chr="chromo",bp="position",p="Beta1.", snp = "SNP", subset(dat, chromo == 5),  xlim = c(16000000, 17000000), highlight1 = h1$V1, logp=FALSE, cex.axis = 1.2, ylim = c(0, max(dat$Beta1.)),
  #manhattan(chr="chromo",bp="position",p="Beta1.", snp = "SNP", dat, highlight1 = h1$V1, logp=FALSE, cex.axis = 1.2, ylim = c(0, max(dat$Beta1.)),
            col=c("grey50","black"),genomewideline=F, suggestiveline=F,
            ylab=expression(beta~score), cex.lab=1.5)   
}

format_manhattan('NB_NC_035784.1.maf005.m005.out', 5)
format_manhattan('NB_NC_035784.1.maf005.m015.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.nofold.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w100.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w200.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w500.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w1000.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w2000.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w5000.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w10000.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w50000.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w100000.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w200000.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w100000.p20.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w100000.p4.out',5)
format_manhattan('NB_NC_035784.1.maf005.m015.w10000.p20.out',5)

# format_mafs <- function(headname){
#   dat <-read.table(paste0(headname,"_all_minq20_minmq30_CV30_masked.mafs"), header = T)
#   for(i in c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')){ 
#     DT = dat[which(dat$chromo %in% i),]
#     allele_total = DT$nInd*2
#     output = cbind(DT$position, round(DT$knownEM*allele_total), allele_total)
#     write.table(output, file = paste0(headname, "_",i,".txt"), sep = "\t", quote = FALSE,
#                 row.names = FALSE, col.names = FALSE)
#   }
# }





library(httr)
library(jsonlite)
library(xml2)
library(tibble)

getGoDetails <- function(go_id) {
  server <- "https://rest.ensembl.org/ontology/id/"
  
  r <- GET(paste(server, go_id, sep = ""), 
           content_type("application/json"))
  
  stop_for_status(r)
  
  res <- fromJSON(toJSON(content(r)))
  tibble(id = res$accession, 
         name = res$name, 
         definition = res$definition)
}
getGoDetails("GO:0005622")
