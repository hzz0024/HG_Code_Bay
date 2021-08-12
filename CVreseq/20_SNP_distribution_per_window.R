library(ggplot2)
headname = "CS_HC_noinvers."
titlename = "CS_HC"
snp_density <- function(headname, titlename){
  for(win in c(1000)){
    name = paste0(headname, win, "bp.", "s", win/5, ".csv")
    DT = read.delim(name, header = TRUE, sep=',')
    mid_pos <- round((DT$start + DT$end)/2)
    id = paste0(DT$scaffold,'_',mid_pos)
    DT <- as.data.frame(cbind(DT,mid_pos, id))
    #DT <- DT[complete.cases(DT), ]
    #DT[,9][DT[,9]<0] = 0.000001 #@chnage
    #zfst <- (DT[,9] - mean(DT[,9]))/(sd(DT[,9], na.rm = FALSE))# @change
    dat <- data.frame(chr=DT$scaffold, start=DT$start, end=DT$end, mid_pos=DT$mid_pos, SNP=DT$id, sites=DT$sites) # @change
    dat$chr <- as.numeric(dat$chr)
    dat$mid_pos <- as.numeric(dat$mid_pos)
    dat$sites <- as.numeric(dat$sites)
    jpeg(paste0(headname,"sliding.snpdensity.jpg"), width = 16, height = 9, units = 'in', res = 300)
    par(mfrow=c(1,1))
    p<-ggplot(dat, aes(x=sites)) + 
      geom_histogram(color="black", fill="white")
    p
    dev.off()
    mid = mean(dat$sites)
    print(paste0("Mean sites is ", mid))
    med = median(dat$sites)
    print(paste0("median sites is ", med))
  }
}


############################## formal run ##############################
############################## formal run ##############################
############################## formal run ##############################

snp_density("SL_OBOYS2_noinvers.", "SL_OBOYS2")
snp_density("SL_LOLA_noinvers.", "SL_LOLA")
snp_density("NEH_UMFS_noinvers.", "NEH_UMFS")
snp_density("CS_UMFS_noinvers.", "CS_UMFS")
snp_density("CS_NEH_noinvers.", "CS_NEH")
snp_density("CS_DEBY_noinvers.", "CS_DEBY")
snp_density("CL_OBOYS2_noinvers.", "CL_OBOYS2")
snp_density("CS_SL_noinvers.", "CS_SL")

snp_density("CS_HC_noinvers.", "CS_HC")
snp_density("HC_CLP_noinvers.", "HC_CLP")
snp_density("HCVA_CLP_noinvers.", "HCVA_CLP")
snp_density("CS_HCVA_noinvers.", "CS_HCVA")