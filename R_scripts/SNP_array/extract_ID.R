make_ID <- function(file_name, global_name){
  dat = read.delim(file_name, header = FALSE, sep='\t')
  idx = dat$V3 < 0.05
  message(paste0('0.1: ',length(dat$V3[idx])))
  df = paste0(dat$V1[idx],'\t',dat$V2[idx],'\t',dat$V3[idx])
  #write.table(df, paste0(global_name, 'out_0.1.txt'), quote = FALSE, row.names = FALSE)
  dat_ID = paste0(dat$V1[idx],'_',dat$V2[idx])
  return(dat_ID)
}

z_a = make_ID('REF-CH-SR-HC_out_all_z.txt','REF-CH-SR-HC_z_')
z_b = make_ID('REF-CH-NB-HC_out_all_z.txt','REF-CH-NB-HC_z_')
z_c = make_ID('SR-REF-COH-ARN_out_all_z.txt','SR-REF-COH-ARN-z_')

z_d <- intersect(z_a, z_b)
length(z_d)
z_e <- intersect(z_a, z_c)
length(z_e)
z_f <- intersect(z_b, z_c)
length(z_f)

f_d <- intersect(f_a, f_b)
length(f_d)
f_e <- intersect(f_a, f_c)
length(f_e)
f_f <- intersect(f_b, f_c)
length(f_f)

m_a <- intersect(z_a, f_a) 
length(m_a)
m_b <- intersect(z_b, f_b) 
length(m_b)


f_a = make_ID('REF-CH-SR-HC_out_all_fish.txt','REF-CH-SR-HC_fish_')
f_b = make_ID('REF-CH-NB-HC_out_all_fish.txt','REF-CH-NB-HC_fish_')
f_c = make_ID('SR-REF-COH-ARN_out_all_fish.txt','SR-REF-COH-ARN_fish_')

setwd("~/Documents/Ryan_workplace/DelBay19_fisher/maf02")
z_a1 = make_ID('REF-CH-NB-HC_out_all_z.txt','REF-CH-SR-HC_z_')
setwd("~/Documents/Ryan_workplace/DelBay19_fisher/maf005")
z_a2 = make_ID('REF-CH-NB-HC_out_all_z.txt','REF-CH-SR-HC_z_')
z1 <- intersect(z_a1, z_a2)
length(z1)

setwd("~/Documents/Ryan_workplace/DelBay19_fisher/maf02")
f_a1 = make_ID('REF-CH-NB-HC_out_all_fish.txt','REF-CH-SR-HC_fjsh_')
setwd("~/Documents/Ryan_workplace/DelBay19_fisher/maf005")
f_a2 = make_ID('REF-CH-NB-HC_out_all_fish.txt','REF-CH-SR-HC_fish_')
f1 <- intersect(f_a1, f_a2)
length(z1)


extract_ID <- function(file, target, output_name){
  dat = read.delim(file, header = TRUE, sep='\t')
  names = paste0(dat$chromo,'_',dat$position)
  datt = dat[names %in% target,]
  df = paste0(datt$chromo,'\t',datt$position,'\t','[',datt$major,'/',datt$minor,']','\t',datt$anc)
  write.table(df, paste0(output_name, 'fdr_0.05.txt'), quote = FALSE, row.names = FALSE)
}

extract_ID('REF_maf0.05_pctind0.7_cv30.mafs',f_a, 'REF-CH-SR-HC_fish_')
extract_ID('REF_maf0.05_pctind0.7_cv30.mafs',f_b, 'REF-CH-NB-HC_fish_')
extract_ID('REF_maf0.05_pctind0.7_cv30.mafs',f_c, 'SR-REF-COH-ARN_fish_')
extract_ID('REF_maf0.05_pctind0.7_cv30.mafs',z_d, 'REF-CH-SR-HC_REF-CH-NB-HC_z_fdr0.1_')

reseq_id <- c("1_55485508","2_2405242","2_59641818","5_49795690","7_53403","8_36715103","9_50117923","2_656131","3_54596543","8_63266137","9_12057031","3_69917769","5_25337621","5_79382471","8_10601236","9_52331392")
reseq_id <- sort(reseq_id)
extract_ID('REF_maf0.05_pctind0.7_cv30.mafs',reseq_id, 'reseq_dom_vs_wild_')
