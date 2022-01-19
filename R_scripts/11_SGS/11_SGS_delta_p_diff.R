setwd("~/Dropbox/Mac/Documents/HG/DelBay20_adult/11_SGS")
pname = "ps_Del20_challenge.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4

pname = "ps_Del19_challenge.txt"
dat1 = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat1$V5[which(dat1$V5>=0)]))
dat1$delta_p <- dat1$V3-dat1$V4

data <- data.frame(
  type = c( rep("2020", length(dat$delta_p)), rep("2019",  length(dat1$delta_p))),
  Delta_P = c(dat$delta_p, dat1$delta_p))

p <- data %>%
  ggplot( aes(x=Delta_P, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  labs(fill="")
p


ks.test(dat1$delta_p, dat$delta_p)
