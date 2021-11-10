setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/Manuscript")
library(ggplot2)
library(dplyr)
library(hrbrthemes)

input = "Long_term_monitoring.txt"
input = "short_term_monitoring.txt"
dat = read.delim(input, header = TRUE, sep='\t')
dat$Month[dat$Month == 1] = '01'
dat$Month[dat$Month == 2] = '02'
dat$Month[dat$Month == 3] = '03'
dat$Month[dat$Month == 4] = '04'
dat$Month[dat$Month == 5] = '05'
dat$Month[dat$Month == 6] = '06'
dat$Month[dat$Month == 7] = '07'
dat$Month[dat$Month == 8] = '08'
dat$Month[dat$Month == 9] = '09'
dat$yearmon <- paste0(dat$Year, "_", dat$Month)
head(dat)
dim(dat)
new = levels(factor(dat$yearmon))
new_order = new[order(as.numeric(gsub("(\\d+)_(\\d+)", "\\1\\2", new)))]
dat$yearmon <- factor(dat$yearmon, levels = new_order)
cbPalette <- c("#4461A8", "#0BC0B3", "#EAC728", "#E97302", "#A71B4B")

# dat$month = factor(dat$Month, levels = month.abb)

ggplot(dat, aes(x=yearmon, y=RecentMortality, fill=Bed)) +
  geom_bar(aes(color=Bed), stat="identity", color="#A9A9A9", alpha = 0.9, position=position_dodge())+
  #geom_line(aes(group=Bed, color=Bed), alpha=0.2) +
  #geom_point(aes(color=Bed, group = Bed), alpha=0.8, cex=1) +
  #scale_shape_manual(rep(c(15),5), breaks=c("Hope Creek", "Arnolds", "Cohansey", "Shell Rock", "New Beds")) + 
  scale_colour_manual(values=cbPalette, breaks=c("Hope Creek", "Arnolds", "Cohansey", "Shell Rock", "New Beds"))+
  scale_fill_manual(values=cbPalette, breaks=c("Hope Creek", "Arnolds", "Cohansey", "Shell Rock", "New Beds"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

graph2ppt(file="Long-term.pptx", width=10, height=8) 
graph2ppt(file="short-term.pptx", width=10, height=8) 

ggplot(dat, aes(x=yearmon, y=Salinity, fill=Bed)) +
  geom_line(aes(group=Bed, color=Bed), alpha=0.5) +
  geom_point(aes(color=Bed, group = Bed), alpha=0.8, cex=0.1) +
  #scale_shape_manual(rep(c(15),5), breaks=c("Hope Creek", "Arnolds", "Cohansey", "Shell Rock", "New Beds")) + 
  scale_colour_manual(values=cbPalette, breaks=c("Hope Creek", "Arnolds", "Cohansey", "Shell Rock", "New Beds"))+
  scale_fill_manual(values=cbPalette, breaks=c("Hope Creek", "Arnolds", "Cohansey", "Shell Rock", "New Beds"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

graph2ppt(file="short-salinity.pptx", width=10, height=8) 
graph2ppt(file="long-salinity.pptx", width=10, height=8)
ggplot(dat, aes(x=Date, y=Salinity, fill=Bed)) + 
  geom_bar(stat="identity", position=position_dodge())+
  scale_colour_manual(values=cbPalette, breaks=c("Hope Creek", "Arnolds", "Cohansey", "Shell Rock", "New Beds"))+
  theme_classic()
