# format env file for baypass analysis
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/01_Salinity_index")
info_pop<-read.table("./info_pop_env.txt", header =T)
head(info_pop)
#env vector and its transpose form
env <- t(info_pop[,5:10])
write.table (env, "./env.txt", sep="\t", quote=F, row.names=F, col.names=F)
