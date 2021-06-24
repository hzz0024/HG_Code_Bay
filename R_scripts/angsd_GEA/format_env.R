
ENV_FILE<-"env_pcaenv_16pop.txt"
VAR<-c(2,3,4,9)

env_all<-read.table(ENV_FILE, header=T, row.names = 1)
env<-env_all[,1]
env<-env_all[,VAR-1]

#check any remaining correlation
cormat<-cor(env)
p_cor<-cor.test.p (env)
col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                           "#D6604D", "#B2182B", "#67001F")) (200)
corrplot((cormat),   method="square",col=col2, tl.col=1, cl.pos="b", p.mat=p_cor, insig="label_sig", pch.cex=1)

#prepare files for env assoc
write.table(t(env), "env_baypass", quote=F, row.names=F, col.names=F)
write.table(env, "env_lfmm.env", quote=F, row.names=F, col.names=F)
write.table(env, "env_rda", quote=F)
