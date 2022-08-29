setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/02_baypass/Simulation")

#source functions from Gautier file
source("./baypass_utils.R")

#load omega matrix
omega = as.matrix(read.table("prunedsnps.output_mat_omega.out"))

#load beta parameters
pi.beta.coef=read.table("prunedsnps.output_summary_beta_params.out", h=T)$Mean

#make geno matrix
#lobster.data <- geno2YN(args[3])

#create the POD - simulate - don't forget to change suffix name
simu.bta <- simulate.baypass(omega.mat=omega, nsnp=10000,beta.pi=pi.beta.coef, suffix="simulates_pods")
