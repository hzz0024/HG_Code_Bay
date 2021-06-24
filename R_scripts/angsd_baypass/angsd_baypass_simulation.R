setwd("~/Documents/Ryan_workplace/DelBay_adult/13_env_gen_association/angsd_results")
args = commandArgs(trailingOnly=TRUE)

#source functions from Gautier file
source("./baypass_utils.R")

#load omega matrix
omega = as.matrix(read.table("prunedsnps.output_mat_omega.out")) # this is the mat_omega.out from prunned SNP run

#load beta parameters
pi.beta.coef=read.table("allsnps.output_summary_beta_params.out", h=T)$Mean  # this is the summary_beta_params.out from basic run

#create the POD - simulate
simu.bta <- simulate.baypass(omega.mat=omega, nsnp=5000,beta.pi=pi.beta.coef, suffix="simulates_pods")
