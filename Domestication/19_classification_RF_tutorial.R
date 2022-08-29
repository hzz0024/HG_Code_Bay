setwd("~/Dropbox/Mac/Documents/HG/Domestication/19_random_forest/doi_10.5061_dryad.k55hh8f__v1")
# Introductory tutorial for conducting Random Forest on data with a categorical response variable 
# (i.e. Random Forest using classification trees)
randomForest(x = X, y = as.factor(y), importance=TRUE ,proximity=TRUE, mtry=1000, ntree=1000, strata=as.factor(y), sampsize=32)
## NOTE: RUNNING THIS CODE IN ITS ENTIRETY TAKES APPROXIMATELY 4-8 HOURS ON TYPICAL LAPTOP AND DESKTOP COMPUTERS
## As a point of reference, a similar, empirical RF analysis for 400 individuals, 9000 loci, and a continuous trait 
## took approximately three weeks to run on a desktop computer with 24 GB of RAM and a processor speed of 3.47 GHz)

# The code relies on the package randomForest (Liaw and Wiener 2002)
install.packages("randomForest")
library(randomForest)

library(RcppCNPy)
X = npyLoad('/Users/ryan/PRF/X.npy')
y = npyLoad('/Users/ryan/PRF/y.npy')

# Import the example data set, which comprises 402 individuals genotyped at 1000 biallelic loci, where 0=homozygote 1, 1=heterozygote, 2=homozygote 2
# Each individual also has a binary phenotype - resistance to a disease - where 0=did not survive and 1=survived
# The objective is to identify loci associated with disease resistance 

class_data <- read.csv("data_classification_RF_tutorial.csv",row.names=1) #This is the simulated data from Table S8

# First, explore the overall distribution of the phenotype
hist(class_data$resistance)   

# Now examine distributions within the two populations
hist(class_data[which(class_data$Pop==1),2])
length(which(class_data$Pop==1 & class_data$resistance==0)) # 131 susceptible individuals in Pop 1
length(which(class_data$Pop==1 & class_data$resistance==1)) # 70 resistant individuals in Pop 1

hist(class_data[which(class_data$Pop==2),2])
length(which(class_data$Pop==2 & class_data$resistance==0)) # 68 susceptible individuals in Pop 2
length(which(class_data$Pop==2 & class_data$resistance==1)) # 133 resistant individuals in Pop 2

# Susceptibility/resistance to a disease appears to differ between the populations, and we should correct
# for this stratification before conducting RF to minimize the risk of false positive associations.
# Specifically, we will correct the genotypes using the approach of Zhao et al.(2012). 
# We will not correct the phenotype, as it is binary and we want to maintain its categorical distribution.

class_data_corrected <- class_data # create another data frame for the corrected genotypes
class_data_corrected[,3:1002] <- NA  # Keep columns 1-2 with population ID and phenotype, but then replace with NA's over which you can write the residuals.

# Now correct the genotypes using the regression/residual method. We're using a standard linear regression because Zhao et al. 2012 found that the correction procedure is robust to selection of the link function
for (i in 3:ncol(class_data)){
  LM_SNP_i <- lm(class_data[,i] ~ factor(class_data$Pop)) # apply linear model to all loci and the response
  class_data_corrected[,i] <- LM_SNP_i$residuals
  colnames(class_data_corrected)[i]<-colnames(class_data)[i] 
  if(i%%50==0) print(i)
}

# Verify that the residuals have been written to the data frame properly, using the last column as an example
class_data_corrected[,1002]-LM_SNP_i$residuals  #Should all be zero if correct

# Export a copy of the corrected data for future reference (This corresponds to the data in Table S9)
write.csv(class_data_corrected,file="data_classification_RF_tutorial_corrected.csv",row.names=FALSE)

# Before running Random Forest, let's also check for an imbalance in the response variable 
# because - as discussed in the manuscript - any imbalances can bias the results.
length(which(class_data$resistance==0)) # 199 susceptible individuals total
length(which(class_data$resistance==1)) # 203 resistant individuals total

# So the phenotypes are evenly distributed. However, if they were imbalanced, we would correct for the imbalance by balancing the representation of
# susceptible and resistant individuals within the training data. To do this, we would over-sample the underrepresented class
# and under-sample the overrepresented class. This would be performed by setting the sampsize parameter to 2/3
# of the class with the lower sample size and using those sample sizes in conjunction with the strata option 

# We will demonstrate how to do this even though the phenotypes are evenly distributed. 

# There are 199 susceptible individuals. 199*(2/3) = 133 individuals, so we will sample 133 susceptible
# individuals and 133 resistant individuals for the training data set, from which the trees are grown.
sample_size <- c(133,133)

# We will use these sample sizes in conjunction with the strata option (see below)

###########################################################################################################################################
###########################################################################################################################################

# Now run Random Forest analysis. Since this is a binary trait, we need to conduct a classification RF

# First, we need to optimize mtry by running different values of mtry at different values of ntree. 

# We will run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, and p, where p is the number of loci
# We will initially run each of these mtry values at ntree=100 to 1000 (by increments of 100). 
# We are looking for a plateau where the out-of-bag error rate (OOB-ER) stops decreasing with larger values of ntree
# Once we reach the plateau, we will choose the mtry value that minimizes the OOB-ER.

results_optimization <- matrix(data=NA , nrow = 0, ncol = 3)
for (i in seq(from = 100, to = 1000 , by = 100)){  # values of ntree
  print(i)
  for (j in c(32,64,100,200,333,1000)){    #values of mtry based on 1000 total loci
    rf_ij <- randomForest(x = class_data_corrected[,3:1002], y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, ntree=i, mtry=j, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$err.rate,1)[1]))
  }
}

# Clean up the file format
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","OOB_ER")

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 32],results_optimization$OOB_ER[results_optimization$mtry == 32], type="l", col="black", xlab="ntree",ylab="OOB-ER",ylim=c(0,1))
lines(results_optimization$ntree[results_optimization$mtry == 64],results_optimization$OOB_ER[results_optimization$mtry == 64], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 100],results_optimization$OOB_ER[results_optimization$mtry == 100], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 200],results_optimization$OOB_ER[results_optimization$mtry == 200], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 333],results_optimization$OOB_ER[results_optimization$mtry == 333], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 1000],results_optimization$OOB_ER[results_optimization$mtry == 1000], col="red")

# This plot shows that mtry=p is the best in terms of OOB-ER (although mtry=0.2p and p/3 are very similar), and that the OOB-ER has reached a plateau. 
# Therefore, we will use mtry=p for our Random Forest analyses. 
# Note that this plot differs from Figure 3 in the manuscript and thus demonstrates that optimal parameter values will vary based on each data set.

###########################################################################################################################################
###########################################################################################################################################
# Now begin the full Random Forest analyses

# Recall that even though we optimized mtry, we must now run a larger number of trees in order to achieve convergence of importance values between forests.
# As a starting point, we will grow 25,000 trees and increase if necessary. We do not need to worry about this increase in ntree affecting our mtry optimization,
# since the OOB-ER reached a plateau for a given mtry value after about 400 trees.

rf_all_1 = randomForest(x = class_data_corrected[,3:1002], y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=1000, ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_all_1,file="rf_all_1.Rdata")

rf_all_2 = randomForest(x = class_data_corrected[,3:1002], y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=1000, ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_all_2,file="rf_all_2.Rdata")

#Check correlation of locus importance values between forests 
importance_rf_all_1<-data.frame(importance(rf_all_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_all_1)<-c("importance")
importance_rf_all_2<-data.frame(importance(rf_all_2,type=1))
colnames(importance_rf_all_2)<-c("importance")

cor(importance_rf_all_1,importance_rf_all_2) # A correlation of 0.98 for locus importance values between forests is extremely good, so we'll use 25,000 trees for the remaining forests

rf_all_3 = randomForest(x = class_data_corrected[,3:1002], y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=1000, ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_all_3,file="rf_all_3.Rdata")
importance_rf_all_3<-data.frame(importance(rf_all_3,type=1))
colnames(importance_rf_all_3)<-c("importance")

############################################################################################################################################
############################################################################################################################################

# The predictive ability of classification trees is measured by the out-of-bag error rate.  An error rate is calculated for each tree within a forest. 
# We will use the error rate from the last tree in the forest, which takes all previous trees into account and thus represents the error rate after the model stabilizes/converges

rf_all_1_err.rate <- rf_all_1$err.rate[1000]
rf_all_2_err.rate <- rf_all_2$err.rate[1000]
rf_all_3_err.rate <- rf_all_3$err.rate[1000]

#Combine importance (mean decrease in accuracy) values of each locus across the three forests
importance_rf_all <-cbind(rownames(importance_rf_all_1),importance_rf_all_1,importance_rf_all_2, importance_rf_all_3)
colnames(importance_rf_all)<-c("Variable","Importance1","Importance2", "Importance3")

# Export importance values for future reference
write.csv(importance_rf_all,file="rf_importance_values_all_loci_classification_tutorial.csv",row.names=FALSE)
read.csv("rf_importance_values_all_loci_classification_tutorial.csv")
############################################################################################################################################
############################################################################################################################################

# Now conduct RF on subsets of the data to identify a group of loci that may be predictive of disease resistance. 
# For each subset, we will use mtry=p since that is the optimal setting that we previously found.

##### Best 2% 

names_best_2perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.98))]
names_best_2perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.98))]
names_best_2perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.98))]
names_best_2perc_unique<-unique(c(names_best_2perc_1,names_best_2perc_2,names_best_2perc_3))

# Extract genotypes 
genotypes_2perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_2perc_unique]

# Now conduct RF on this subset
rf_2perc_1 = randomForest(x = genotypes_2perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_2perc_1,file="rf_2perc_1.Rdata")

rf_2perc_2 = randomForest(x = genotypes_2perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_2perc_2,file="rf_2perc_2.Rdata")

rf_2perc_3 = randomForest(x = genotypes_2perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_2perc_3,file="rf_2perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_2perc_1_err.rate <- rf_2perc_1$err.rate[1000]
rf_2perc_2_err.rate <- rf_2perc_2$err.rate[1000]
rf_2perc_3_err.rate <- rf_2perc_3$err.rate[1000]

rm(rf_2perc_1,rf_2perc_2,rf_2perc_3) # remove the objects to save memory in R 


##### Best 3% 

names_best_3perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.97))]
names_best_3perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.97))]
names_best_3perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.97))]
names_best_3perc_unique<-unique(c(names_best_3perc_1,names_best_3perc_2,names_best_3perc_3))

# Extract genotypes
genotypes_3perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_3perc_unique]

# Now conduct RF on this subset
rf_3perc_1 = randomForest(x = genotypes_3perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_3perc_1,file="rf_3perc_1.Rdata")

rf_3perc_2 = randomForest(x = genotypes_3perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_3perc_2,file="rf_3perc_2.Rdata")

rf_3perc_3 = randomForest(x = genotypes_3perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_3perc_3,file="rf_3perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_3perc_1_err.rate <- rf_3perc_1$err.rate[1000]
rf_3perc_2_err.rate <- rf_3perc_2$err.rate[1000]
rf_3perc_3_err.rate <- rf_3perc_3$err.rate[1000]

rm(rf_3perc_1,rf_3perc_2,rf_3perc_3) # remove the objects to save memory in R 

##### Best 4% 

names_best_4perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.96))]
names_best_4perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.96))]
names_best_4perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.96))]
names_best_4perc_unique<-unique(c(names_best_4perc_1,names_best_4perc_2,names_best_4perc_3))

# Extract genotypes
genotypes_4perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_4perc_unique]

# Now conduct RF on this subset
rf_4perc_1 = randomForest(x = genotypes_4perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_4perc_1,file="rf_4perc_1.Rdata")

rf_4perc_2 = randomForest(x = genotypes_4perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_4perc_2,file="rf_4perc_2.Rdata")

rf_4perc_3 = randomForest(x = genotypes_4perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_4perc_3,file="rf_4perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_4perc_1_err.rate <- rf_4perc_1$err.rate[1000]
rf_4perc_2_err.rate <- rf_4perc_2$err.rate[1000]
rf_4perc_3_err.rate <- rf_4perc_3$err.rate[1000]

rm(rf_4perc_1,rf_4perc_2,rf_4perc_3) # remove the objects to save memory in R 

##### Best 5% 

names_best_5perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.95))]
names_best_5perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.95))]
names_best_5perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.95))]
names_best_5perc_unique<-unique(c(names_best_5perc_1,names_best_5perc_2,names_best_5perc_3))

# Extract genotypes
genotypes_5perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_5perc_unique]

# Now conduct RF on this subset
rf_5perc_1 = randomForest(x = genotypes_5perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_5perc_1,file="rf_5perc_1.Rdata")

rf_5perc_2 = randomForest(x = genotypes_5perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_5perc_2,file="rf_5perc_2.Rdata")

rf_5perc_3 = randomForest(x = genotypes_5perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_5perc_3,file="rf_5perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_5perc_1_err.rate <- rf_5perc_1$err.rate[1000]
rf_5perc_2_err.rate <- rf_5perc_2$err.rate[1000]
rf_5perc_3_err.rate <- rf_5perc_3$err.rate[1000]

rm(rf_5perc_1,rf_5perc_2,rf_5perc_3) # remove the objects to save memory in R 

##### Best 10% 

names_best_10perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.90))]
names_best_10perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.90))]
names_best_10perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.90))]
names_best_10perc_unique<-unique(c(names_best_10perc_1,names_best_10perc_2,names_best_10perc_3))

# Extract genotypes
genotypes_10perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_10perc_unique]

# Now conduct RF on this subset
rf_10perc_1 = randomForest(x = genotypes_10perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_10perc_1,file="rf_10perc_1.Rdata")

rf_10perc_2 = randomForest(x = genotypes_10perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_10perc_2,file="rf_10perc_2.Rdata")

rf_10perc_3 = randomForest(x = genotypes_10perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_10perc_3,file="rf_10perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_10perc_1_err.rate <- rf_10perc_1$err.rate[1000]
rf_10perc_2_err.rate <- rf_10perc_2$err.rate[1000]
rf_10perc_3_err.rate <- rf_10perc_3$err.rate[1000]

rm(rf_10perc_1,rf_10perc_2,rf_10perc_3) # remove the objects to save memory in R 

##### Best 20% 

names_best_20perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.80))]
names_best_20perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.80))]
names_best_20perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.80))]
names_best_20perc_unique<-unique(c(names_best_20perc_1,names_best_20perc_2,names_best_20perc_3))

# Extract genotypes
genotypes_20perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_20perc_unique]

# Now conduct RF on this subset
rf_20perc_1 = randomForest(x = genotypes_20perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_20perc_1,file="rf_20perc_1.Rdata")

rf_20perc_2 = randomForest(x = genotypes_20perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_20perc_2,file="rf_20perc_2.Rdata")

rf_20perc_3 = randomForest(x = genotypes_20perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_20perc_3,file="rf_20perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_20perc_1_err.rate <- rf_20perc_1$err.rate[1000]
rf_20perc_2_err.rate <- rf_20perc_2$err.rate[1000]
rf_20perc_3_err.rate <- rf_20perc_3$err.rate[1000]

rm(rf_20perc_1,rf_20perc_2,rf_20perc_3) # remove the objects to save memory in R 

##### Best 30% 

names_best_30perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.70))]
names_best_30perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.70))]
names_best_30perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.70))]
names_best_30perc_unique<-unique(c(names_best_30perc_1,names_best_30perc_2,names_best_30perc_3))

# Extract genotypes
genotypes_30perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_30perc_unique]

# Now conduct RF on this subset
rf_30perc_1 = randomForest(x = genotypes_30perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_30perc_1,file="rf_30perc_1.Rdata")

rf_30perc_2 = randomForest(x = genotypes_30perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_30perc_2,file="rf_30perc_2.Rdata")

rf_30perc_3 = randomForest(x = genotypes_30perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_30perc_3,file="rf_30perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_30perc_1_err.rate <- rf_30perc_1$err.rate[1000]
rf_30perc_2_err.rate <- rf_30perc_2$err.rate[1000]
rf_30perc_3_err.rate <- rf_30perc_3$err.rate[1000]

rm(rf_30perc_1,rf_30perc_2,rf_30perc_3) # remove the objects to save memory in R 


# Now combine all of the error rates from the subsets and for all loci to identify a group for the backward purging approach

All_initial_err.rate <- rbind(cbind(rf_all_1_err.rate,rf_all_2_err.rate,rf_all_3_err.rate),
                              cbind(rf_2perc_1_err.rate,rf_2perc_2_err.rate,rf_2perc_3_err.rate),
                              cbind(rf_3perc_1_err.rate,rf_3perc_2_err.rate,rf_3perc_3_err.rate),
                              cbind(rf_4perc_1_err.rate,rf_4perc_2_err.rate,rf_4perc_3_err.rate),
                              cbind(rf_5perc_1_err.rate,rf_5perc_2_err.rate,rf_5perc_3_err.rate),
                              cbind(rf_10perc_1_err.rate,rf_10perc_2_err.rate,rf_10perc_3_err.rate),
                              cbind(rf_20perc_1_err.rate,rf_20perc_2_err.rate,rf_20perc_3_err.rate),
                              cbind(rf_30perc_1_err.rate,rf_30perc_2_err.rate,rf_30perc_3_err.rate))

# Plot error rates for the various subsets
All_initial_err.rate<-data.frame(All_initial_err.rate)
All_initial_err.rate$Number_loci<-c(1000,length(names_best_2perc_unique),length(names_best_3perc_unique),length(names_best_4perc_unique),length(names_best_5perc_unique),length(names_best_10perc_unique),length(names_best_20perc_unique),length(names_best_30perc_unique))
rownames(All_initial_err.rate)<-c("All","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
All_initial_err.rate$Average<-apply(All_initial_err.rate[,1:3],1,mean)

# Write error rates to file for future reference
write.csv(All_initial_err.rate,file="All_initial_err_rate_classification_tutorial.csv")

# Plot error rates as well
par(mar=c(5,6,3,3))
plot(All_initial_err.rate$Number_loci,All_initial_err.rate$Average,log="x", pch=19,xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5)

# Based on this table and plot, the best 5% of loci have the lowest error rate
# As a conservative measure, I'll run backward purging RF with the best 10% loci

#################### Backward purging approach
names_purging <- names_best_10perc_unique

genotypes_purging<-class_data_corrected[,colnames(class_data_corrected) %in% names_purging]

rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_purging_1,file="rf_purging_1.Rdata")
rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_purging_2,file="rf_purging_2.Rdata")
rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_purging_3,file="rf_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[length(names_purging)]]<-names_purging
error_rate_best<-data.frame(V1=1:length(names_purging),V2=1:length(names_purging),V3=1:length(names_purging))
rownames(error_rate_best)<-1:length(names_purging)
error_rate_best[length(names_purging),] <- c(rf_purging_1$err.rate[1000],rf_purging_2$err.rate[1000],rf_purging_3$err.rate[1000])


for (i in 1:(length(names_purging)-2)){  # RF cannot be conducted with 1 locus, which is why the loop is from 1:length(names_purging)-2
  print(i)
  imp_purging_1<-data.frame(importance(rf_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging)
  colnames(imp_purging_1)<-"Mean_Decrease_Accuracy1"
  rownames(imp_purging_2)<-colnames(genotypes_purging)
  colnames(imp_purging_2)<-"Mean_Decrease_Accuracy2"
  rownames(imp_purging_3)<-colnames(genotypes_purging)
  colnames(imp_purging_3)<-"Mean_Decrease_Accuracy3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[length(names_purging)-i]]<-names_keep
  genotypes_purging<-class_data_corrected[,colnames(class_data_corrected) %in% names_keep]
  rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
  rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
  rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=1000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
  error_rate_best[length(names_purging)-i,] <- c(rf_purging_1$err.rate[1000],rf_purging_2$err.rate[1000],rf_purging_3$err.rate[1000])
}


error_rate_best$Average<-apply(error_rate_best,1,mean)
write.csv(error_rate_best, file="Backward_purging_OOB-ER_classification_tutorial.csv") # Save the error rates

# Now plot the backward purging results. Omit error rates from one locus since RF cannot be conducted with just one locus
plot(seq(2,nrow(error_rate_best),1),error_rate_best$Average[-c(1)],xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5,pch=16)

# Which group of loci yields the lowest error rate?
which(error_rate_best$Average==min(error_rate_best$Average[-c(1)])) #34 loci have the lowest OOB-ER

# Export the names of the predictor loci
write.csv(names_all_iterations[[34]],file="Predictor_loci_classification_tutorial.csv")
