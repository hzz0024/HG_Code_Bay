#path to bam list
CH="/workdir/hz269/02_info/ch_50.list"
REF="/workdir/hz269/02_info/ref_48.list"
CHR="/workdir/hz269/02_info/challenge_98.list"
WILD="/workdir/hz269/02_info/wild_235.list"
HC="/workdir/hz269/02_info/HC_48.list"
ARN="/workdir/hz269/02_info/ARN_47.list"
COH="/workdir/hz269/02_info/COH_44.list"
SR="/workdir/hz269/02_info/SR_48.list"
NB="/workdir/hz269/02_info/NB_48.list"
HC_ARN="/workdir/hz269/02_info/HC_ARN_95.list"
ARN_COH="/workdir/hz269/02_info/ARN_COH_91.list"
COH_SR="/workdir/hz269/02_info/COH_SR_92.list"
SR_NB="/workdir/hz269/02_info/SR_NB_96.list"
test="/workdir/hz269/02_info/test.list"

#path to the anc genome
ANC="/workdir/hz269/genome/cv30.fa"
ANC_MASKED="/workdir/hz269/genome/cv30_masked.fasta"

#path to bam folder
#BAM_PATH=../02_info

#path to pcaangsd
#PCA_ANGSD_PATH="/project/lbernatchez/programs/pcangsd "

#filter : will keep SNP above this allele frequency (over all individuals)
MIN_MAF=0.05 

#filter : will keep SNP with at least one read for this percentage of individuals (over all individuals in step 03, and within each pop at step 07)
PERCENT_IND=0.7

#filter: will keep SNP with at least a coverage of this factor multiplied by the number of ind - across all ind. usually set 2-4
#times the coverage to remove repeated regions
#MAX_DEPTH_FACTOR=3

#window size for sliding window FST & Thetas
WINDOW=200 

#window step
WINDOW_STEP=200 

#min nb of pop to consider for NGS admix
K_MIN=2 

#maximum nb of pop to consider for NGS admix
K_MAX=5 
