##############################
#  check the maf min and max #
##############################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/05_maf")
# load reference file with header in it. note here shared does not mean anything
pname = 'CHR20_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat <- read.delim(pname, header = TRUE, sep='\t')
max(dat$knownEM)
min(dat$knownEM)
hist(dat$knownEM)

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/05_maf")
# load reference file with header in it. note here shared does not mean anything
pname = 'HC_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat <- read.delim(pname, header = TRUE, sep='\t')
max(dat$knownEM)
min(dat$knownEM)
hist(dat$knownEM)

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/05_maf")
# load reference file with header in it. note here shared does not mean anything
pname = 'NB_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat <- read.delim(pname, header = TRUE, sep='\t')
max(dat$knownEM)
min(dat$knownEM)
hist(dat$knownEM)
