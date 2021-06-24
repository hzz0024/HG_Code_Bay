sbatch --job-name=C_both --output=/local/workdir/hz269/DelBay_all_angsd/log/08_R_REF19_REF20_ARN_COH.log --nodes=1 --ntasks=10 --mem=40000 R_REF19_REF20_ARN_COH.sh
sbatch --job-name=C_batch --output=/local/workdir/hz269/DelBay_all_angsd/log/08_R_REF19_REF20_CHR19_CHR20.log --nodes=1 --ntasks=10 --mem=40000 R_REF19_REF20_CHR19_CHR20.sh
sbatch --job-name=C_19 --output=/local/workdir/hz269/DelBay_all_angsd/log/08_R_REF19_SR_ARN_COH.log --nodes=1 --ntasks=10 --mem=40000 R_REF19_SR_ARN_COH.sh
sbatch --job-name=C_20 --output=/local/workdir/hz269/DelBay_all_angsd/log/08_R_REF20_SR_ARN_COH.log --nodes=1 --ntasks=10 --mem=40000 R_REF20_SR_ARN_COH.sh
sbatch --job-name=F_1920 --output=/local/workdir/hz269/DelBay_all_angsd/log/08_R_REF19_CHR19_REF20_CHR20.log --nodes=1 --ntasks=10 --mem=40000 R_REF19_CHR19_REF20_CHR20.sh
sbatch --job-name=19_SR --output=/local/workdir/hz269/DelBay_all_angsd/log/08_R_REF19_CHR19_SR_HC.log --nodes=1 --ntasks=10 --mem=40000 R_REF19_CHR19_SR_HC.sh
sbatch --job-name=19_NB --output=/local/workdir/hz269/DelBay_all_angsd/log/08_R_REF19_CHR19_NB_HC.log --nodes=1 --ntasks=10 --mem=40000 R_REF19_CHR19_NB_HC.sh
sbatch --job-name=20_SR --output=/local/workdir/hz269/DelBay_all_angsd/log/08_R_REF20_CHR20_SR_HC.log --nodes=1 --ntasks=10 --mem=40000 R_REF20_CHR20_SR_HC.sh
sbatch --job-name=20_NB --output=/local/workdir/hz269/DelBay_all_angsd/log/08_R_REF20_CHR20_NB_HC.log --nodes=1 --ntasks=10 --mem=40000 R_REF20_CHR20_NB_HC.sh

