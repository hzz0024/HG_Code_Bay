# angsd_pipeline for population genomics with low-coverage whole-genome sequencing

This pipeline was built by Claire Mérot with help and scripts from Eric Normandeau and many inputs from Anne-Laure Ferchaud and Amanda Xuereb.
It has been used for the analyses in the paper reference below, and it is flexible to be used for similar kind of analysis on different datasets.

Claire Mérot, Emma Berdan, Hugo Cayuela, Haig Djambazian, Anne-Laure Ferchaud, Martin Laporte, Eric Normandeau, Jiannis Ragoussis, Maren Wellenreuther, Louis Bernatchez, Locally adaptive inversions modulate genetic variation at different geographic scales in a seaweed fly, Molecular Biology and Evolution, 2021;, msab143
https://doi.org/10.1093/molbev/msab143



IMPORTANT: run all commands from the angsd_pipeline folder

# overview of the pipeline 
![overview_angsd_pipeline](https://github.com/clairemerot/angsd_pipeline/blob/master/Angsd_possibilities.jpg)

## 00_DEPENDANCIES
install angsd & associated programs
http://www.popgen.dk/angsd/index.php/ANGSD

add angsd to the path in .bashrc

add the misc folder (containing RealSFS, theta stat etc to the path in .bashrc

export PATH="/home/camer78/Softwares/angsd2/angsd:$PATH"
export PATH="/home/camer78/Softwares/angsd2/angsd/misc:$PATH"

install NGSAdmix (maybe in the misc folder, else export its path)
http://www.popgen.dk/software/index.php/NgsAdmix

install pcangsd (maybe in the misc folder) & check if you have python2
http://www.popgen.dk/software/index.php/PCAngsd
copy the path into 01_config.sh PCA_ANGSD_PATH=~/Softwares/pcangsd

for all script file, you may edit the header to put your email adress and adjust cpu/memory/time/allocation and slurm partition 

## 01_PREPARE_DATA

