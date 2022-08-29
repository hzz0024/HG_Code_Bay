# Data processing for population genomics with low-coverage whole-genome sequencing

This pipeline was built by Honggang Zhao with help and scripts from Nina's lab and many inputs from Claire Mérot github site.

# overview of the pipeline 
![overview_angsd_pipeline](https://github.com/clairemerot/angsd_pipeline/blob/master/Angsd_possibilities.jpg)

## 00_DEPENDANCIES

A folder named DBW_NEH_crossing_data_processing, and run all commands from the DBW_NEH_crossing_data_processing folder

add angsd to the path in .bashrc

add the misc folder (containing RealSFS, theta stat etc to the path in .bashrc

export PATH="/home/camer78/Softwares/angsd2/angsd:$PATH"
export PATH="/home/camer78/Softwares/angsd2/angsd/misc:$PATH"

## 01_PREPARE_DATA

