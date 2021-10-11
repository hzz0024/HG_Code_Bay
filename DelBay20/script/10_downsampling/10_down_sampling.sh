#!/bin/bash

BASEDIR=/workdir/hz269/DelBay_all_angsd/10_downsamping/ # Path to base directory.
TARGETDIR=/workdir/hz269/DelBay_all_angsd/10_downsamping/down_sampling_1x_final/ # Path to target directory for data store.
BAMLIST=/workdir/hz269/DelBay_all_angsd/10_downsamping/sample_lists/bam_list.txt # Path to a list of bam files for downsampling (after overlap clipping and realignment around indels).
SAMTOOLS=/programs/samtools-1.11/bin/samtools # Path to Samtools
PICARD=/programs/picard-tools-2.19.2/picard.jar # Path to GATK
SUM=DelBay20_dedup_overlapclipped_realign_depth_per_position_per_sample_summary.tsv # Path to the coverage summary of bam files (before downsampling)
Del19_cvg=0.777531174804435 # This is the average coverage across 2019 samples. I used the genome-wide coverage as it is more comparable than the realized coverage.

for SAMPLEBAM in `cat $BAMLIST`; do

    SAMPLEPREFIX=`echo ${SAMPLEBAM%.bam}`
    echo "Sample: $SAMPLEPREFIX"
    SAMPLE=`echo $SAMPLEBAM | cut -d$'_' -f 1-2`
    DEPTH=`grep $SAMPLE $SUM | cut -f 2`
    pct=`awk "BEGIN {printf \"%.2f\n\", $Del19_cvg/$DEPTH}"`
    
    echo "start downsampling"
    if (( $(echo "$pct >= 1" |bc -l) )); 
    then
        echo "Sample $SAMPLE has smaller coverage than the Del19 target coverage of $Del19_cvg"
        cp $SAMPLEBAM $TARGETDIR$SAMPLEPREFIX'.1.bam'
        echo "--------------------------------------"
    else
        java -jar $PICARD DownsampleSam \
            I=$BASEDIR$SAMPLEBAM \
            O=$TARGETDIR$SAMPLEPREFIX'.'$pct'.bam' \
            STRATEGY=Chained \
            RANDOM_SEED=1 \
            P=$pct \
            ACCURACY=0.0001
        echo "--------------------------------------"
    fi
done
