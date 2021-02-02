#!/bin/bash
start=`date +%s`
BAMLIST=/workdir/hz269/DelBay20/sample_lists/NY_dedup_overlapclipped.list # Path to a list of merged, deduplicated, and overlap clipped bam files. Full paths should be included.
BASEDIR=/workdir/hz269/DelBay20
REFERENCE=$BASEDIR/reference/CV30_masked.fasta # Path to reference fasta file and file name
SAMTOOLS=/programs/samtools-1.11/bin/samtools # Path to samtools
GATK=/programs/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar # Path to GATK

cd bam

## Loop over each sample
for SAMPLEBAM in `cat $BAMLIST`; do

if [ -e $SAMPLEBAM'.bai' ]; then
	echo "the file already exists"
else
	## Index bam files
	$SAMTOOLS index $SAMPLEBAM
fi
done

## Realign around in-dels
# This is done across all samples at once

## Use an older version of Java
export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH

## Create list of potential in-dels
if [ ! -f $BASEDIR'/bam/NY_for_indel_realigner.intervals' ]; then
		java -Xmx120g -jar $GATK \
			-T RealignerTargetCreator \
			-R $REFERENCE \
			-I $BAMLIST \
			-o $BASEDIR'/bam/NY_for_indel_realigner.intervals' \
			-drf BadMate
fi

## Run the indel realigner tool
java -Xmx120g -jar $GATK \
	-T IndelRealigner \
	-R $REFERENCE \
	-I $BAMLIST \
	-targetIntervals $BASEDIR'/bam/NY_for_indel_realigner.intervals' \
	--consensusDeterminationModel USE_READS  \
	--nWayOut _realigned.bam

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
