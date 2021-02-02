#!/bin/bash
start=`date +%s` 
## This script is used to deduplicate bam files and clipped overlapping read pairs for paired end data. It can process both paired end and single end data.
BAMLIST=/workdir/hz269/DelBay20/sample_lists/bam_list_1.txt # Path to a list of merged, deduplicated, and overlap clipped bam files. Full paths should be included.
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
PICARD=/programs/picard-tools-2.19.2/picard.jar # Path to picard tools
BAMUTIL=/programs/bamUtil/bam # Path to bamUtil
cd bam

## Loop over each sample
for SAMPLEBAM in `cat $BAMLIST`; do

	## Extract the file name prefix for this sample
	SAMPLESEQID=`echo $SAMPLEBAM | sed 's/_bt2_.*//' | sed -e 's#.*/bam/\(\)#\1#'`
	SAMPLEPREFIX=`echo ${SAMPLEBAM%.bam}`

	## Remove duplicates and print dupstat file
	# We used to be able to just specify picard.jar on the CBSU server, but now we need to specify the path and version
	java -Xmx60g -jar $PICARD MarkDuplicates I=$SAMPLEBAM O=$SAMPLEPREFIX'_dedup.bam' M=$SAMPLEPREFIX'_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

	## Clip overlapping paired end reads (only necessary for paired end data)
	$BAMUTIL clipOverlap --in $SAMPLEPREFIX'_dedup.bam' --out $SAMPLEPREFIX'_dedup_overlapclipped.bam' --stats

done

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
