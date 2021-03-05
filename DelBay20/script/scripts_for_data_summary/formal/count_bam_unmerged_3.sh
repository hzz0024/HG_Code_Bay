#!/bin/bash
start=`date +%s`
#!/bin/bash
BASEDIR=/workdir/hz269/DelBay_test/
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_1.txt
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt
REFNAME=CV30_masked
SAMTOOLS=/programs/samtools-1.11/bin/samtools
MINQ=20 # Minimum mapping quality filter used in the sorting step. When left undefined, the sorted bam file will not be counted (e.g. when reads are not filtered in the sorting step).

if [ -z "$MINQ" ]; then
	printf 'sample_id\tmapped_bases\n'
else
	printf 'sample_id\tmapped_bases\tqual_filtered_mapped_bases\n'
fi

for SAMPLEFILE in `cat $SAMPLELIST`; do

	# Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
	POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
	SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
	SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
	LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	SAMPLE_SEQ_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID

	## Extract data type from the sample table
	DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`

	## Count raw mapped bases
	RAWBAMFILE=$BASEDIR'bam/'$SAMPLE_SEQ_ID'_'$DATATYPE'_bt2_'$REFNAME'.bam'
	MAPPEDBASES=`$SAMTOOLS stats $RAWBAMFILE | grep ^SN | cut -f 2- | grep "^bases mapped (cigar)" | cut -f 2`

	if [ -z "$MINQ" ]; then
		printf "%s\t%s\n" $SAMPLE_SEQ_ID $MAPPEDBASES
	else
		## Count quality filtered mapped bases
		QUALFILTBAMFILE=$BASEDIR'bam/'$SAMPLE_SEQ_ID'_'$DATATYPE'_bt2_'$REFNAME'_minq'$MINQ'_sorted.bam'
		QUAFILTBASES=`$SAMTOOLS stats $QUALFILTBAMFILE | grep ^SN | cut -f 2- | grep "^bases mapped (cigar)" | cut -f 2`
		printf "%s\t%s\t%s\n" $SAMPLE_SEQ_ID $MAPPEDBASES $QUAFILTBASES
	fi
done

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
