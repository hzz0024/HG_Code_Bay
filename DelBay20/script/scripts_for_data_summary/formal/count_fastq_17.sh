#!/bin/bash
start=`date +%s`
#!/bin/bash

## This script is used to count number of bases in raw, adapter clipped, and quality filtered fastq files. The result of this script will be stored in a nohup file.
BASEDIR=/workdir/hz269/DelBay_test/
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_17.txt
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt
RAWFASTQDIR=/workdir/hz269/DelBay_test/raw_fastq/
QUALFILTERED=true # Whether the sample has gone through quality filtering. true or false

# Create headers for the output
if $QUALFILTERED; then
	printf 'sample_seq_id\traw_reads\traw_bases\tadapter_clipped_bases\tqual_filtered_bases\n'
else
	printf 'sample_seq_id\traw_reads\traw_bases\tadapter_clipped_bases\n'
fi

# Loop over each sample in the sample table
for SAMPLEFILE in `cat $SAMPLELIST`; do
	RAWFASTQFILES=$RAWFASTQDIR'/'$SAMPLEFILE'*.gz'  # The input path and file prefix

	# Count the number of reads in raw fastq files. We only need to count the forward reads, since the reverse will contain exactly the same number of reads. fastq files contain 4 lines per read, so the number of total reads will be half of this line number.
	RAWREADS=`zcat $RAWFASTQFILES | wc -l`

	# Count the number of bases in raw fastq files. We only need to count the forward reads, since the reverse will contain exactly the same number of bases. The total number of reads will be twice this count.
	RAWBASES=`zcat $RAWFASTQFILES | awk 'NR%4==2' | tr -d "\n" | wc -m`

	# Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
	POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
	SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
	SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
	LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	SAMPLE_SEQ_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID

	# Find all adapter clipped fastq files corresponding to this sample and store them in the object ADAPTERFILES.
	ADAPTERFILES=$BASEDIR'adapter_clipped/'$SAMPLE_SEQ_ID'*.gz'

	# Count all bases in adapter clipped files.
	ADPTERCLIPBASES=`zcat $ADAPTERFILES | awk 'NR%4==2' | tr -d "\n" | wc -m`

	# If reads are quality filtered, count quality filtered files.
	if $QUALFILTERED; then

		# Find all quality trimmed fastq files corresponding to this sample and store them in the object QUALFILES.
		QUALFILES=$BASEDIR'qual_filtered/'$SAMPLE_SEQ_ID'*.gz'

		# Count bases in quality trimmed files.
		QUALFILTPBASES=`zcat $QUALFILES | awk 'NR%4==2' | tr -d "\n" | wc -m`

		# Write the counts in appropriate order.
		printf "%s\t%s\t%s\t%s\t%s\n" $SAMPLE_SEQ_ID $((RAWREADS/4)) $RAWBASES $ADPTERCLIPBASES $QUALFILTPBASES

		# When reads are not quality filtered, directly write the output
	else

		# Write the counts in appropriate order.
		printf "%s\t%s\t%s\t%s\s\n" $SAMPLE_SEQ_ID $((RAWREADS/4)) $RAWBASES $ADPTERCLIPBASES

	fi

done

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
