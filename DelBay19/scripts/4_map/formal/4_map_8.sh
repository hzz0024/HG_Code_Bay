start=`date +%s` 
BASEDIR=/workdir/hz269/DelBay20
BOWTIE=/programs/bowtie2-2.3.5.1-linux-x86_64/bowtie2
SAMTOOLS=/programs/samtools-1.11/bin/samtools
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_8.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
FASTQDIR=$BASEDIR/adapter_clipped/ # Path to the directory where fastq file are stored.
FASTQSUFFIX1=_adapter_clipped_f_paired.fastq.gz # Suffix to fastq files. Use forward reads with paired-end data.
FASTQSUFFIX2=_adapter_clipped_r_paired.fastq.gz # Suffix to fastq files. Use reverse reads with paired-end data.
MAPPINGPRESET=very-sensitive # The pre-set option to use for mapping in bowtie2 (very-sensitive for end-to-end (global) mapping [typically used when we have a full genome reference], very-sensitive-local for partial read mapping that allows soft-clipping [typically used when mapping genomic reads to a transcriptome]
REFERENCE=$BASEDIR/reference/CV30_masked.fasta # Path to reference fasta file and file name
REFNAME=CV30_masked # Reference name to add to output files, e.g. gadMor2

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do

	## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
	SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
	POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
	SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
	LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely

	## Extract data type from the sample table
	DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`

	## The input and output path and file prefix
	SAMPLETOMAP=$FASTQDIR$SAMPLE_UNIQ_ID
	SAMPLEBAM=$BASEDIR'/bam/'$SAMPLE_UNIQ_ID

	## Define platform unit (PU), which is the lane number
	PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`

	## Define reference base name
	REFBASENAME="${REFERENCE%.*}"

	## Map reads to the reference
	echo $SAMPLE_UNIQ_ID

	# Map the paired-end reads
	if [ $DATATYPE = pe ]; then
	# We ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (loose their mate so they become single-end), these can be mapped in a separate step and the resulting bam files merged with the paired-end mapped reads.
	$BOWTIE -q --phred33 --$MAPPINGPRESET -p 1 -I 0 -X 1500 --fr --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

	# Map the single-end reads
	elif [ $DATATYPE = se ]; then
	$BOWTIE -q --phred33 --$MAPPINGPRESET -p 1 --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -U $SAMPLETOMAP$FASTQSUFFIX1 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

	fi

	## Convert to bam file for storage (including all the mapped reads)
	$SAMTOOLS view -bS -F 4 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam' > $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam'
	rm -f $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

	## Filter the mapped reads (to onky retain reads with high mapping quality)
	# Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20)
	$SAMTOOLS view -h -q 20 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam' | $SAMTOOLS view -buS - | $SAMTOOLS sort -o $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'_minq20_sorted.bam'

done

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
