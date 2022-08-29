start=`date +%s`  ## date at start
BASEDIR=/workdir/hz269/DelBay19
TRIMMOMATIC=/programs/trimmomatic/trimmomatic-0.39.jar
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_5.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
RAWFASTQDIR=$BASEDIR/raw_fastq/ # Path to raw fastq files.
RAWFASTQSUFFIX1=_R1.fastq.gz # Suffix to raw fastq files. Use forward reads with paired-end data.
RAWFASTQSUFFIX2=_R2.fastq.gz # Suffix to raw fastq files. Use reverse reads with paired-end data.
ADAPTERS=$BASEDIR/reference/NexteraPE-PE.fa # Path to a list of adapter/index sequences, copied from /programs/bbmap-38.86/resources/adapters.fa or /programs/trimmomatic/adapters/TruSeq3-PE-2.fa

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do
	## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
	SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
	POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
	SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
	LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely
	echo "Sample: $SAMPLE_UNIQ_ID"
	## Extract data type from the sample table
	DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`

	## The input and output path and file prefix
	RAWFASTQ_ID=$RAWFASTQDIR$SAMPLEFILE
	SAMPLEADAPT=$BASEDIR'/adapter_clipped/'$SAMPLE_UNIQ_ID

	## Adapter clip the reads with Trimmomatic
	# The options for ILLUMINACLIP are: ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
	# The MINLENGTH drops the read if it is below the specified length in bp
	# For definitions of these options, see http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

	if [ $DATATYPE = pe ]; then
		java -jar $TRIMMOMATIC PE -threads 1 -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $RAWFASTQ_ID$RAWFASTQSUFFIX2 $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_f_unpaired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_unpaired.fastq.gz' 'ILLUMINACLIP:'$ADAPTERS':2:30:10:1:true MINLENGTH:80'

	elif [ $DATATYPE = se ]; then
		java -jar $TRIMMOMATIC SE -threads 1 -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $SAMPLEADAPT'_adapter_clipped_se.fastq.gz' 'ILLUMINACLIP:'$ADAPTERS':2:30:10 MINLENGTH:40'
		fi

done

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
