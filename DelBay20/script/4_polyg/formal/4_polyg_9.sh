#!/bin/bash
start=`date +%s`  ## date at start
## This script is used to quality filter and trim poly g tails. It can process both paired end and single end data.
BASEDIR=/workdir/hz269/DelBay20
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_9.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table. 
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
FILTER=polyg # Type of filtering. Values can be: polyg (forced PolyG trimming only), quality (quality trimming, PolyG will be trimmed as well if processing NextSeq/NovaSeq data), or length (trim all reads to a maximum length)
THREAD=1 # Number of thread to use. Default is 10
FASTP=/programs/fastp-0.20.0/bin/fastp ## Path to the fastp program. The default path is /workdir/programs/fastp_0.19.7/fastp
MAXLENGTH=100 # Maximum length. This input is only relevant when FILTER=length, and its default value is 100.

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do

	## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
	SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
	POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
	SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
	LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID
	echo "Sample: $SAMPLE_UNIQ_ID"
	## Extract data type from the sample table
	DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`
	## The input and output path and file prefix
	SAMPLEADAPT=$BASEDIR'/adapter_clipped/'$SAMPLE_UNIQ_ID
	SAMPLEQUAL=$BASEDIR'/qual_filtered/'$SAMPLE_UNIQ_ID

	## Trim polyg tail or low quality tail with fastp.
	# --trim_poly_g forces polyg trimming, --cut_right enables cut_right quality trimming
	# -Q disables quality filter, -L disables length filter, -A disables adapter trimming
	# Go to https://github.com/OpenGene/fastp for more information
	if [ $DATATYPE = pe ]; then
		if [ $FILTER = polyg ]; then
			$FASTP --trim_poly_g --cut_right -L -A --thread $THREAD -i $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' -I $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_f_paired.fastq.gz' -O $SAMPLEQUAL'_adapter_clipped_qual_filtered_r_paired.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
		elif [ $FILTER = quality ]; then
			$FASTP -L -A --thread $THREAD -i $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' -I $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_f_paired.fastq.gz' -O $SAMPLEQUAL'_adapter_clipped_qual_filtered_r_paired.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
		elif [ $FILTER = length ]; then
			$FASTP --max_len1 $MAXLENGTH -Q -L -A --thread $THREAD -i $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' -I $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_f_paired.fastq.gz' -O $SAMPLEQUAL'_adapter_clipped_qual_filtered_r_paired.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
		fi
	fi
done

end=`date +%s` ## date at end
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
