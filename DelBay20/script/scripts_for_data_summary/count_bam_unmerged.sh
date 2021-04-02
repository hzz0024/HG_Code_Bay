#!/bin/bash
BASEDIR=/workdir/hz269/DelBay_test/
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_1.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table. An example of such a sample list is /workdir/cod/greenland-cod/sample_lists/sample_list_pe_1.tsv
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se. An example of such a sample table is: /workdir/cod/greenland-cod/sample_lists/sample_table.tsv
REFNAME=CV30_masked # Reference name to add to output files, e.g. gadMor2
SAMTOOLS=/programs/samtools-1.11/bin/samtools # Path to samtools
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