cat haplotig_masked_genome.fasta |\
perl -pe '/^>/ ? s/.*// : s/[\n\r]//g' |\
tail -n +2 |\
perl -pe 's/.*/length($&)/e'
