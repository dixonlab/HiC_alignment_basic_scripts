#!/bin/sh

### $1 is read1
### $2 is read2
### $3 is reference genome
### $4 is number of threads
### $5 is the output name

if [ $# != 5 ]; then

echo "Usage: ./hic_align_pipeline.sh <read1> <read2> <reference genome> <number of threads> <output file name - will be *.bam>"

else 

./bwa_mem_hic_aligner.pl $1 $3 $4 |\
samtools view -bS -o $5\_read1.bam -

./bwa_mem_hic_aligner.pl $2 $3 $4 |\
samtools view -bS -o $5\_read2.bam -

./two_read_bam_combiner.pl --qual 30 --single --file1 $5\_read1.bam --file2 $5\_read2.bam |\
samtools view -u -o - - |\
samtools sort -T $5 -o $5.bam -

fi
