These are some scripts that we use in the lab for Hi-C alignment. 

There are three scripts for alignment. One is a shell scripts that just runs the other two in a single pipeline (hic_align_pipeline.sh)

One of the two perl scripts will align a single fastq file to the reference genome using BWA mem and perform some simple post-processing of the alignment (bwa_mem_hic_aligner.pl)
The other will take two bam files (one for read1, one for read 2), and manually merge them into a single bam file with paired information retained for each read (two_read_bam_combiner.pl).

For everything to work, you need bwa and samtools installed and in your path.

You can run the whole pipeline as follows:
./hic_align_pipeline.sh <read1> <read2> <reference genome> <number of threads> <output file name - will be *.bam>

Here you provide two fastq files for read1 and read2 (can also be gzipped), a reference genome that is indexed by bwa (i.e.: bwa index -a bwtsw my_reference_genome.fa), the number of threads to use for alignment, and the prefix for the output file name
These will all be passed in the shell script to the bwa_mem_hic_aligner.pl and two_read_bam_combiner.pl perl scripts.

You can also run them separately or as part of your own shell script.
To run the bwa_mem_hic_aligner.pl script, you run the command as follows:
./bwa_mem_hic_aligner.pl <fastq file> <reference genome> <threads>
In this case, you are giving the script one fastq file (either read1 or read2, can be gzipped), a refernece genome indexed by bwa, and the number of threads to use for alignment. The output will go to standard output in sam format. In reality we never actually save it in sam format, and instead immediately pipe the results to samtools to conver to bam, such as:
./bwa_mem_hic_aligner.pl <fastq file> <reference genome> <threads> | samtools view -bS -o output_file_name_here.bam -

To run the two_read_bam_combiner.pl perl script, run the following command:
./two_read_bam_combiner.pl --file1 read1_bam_file --file2 read2_bam_file
This script has a few other parameters you can invoke depending on what you want your output file to contain. They are as follows:

--qual N    //This will only keep reads that align with a mapping quality greater than or equal to N

--single    //This will retain single reads as well. In other words, in the event that one read aligns with qual>=N and its pair does not, this will keep the read as a single stranded read in the bam file. The frequency of such reads can increase near unmappable regions such as the centromere or telomere, so it may be a good idea to keep them for bias/normalization purposes. Also, if using Hi-C for phasing these can occaisonally still provide useful information if the read contains two heterozygous SNPs.

--bed_file <bed file> //If you provide a bed file, the script will only keep reads that align within this bed file

--paired_bed <6 column bed-like file> //If provided, will filter any reads that align between the two regions. One of the only reasons I have ever used this is if there is concern about interactions between highly similar regions like segmental duplications

--frag_file <bed file of restriction enzyme fragments> //If provided, will add a tag in the bed file containing the restriction enzyme fragment. Useful in the event you want to know if analyzing fragment level data

--chr_file <list of chromosomes> //If provided, will only consider chromosomes in this file as mapped. Perhaps useful if you want to filter out reads aligning to unplaced contigs. 


Notes on why we use this scripts (as opposed to other bwa based alignment methods):
The short answer is that it is historical, we have been using these scripts or versions of them for a while.
One particular reason to use this for the hic_breakfinder tool is filtering out of poorly mapped reads. This can make a difference in terms of Hi-C based SV calls. Poorly mapped reads are not inherently filtered out by using "bwa mem -SP5M", but can be with other post-processing scripts.

There are two additional scripts for generating a contacts file from the bam files, and then for also merging contact files (i.e. from replicate experiments).

