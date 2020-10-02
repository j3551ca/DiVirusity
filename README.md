# DiVirusity
Quantify viral genetic diversity within hosts from NGS data (paired-end raw reads) without needing to reconstruct full-length haplotypes.

3 diversity metrics are used (Shannon Diversity, Nucleotide (pi) Diversity, and Alternate Allele Frequency) to quanitfy the amount of genetic variation ocurring at each locus in the minority of a viral population infecting a host. If this is not your end goal, this script also contains code for 1) read trimming and reference-guided assembly (.fastq x 2 --> .bam) 2) generation of consensus sequences (.fasta) from alignment files (.bam) (i.e. a sequence of the most common allele at each position in an assembly) or 3) minority variant calling.  

This requires 2 other programs LoFreq (https://github.com/CSB5/lofreq) and SNPGenie (https://github.com/chasewnelson/SNPGenie).

Input: fasta reference genome (.fasta), fastq raw sequence reads (_R1.fastq, _R2.fastq), gtf file containing segments/ chromosomes matching reference (.gtf).

Useage: run in this folder (./DiVirusity/master) with raw sequencing data fastq files, reference genome, and gtf file as: ./diVirusity.sh reference gtf BQ MQ

Output: text files. 3 different diversity statistics for each position in the open reading frames of the virus. 

![https://github.com/j3551ca/DiVirusity/blob/main/divirusity.pdf]
