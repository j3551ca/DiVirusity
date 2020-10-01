#!/usr/bin/env bash

#this script determines *potential* coinfections by viruses included in bait design within samples. It works by mapping reads to baits/ genomes of viruses suspected to be present in samples and counting the number of occurences of each reference in the header of reads. As such, alignment parameters, sequencing type, and naming of baits are important. This does not guarantee a read originates from the virus it aligned to or that the number of reads is accurate (may be slightly off due to reference names existing in sam header if aligner used is different, for example), but rather provides a general, brief overview of putative co-infections which can then be more thoroughly analysed. 
#requires multifasta of viruses being searched as well as sample fastq files 
#dependencies: samtools, TrimmomaticPE, bwa


echo -e "\n	 -----------------------
	      CO-INFchecked
	 -----------------------\n"

#read -p "Enter directory containing raw read data (.fastq) & multifasta of putative co-infecting viruses (.fasta):" root_dir
read -p "Enter name of multifasta containing sequences of putative co-infecting viruses:" multi_fasta

#cd ${work_dir}

fastq=(${root_dir}/*.fastq)
if [[ ${#fastq[@]} -lt 2 ]]; 
then 
	echo -e "ERROR: insufficient fastq files. Exiting...\a"; 
	exit; 
fi

fasta=(${root_dir}/*.fasta)
# if [[ ${#fasta[@]} -gt 1 ]] ;
	
if [[ ! -f "${multi_fasta}" ]];
then 
	echo -e "ERROR: no reference genome found. Ensure correct spelling and location. Exiting...\a"; 
	exit; 
elif [[ ${#fasta[@]} -lt 1 ]]; 
then 
	echo -e "ERROR: no reference genome found. Ensure .fasta extension. Exiting...\a"; 
	exit; 
else
	echo -e "Multi-fasta found. Commencing reference-guided alignment..."
fi


#retrieve names of viruses in fasta baits file
grep '>' ${multi_fasta} | cut -d '_' -f 1 | cut -c 2- | cut -d '-' -f 1 | sort -u >>virus_id.txt

#aligning to baits of multiple viruses
bwa index ${multi_fasta}

#this loop may already be done by QuasDivGui.sh - unnecessary inclusion
# for r1_file in *_R1_*.fastq;
# do
# 	r2_file=${r1_file/_R1_/_R2_};
# 	TrimmomaticPE -threads 2 ${root_dir}/$r1_file ${root_dir}/$r2_file ${root_dir}/${r1_file%.*}_paired.fq \
# 	${root_dir}/${r1_file%.*}_unpaired.fq ${root_dir}/${r2_file%.*}_paired.fq ${root_dir}/${r2_file%.*}_unpaired.fq \
# 	ILLUMINACLIP:/home/jessica/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;
# done

#can start here... 
# for r1_file in *_paired.fq;
# do
# 	/home/gideon/FastQC/fastqc ${root_dir}/${r1_file};
# done

#...or here. Take trimmed reads and mass-align them to the multifasta to see where they map to..coinfection starting point takes only
#uniquely aligned reads..though in reality, a read can be from another co-infecting virus even if it maps to multiple different genomes #homology #recombination 
#but it's a start.. 
for refine_reads_for in *_R1_*_paired.fq; 
do 
	refine_reads_rev=${refine_reads_for/_R1_/_R2_};
	prefix=${refine_reads_for/_L001_R1_001_paired.fq/_}; 
	bwa mem -t 6 ${root_dir}/${multi_fasta} ${root_dir}/${refine_reads_for} ${root_dir}/${refine_reads_rev} | \
		samtools view -bS -F 0x04 - > ${root_dir}/${prefix}coinf.bam; 
done

#try samtools view -bS -F 2308 - > coinf.bam (do not include unmapped reads, reads that are secondary alignment, or reads that are part of supplementary alignment)
#from https://broadinstitute.github.io/picard/explain-flags.html 


for coinf_aln in *_coinf.bam; 
do 
	samtools sort -o ${coinf_aln%.*}.sort.bam ${coinf_aln}; 
	samtools index -b ${coinf_aln%.*}.sort.bam; 
	#samtools view -h ${coinf_aln%.*}.sort.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' > ${coinf_aln%.*}_unique.sort.sam
	samtools view ${coinf_aln%.*}.sort.bam | grep 'XT:A:U' | samtools view -S -T ${multi_fasta} - ${coinf_aln%.*}uniqueMap.sam
done 

rm *_coinf.bam *_coinf.sort.bam

#samtools view -O SAM -o ${r3_file%.*}.sort.sam ${r3_file%.*}.sort.bam; <<<---original
#samtools view ${r3_file%.*}.sort.bam | grep 'XT:A:U' | samtools view -bS -T referenceSequence.fa - > reads.uniqueMap.bam
#samtools view -h ${r3_file%.*}.sort.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > ${r3_file%.*}_unique_map.sort.bam

#distinguish virus name input from output
cp ./virus_id.txt ./coinfection_count.txt

#output and print # reads matching specified virus:
for sam in *.sam; 
do 
	echo -e "\n\r${sam}" >> coinfection_count.txt; 
	while read -r line; 
	do 
		virus=$(grep -ic "${line}" ${sam}); 
		echo -e "\t${virus}" >> coinfection_count.txt; 
	done < virus_id.txt; 
done 

#transpose
awk '{ for (i=1;i<=NF;i++ ) printf $i "\t" }' coinfection_count.txt > coinfection_count_final.txt
sed -i "1s/^/\t/" coinfection_count_final.txt

echo -e "See coinfection_count_final.txt for results."



#declare -a virus
#vir_len=`wc -l virus_id.txt`
#for sam in *.sam; do for i in vir_len; do virus[i]=$(grep -ic "${line}" ${sam}); 
