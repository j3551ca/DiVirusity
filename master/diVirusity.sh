#!/usr/bin/env bash

#diVirusity#
#Description: script for analysis of Illumina short read viral sequences corrected for sequencing error to produce normalized shannon diversity, nucleotide frequency, and % non-ref base mutations at each position in viral genome. Input files are fastq and fasta. Output files are txt, csv, and python fig. 
#Dependencies: bcftools >v1.9, bedtools, bwa, LoFreq, samtools >v1.9, SNPGenie, TrimmomaticPE, python >= v2.7, **Vince Montoya's python scripts
#Input: fasta reference genome, fastq raw sequence reads, gtf file containing segments/ chromosomes matching reference and corresponding lengths (see Test_data file),
#Useage: run in this folder with raw seq data fastq files, reference genome, and gtf file as: ./diVirusity.sh reference gtf BQ MQ

#bash continually check exit status of commands and exit program if failed. 
#set -e

echo -e "\n	 -----------------------
	    Welcome to DiVirusity
	 -----------------------\n"

#Instead of specifying working directory..mv files into same directory as script..then can change to shannon div sub-folder and ntdiv-folder more easily (less chance for error)


#renaming positional parameters
ref_fasta=${1}
gtf=${2}
BQ=${3:-30}
MQ=${4:-20}

if [[ $# -lt 4 ]]; 
then
	echo -e "ERROR: missing at least one argument. Ensure input is: script reference gtf BQ MQ. \nTo use default BQ or MQ enter << - >> in place of integer. Exiting...";
	exit;
fi
	
#absolute working directory
#curr_loc=`realpath -e "$0"`
#root_dir=`dirname "${curr_loc}"`
root_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )
echo -e "Commencing analysis in ${root_dir}. Ensure raw read data (.fastq) & reference genome (.fasta/.fa) are located here.\n"

#SNPGenie directory absolute path*******************************how? find / -type d -name SNPGenie-master 2> /dev/null
read -p "Enter absolute path to SNPGenie-master:" snp_gen_dir

if [[ -z "$snp_gen_dir" ]];
then 
	echo -e "\nAbsolute path to SNPGenie-master not entered. Searching directories.\n 
This may take some time...\n" 
	snp_gen_dir=`find / -type d -name SNPGenie-master 2>/dev/null`
elif [[ -z "$snp_gen_dir" ]];
then
	echo -e "\nCould not locate SNPGenie-master directory. Ensure correct installation and spelling of program directory.\n
	Exiting..."
	exit;
else
	echo -e "Absolute path to SNPGenie-master set."
fi


#cd ${root_dir}

shopt -s nullglob dotglob 

#expand zipped files if they exist; dumping error if they don't. 
gunzip *.fastq.gz 2>/dev/null

#testing for existance of necessary files
fastq=(${root_dir}/*.fastq)


if [[ ${#fastq[@]} -lt 2 ]];
then 
	echo -e "ERROR: insufficient fastq files. Exiting...\a"; 
	exit;
elif [[ ! -f "${ref_fasta}" ]];
then
	echo -e "ERROR: no reference genome found. Ensure correct spelling and location. Exiting...\a";
	exit;
else
	echo -e "Commencing intra-host genetic diversity analysis.\n"
fi

#make sra fastq files compatible with rest of downstream analysis. Above line searches for given patterns. If non-sra-derived reads are named similarly, they 
#will simply be renamed. Redundant, but not fatal to program. 
# # sra_fastq=(${root_dir}/*_R1.fastq)
#
# if [[ ${#sra_fastq[@]} -ge 1 ]];
# then
# for file_for in ${sra_fastq[@]};
# 	do
# 		file_rev=${file_for/_R1/_R2};
# 		mv ${file_for} ${file_for/_R1/_L001_R1_001};
# 		mv ${file_rev} ${file_rev/_R2/_L001_R2_001};
# 	done
# fi

#alignment with new samtools (v1.9)
bwa index ${root_dir}/${ref_fasta}

echo "-----------Processing reads...-----------"


for raw_read_for in *_R1_*.fastq; 
do 
	raw_read_rev=${raw_read_for/_R1_/_R2_}; 
	TrimmomaticPE -threads 2 ${root_dir}/${raw_read_for} ${root_dir}/${raw_read_rev} \
		${root_dir}/${raw_read_for%.*}_paired.fq ${root_dir}/${raw_read_for%.*}_unpaired.fq \
			${root_dir}/${raw_read_rev%.*}_paired.fq ${root_dir}/${raw_read_rev%.*}_unpaired.fq \
				ILLUMINACLIP:./TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36; 
done

echo "-----------Performing reference-guided alignment...-----------"

for filt_reads_for in *_R1_*_paired.fq; 
do 
	filt_reads_rev=${filt_reads_for/_R1_/_R2_}; 
	prefix=${filt_reads_for/_L001_R1_001_paired.fq/_}; 
	bwa mem -t 6 ${root_dir}/${ref_fasta} ${root_dir}/${filt_reads_for} ${root_dir}/${filt_reads_rev} | \
		samtools view -bS -F 0x04 - > ${root_dir}/${prefix}.bam; 
done

for aln_bam in *_.bam; 
do 
	samtools sort -o ${aln_bam%.*}.sort.bam ${aln_bam}; 
	samtools index -b ${aln_bam%.*}.sort.bam; 
done 

#cp *_.bam ${ref_fasta} ./shannon_div/

#create bed file from gtf:
cut -f1,4,5 ${gtf} > ${gtf%.*}.bed

#generate coverage csv files for read depth coverage plots and read_info.sh. Command includes only higher quality bases from ORFs. 
for sort_aln_bam in *.sort.bam; 
do 
	samtools depth -a -b ${gtf%.*}.bed -d 0 -q ${BQ} -Q ${MQ} ${sort_aln_bam} > ${sort_aln_bam%.*}.cov; 
done

echo -e "-----------Obtaining Mode Nucleotide Distributions of Alignment Files-----------"
for aln in *.sort.bam; 
do 
	samtools stats ${aln} | grep ^GCC | cut -f 2- >> ACGT_distr.txt ; 
done


#remove - in filename and replace with _
#for file in *.cov; do mv "${file}" "${file/-/_}"; done


export root_dir
./read_info_formatted.sh 
#sleep 30
./co_infchecked_formatted.sh

wait

#make file containing segments and respective lengths
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' ${ref_fasta} | \
	awk '{print $1}' | awk 'sub(/^>/, "") 1' | awk '{printf "%s%s",$0,NR%2?"\t":RS}' > seglen.txt

#option to comment this portion out and import Geneious-generated cns seqs
#consensus sequence generation with new v1.9 samtools/bcftools commands 

echo -e "-----------Extracting majority consensus sequences-----------"
for sort_aln_bam in *.sort.bam; 
do 
	prefix=${sort_aln_bam%.*}; 
	bcftools mpileup -d1000000 -q20 -Q30 -f ${root_dir}/${ref_fasta} ${sort_aln_bam} | \
		bcftools call -c -O z > ${prefix}_call.vcf.gz; 
	#lofreq call -f ${root_dir}/${prefix}.sort_final_cns.fa -o ${sort_quasi_bam%.*}_rawlofreq.vcf -q 30 -Q 30 ${sort_quasi_bam} --verbose
done

#filter sites to include those where alt allele is majority and obtain sites with 0 coverage
for raw_vcf in *_call.vcf.gz; 
do 
	prefix=${raw_vcf/_call.vcf.gz/_}; 
	bcftools filter -O z -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3])' ${raw_vcf} >${prefix}filt.vcf.gz; 
	bedtools genomecov -bga -i ${raw_vcf} -g ${root_dir}/seglen.txt | \
		grep -w 0$  > ${prefix}nocovreg.bed; 
done

#make raw/ template cns seqs (with bcftools algorithm that normally subs sites of 0 coverage with reference base at that position)
for filt_vcf in *_filt.vcf.gz; 
do 
	prefix=${filt_vcf/.vcf.gz/_}; bcftools index ${filt_vcf}; 
	bcftools consensus -f ${root_dir}/${ref_fasta} ${filt_vcf}  > ${prefix}cns.fa; 
done 

#apply sites of 0 coverage to template cns seqs. This shows N at positions of 0 cov, no IUPAC bases. Caveats = if 2 reads at one site, will
#choose a base (50:50), if there is low coverage at a site (ex. 1 read) may count it as N. Important to know coverage of samples beforehad to
#account for this before this cns-generation stage.  
for raw_cns in *_cns.fa; 
do 
	prefix=${raw_cns/_filt_cns.fa/_}; 
	bedtools maskfasta -fi ${raw_cns} -fo ${prefix}final_cns.fa -bed ${prefix}nocovreg.bed; 
done

#concatente orfs to be used in phylogenetic tree analysis NOTE: bed file ORF starts need to be -1 since bedtools getfasta cuts off start (inclusive)
#for cns_fa in *.sort_final_cns.fa; do bedtools getfasta -fi ${cns_fa} -bed PRV_CDS.bed -split -fo ${cns_fa/.sort_final_cns.fa/orf_concat.fa}; done
#for file in *_orf_concat.fa; do grep -v ">" ${file} |tr '\n' ' ' | sed -e 's/ //g' > ${file/_orf_concat.fa/_final_concat.fa}; done
#NOTE: manually add newline after header: 
#for file in *_final_concat.fa; 
#do sed  '1 i\
#\>fastaseq ' ${file} > ${file%.*}.fasta;
#done
#rm *_final_concat.fa
#for file in *_final_concat.fasta; do sed -i "s/fastaseq/${file}/g" ${file}; done

#for cns_fa in *_geneious_cns.fasta; do sed -n '/^>/p' ${cns_fa} |sed -i 's/.*\ />/' ${cns_fa}; samtools faidx ${cns_fa}; bedtools getfasta -fi ${cns_fa} -bed PRV_CDS.bed -split -fo ${cns_fa/_cns.fasta/_orf_concat.fa}; done

#remove first n characters of every line
#for file in *_afr; do sed -r 's/.{23}//' ${file} > ${file/_afr/_af}; done

#remove everything before NC_...on lines starting with '>'
#sed -i 's/.*\ />/' ${cns_fa}

echo "-----------Re-aligning reads to respective sample consensus sequences...-----------"

for cns_fa in *.sort_final_cns.fa; 
do 
	bwa index ${cns_fa}; 
done

for filt_reads_for in *_R1_*_paired.fq; 
do 
	filt_reads_rev=${filt_reads_for/_R1_/_R2_}; 
	prefix=${filt_reads_for/_L001_R1_001_paired.fq/_}; 
	bwa mem -t 6 ${root_dir}/${prefix}.sort_final_cns.fa ${root_dir}/${filt_reads_for} ${root_dir}/${filt_reads_rev} | \
		samtools view -bS -F 0x04 - > ${root_dir}/${prefix}quasi.bam; 
done

for quasi_bam in *_quasi.bam; 
do 
	samtools sort -o ${quasi_bam%.*}.sort.bam ${quasi_bam}; 
	samtools index -b ${quasi_bam%.*}.sort.bam; 
done 

#Generate vcf using lofreq (sensitive sub-majority variant caller)
for fai_fa in *.sort_final_cns.fa; 
do 
	samtools faidx ${fai_fa}; 
done


echo -e "-----------Calling SNPs...-----------"
#default filters
# for sort_quasi_bam in *_quasi.sort.bam;
# do
# 	prefix=${sort_quasi_bam/_quasi.sort.bam/_};
# 	lofreq call -f ${root_dir}/${prefix}.sort_final_cns.fa -o ${sort_quasi_bam%.*}_rawlofreq.vcf -q 30 -Q 30 ${sort_quasi_bam} --verbose;
# 	bcftools filter -i '(500) < (DP4[0]+DP4[1]+DP4[2]+DP4[3]) & (AF) >= (0.01)' ${sort_quasi_bam%.*}_rawlofreq.vcf >${prefix}_lofreq.vcf;
# 	#doesn't work properly: lofreq filter -c 100 -a 0.01 -i ${sort_quasi_bam%.*}_rawlofreq.vcf -o ${sort_quasi_bam%.*}_lofreq.vcf;
# done

##delete later
for sort_quasi_bam in *_quasi.sort.bam;
do prefix=${sort_quasi_bam%-*};
lofreq call -f ${root_dir}/${prefix}-_geneious_.sort_final_cns.fa -o ${sort_quasi_bam%.*}_raw30lofreq.vcf -q 30 -Q 30 ${sort_quasi_bam} --verbose;
#option to choose read and variant allele frequency minimum:
#bcftools filter -i '(100) < (DP4[0]+DP4[1]+DP4[2]+DP4[3]) & (AF) >= (0.01)' ${sort_quasi_bam%.*}_raw30lofreq.vcf >${prefix}_100_lofreq.vcf;
bcftools filter -i '(500) < (DP4[0]+DP4[1]+DP4[2]+DP4[3]) & (AF) >= (0.01)' ${sort_quasi_bam%.*}_raw30lofreq.vcf >${prefix}_500_lofreq.vcf
done
##

cp *_quasi.bam ${ref_fasta} ./shannon_div/

for file in *_quasi.bam;
do mv ${file} ${file%.*}_.bam;
done


#if u don't choose to use default filters, run without first then apply user-specified filters in following command..use as pair
#lofreq call -f PRVrefgen.fasta -o sample_lofreq_nodefaults.vcf sample.sort.bam --no-default-filter #no defaults
#lofreq filter -v 1000 -a 0.01 -b bonf -q bonf -k bonf --verbose -i sample_lofreq_nodefaults.vcf -o sample_lofreq_nodefaults_fdrfilter.vcf #user-specified filters with fdr correction #type

#extract allele frequencies:
for lofreq_vcf in *_lofreq.vcf; 
do 
	prefix=${lofreq_vcf%.*}; 
	bcftools query -f '%CHROM %POS %AF\n' ${lofreq_vcf} 2>/dev/null >${prefix}_af;
done

# for lofreq_vcf in *_100_lofreq.vcf;
# do
# 	prefix=${lofreq_vcf%.*};
# 	bcftools query -f '%CHROM %POS %AF\n' ${lofreq_vcf} 2>/dev/null >${prefix}_af;
# done

#copy to result file - or move? 
mv *_af ./alt_allele_freq


echo -e "See files ending in _af for Alternate Allele Frequency per site."


name=`cut -f1 ${gtf}`

#split bed file - SNPGenie can only process one segment/ chrm at a time. NOTE: ensure ORF start is position -1. getfasta program chops off first number (inclusive)
# for chrm in ${name};
# do touch ${chrm}_.bed;
# 	while read -r line;
# 	do
# 		if echo "${line}" | grep -q "${chrm}";
# 		then echo "${line}" > ${chrm}_.bed;
# 		fi;
# 	done < ${gtf%.*}.bed;
# done

#split gtf file by segment
for chrm in ${name}; 
do 
	touch ${chrm}_snpgen.gtf; 
	while read -r line; 
	do 
		if echo "${line}" | grep -q "${chrm}"; 
		then echo "${line}" > ${chrm}_snpgen.gtf; 
		fi; 
	done < ${gtf}; 
done

#rm ${gtf%.*}.bed or name other bed files differently so as not to use the original bed file when splitting fastas..

#split reference fasta file by segment 

# for fa in *.sort_final_cns.fa;
# do
# 	prefix=${fa%-*};
# 	for chrm in ${name};
# 	do
# 		bedtools getfasta -fi ${data_dir}/${fa} -bed ${chrm}_.bed -fo ${prefix}_${chrm}_snpgen.fa;
# 	done;
# done

for fa in *.sort_final_cns.fa; 
do 
	prefix=${fa%-*}; 
	awk -v samp="$prefix" '/^>/ {OUT=samp "_" substr($0,2,11)  "_snpgen.fa"}; OUT {print > OUT}' ${fa}; 
done 

	
#split vcf by segment for SNPgenie analysis
#for file in *_lofreq.vcf; do bgzip ${file}; tabix -p vcf ${file}.gz; for bed in *_.bed; do tabix ${file}.gz -fhR ${bed} > ${file%%_*}_${bed%.*}snpgen.vcf; done; done

for file in *_lofreq.vcf;
do
	gzip ${file};
	zgrep "^#" ${file}.gz | gzip -c > header.vcf.gz;
	for chrm in ${name};
	do
		zgrep "^${chrm}" ${file}.gz | gzip -c > tmp.vcf.gz;
		zcat header.vcf.gz tmp.vcf.gz | gzip -c > ${file/_quasi.sort_lofreq.vcf/_}${chrm}_snpgen.vcf.gz;
	done
done


cp *_snpgen* ${snp_gen_dir} 

cd ${snp_gen_dir}
gunzip *.gz

echo -e "-----------Calculating Pi Diversity per genomic position...-----------"
#Run SNPgenie to calculate nucleotide diversity & whether mutations snps result in syn or nonsyn mutations 

for snpreport in *_snpgen.vcf;
do
	prefix=${snpreport%-*};
	for chrm in ${name}; 
	do 
		if [[ $snpreport == *"${chrm}"* ]] && [[ $snpreport == *"${prefix}"* ]]; 
		then 
			./snpgenie_dir_edits.pl --minfreq=0.01 --vcfformat=2 --snpreport=${snpreport} --fastafile=${prefix}_${chrm}_snpgen.fa --gtffile=${chrm}_snpgen.gtf; 
		# else
# 			echo "SNPGenie unsuccessful";
		fi; 
	done; 
done   

#OR1102 files with different naming:
# for snpreport in *OR1102*_snpgen.vcf;
# do
# 	prefix=${snpreport%_100*};
# 	for chrm in ${name};
# 		do
# 			if [[ $snpreport == *"${chrm}"* ]] && [[ $snpreport == *"${prefix}"* ]];
# 				then ./snpgenie_dir_edits.pl --minfreq=0.01 --vcfformat=2 --snpreport=${snpreport} --fastafile=${prefix}_.sort_final_cns.fa_${chrm}_snpgen.fa --gtffile=${chrm}_snpgen.gtf;
# 			fi;
# 	done;
# done


#copy files from various directories into one for data processing
mkdir all_pi_results; 

for chrm in ${name}; 
do 
	for dir in ${snp_gen_dir}/*_${chrm}_snpgen.vcf_SNPGenie_Results; 
	do 
		prefix=${dir##*/}; 
		cp ${dir}/site_results.txt ${dir%/*}/all_pi_results/${prefix/_snpgen.vcf_SNPGenie_Results/_}site_results.txt; 
	done; 
done

############for bicistronic segment..need to run for other orf on same segment, but will overwrite existing..so will rename first file.
##this will need to be edited depending on name of your polycistronic segment##################
for bicistronic in ${snp_gen_dir}/all_pi_results/*NC_036472.1*;
do
	mv ${bicistronic} ${bicistronic/NC_036472.1/NC_036472.1a}
done

rm -r *NC_036472.1_snpgen.vcf_SNPGenie_Results

#rerun with gtf of second orf
for snpreport in *NC_036472.1_snpgen.vcf;
do
	prefix=${snpreport%-*};
	if [[ $snpreport == *"${prefix}"* ]]; 
	then 
		./snpgenie_dir_edits.pl --minfreq=0.01 --vcfformat=2 --snpreport=${snpreport} --fastafile=${prefix}_NC_036472.1_snpgen.fa --gtffile=NC_036472.1b_snpgen.gtf; 
	fi; 
done  

#copy the site_results.txt files into all_pi_results with distinct name from first orf 
for dir in ${snp_gen_dir}/*_NC_036472.1_snpgen.vcf_SNPGenie_Results; 
do 
	prefix=${dir##*/}; 
	cp ${dir}/site_results.txt ${dir%/*}/all_pi_results/${prefix/_snpgen.vcf_SNPGenie_Results/b_}site_results.txt; 
done; 

cp ./all_pi_results/*site_results.txt ${root_dir}/nucleotide_div/ 

rm *_snpgen.vcf
rm -r *_SNPGenie_Results

cd ${root_dir}/nucleotide_div

#extract necessary columns 
for txt in *site_results.txt; 
do 
	prefix=${txt/_site_results.txt/_ntdiv}; 
	cut -f2,3,10,12 ${txt} >${prefix}.txt;
done

echo -e "See files ending in _ntdiv for Nucleotide (pi) Diversity per site."

cp *_ntdiv.txt ${root_dir}

cd ${root_dir}/shannon_div

echo -e "-----------Calculating Shannon Diversity per genomic position...-----------\a"

read -p "Enter minimum base quality [30]:" bq
read -p "Enter minimum mapping quality [20]:" mq
read -p "Enter minimum variant depth [0.01 = 1%]:" vd
read -p "Enter min for/rev read depths to call variant [50]:" sb
read -p "Enter length of longest segment [3934]:" bp 
bq=${bq:-30}
mq=${mq:-20}
vd=${vd:-0.01}
sb=${sb:-50}
bp=${bp:-3934}
#echo $bq

python2.7 PRdiversity.fromBAM.samtools1_9.py ${bq} ${mq} ${vd} ${sb} 

for file in *.mpileup; 
do 
	python2.7 mpileup2snp.dict.py $file ${vd} ${sb} ${bp} > ${file/mpileup/fixedmpileup.snp.csv}; 
done

for shandiv in *.fixedmpileup.snp.csv;
do 
	awk -F, '($3 < 500) {$18=0}1' OFS=, $shandiv | cut -d "," -f1,2,18 > ${shandiv/.fixedmpileup.snp.csv/_shandiv500.csv};
done


echo -e "See files ending in _shandiv500.csv for Shannon Diversity per site."

cp *_shandiv.csv ${root_dir}



