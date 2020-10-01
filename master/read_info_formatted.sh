#!/usr/bin/env bash 


echo -e "\n	 -----------------------
	      ALIGNMENT INFO
	 -----------------------\n"

#set -e 
#read -p "Enter directory containing paired fastq and sorted bam files:" working_dir
#echo -e "\n"
read -p "Enter minimum read depth to search ORF sites for:" rd
echo -e "\a\n"

#cd ${root_dir}

echo -e "\nCalculating read information for samples...This may take some time...\n"


let c=0; 
touch reads_metadata.txt pos_insuff_reads.txt; 
echo -e "sample\t#sites<${rd}reads\ttotal_reads_for\tavg_read_length_for\ttotal_reads_rev\tavg_read_length_rev\ttotal_mapped_reads\tproportion_mapped\tavg_cov" \
	>> reads_metadata.txt; 

for file in *.sort.bam; 
do 
	prefix=${file%_*}; 
	echo -e "\n\r${file}" >> reads_metadata.txt; 
	let c=0; 
	while IFS=$'\t' read -r line; 
	do 
		if [[ ${line##*$'\t'} -lt ${rd} ]]; 
		then 
			let c++; 
			echo -e "${file}\t${line}" >> pos_insuff_reads.txt; 
		else 
			let c=c+0; 
		fi; 
	done < ${prefix}_.sort.cov; 
	echo -e "${c}" >> reads_metadata.txt; 
		for paired in ${prefix}*_paired.fq; 
		do 
			echo $(cat ${paired} | wc -l)/4|bc >> reads_metadata.txt; 
			echo `awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' ${paired}` >> reads_metadata.txt;
		done;
echo `samtools view -c ${file}` >> reads_metadata.txt; 
done 

	
#IFS (internal field separator) can be set to default space: IFS=' ' or tab IFS=$'\t' or semicolon IFS=';'
#for cov in *.sort.cov; do let c=0; while IFS=$'\t' read -r line; do if [[ ${line##*$'\t'} -lt ${rd} ]]; then let c++; else let c=c+0; fi; done < ${cov}; echo -e "\n\r${cov} \n\r${c}" >> read_metadata.txt; done


awk '{ for (i=1;i<=NF;i++ ) printf $i "\t" }' reads_metadata.txt > reads_info.txt

#useless because even if program fails..command is technically still successful.. have to code at each step..
if [[ $? -gt 0 ]]; 
then 
	echo -e "\nProgram failed."; 
else echo -e "\nSee 'reads_info.txt' for results.";
fi

exit 0

























