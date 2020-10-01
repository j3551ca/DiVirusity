
import os
import sys
import datetime
current = datetime.datetime.now()
TimeStamp = current.strftime("%Y-%m-%d-%H.%M")

################BEFORE USING
#@author: Vince Montoya
#This pipeline is optional, you can either ignore it and just run the diversity scripts (mpileup2snp.py and shannonSum.py)...
#Or you can modify it so that it works on your computer
#You need BWA, Samtools, and python (ver2) installed
#Here is what you need to modify to ensure (**hopefully**) it works on your computer (SEE BELOW FOR CONTEXT)
#Update dirPath and suffix ('.fastq') to only include the forward reads
#Ensure that the reverse reads are identical to the suffix except for the _R1 flag
#Specify filename (including path) of reference file
#Specify the number of reads required by both strands to call a variant in the mpileup file 
#Sometimes strand bias is a problem, if not important just enter 0 (default)
################

#######Usage
#python PRV.diversity.py BQ MQ variantDepth strandbiasDepth
#BQ = base quality (good starting point 20)
#MQ = mapping quality (if there are conserved regions in your segments this flag can cause problems, so a lower value is better here. Otherwise, 20 is a good starting point.
#variantDepth is the depth required to call a variant (0.05 = 5%), enter 0 if u want everything
#strandbiasDepth is the number of sequences required for both forward and reverse strands to call a variant, enter 0 if u want everything
#######

try:
	BQ = sys.argv[1]
except:
	print 'Specify base quality score for mpileup gen...'
	sys.exit(1)
try:
	MQ = sys.argv[2]
except:
	print 'Specify mapping quality score for mpileup gen...'
	sys.exit(1)
try:
	variantDepthCutoff = sys.argv[3]
except:
	print 'Specify the depth required to call a variant in the mpileup file (ex: 0.05 = 5%)'
	print "If not important, enter '0'"
	sys.exit(1)
try:
	strandBiasCutoff = sys.argv[4]
except:
	print 'Specify the number of reads required by both strands to call a variant in the mpileup file (ex 5 = minimum 5 reads)'
	sys.exit(1)

work_dir_path = (os.path.dirname(os.path.abspath(__file__)))
print(work_dir_path)
print("Hello, World!")



def unwrapFasta(Ref):
	cmd = 'python unwrapFasta.py %s' %(Ref)
	logfile.write("Unwrapping ref fasta file:\n"+ cmd+"\n")
	os.system(cmd)
	return

def index(infile):
	cmd = 'bwa index %s' %(infile)
	logfile.write("BWA mapping with:\n"+ cmd+"\n")
	os.system(cmd)
	return

def run_BWA(Ref,SAM,forward, reverse):
	cmd = "bwa mem -t 8 -M %s %s %s >%s" % (Ref, forward, reverse, SAM)
	logfile.write("BWA mapping with:\n"+ cmd+"\n")
	os.system(cmd)
	return

def run_BAM(Ref, SAM, BAMsort, BAMindex):
	cmd1 = "samtools view -bS %s | samtools sort -o %s -" % (SAM, BAMsort)
	cmd2 = "samtools index -b %s >%s" %(BAMsort, BAMindex)
	logfile.write("Generating BAM with command:\n"+ cmd1+"\n")
	os.system(cmd1)
	logfile.write("Generating BAM with command:\n"+ cmd2+"\n")
	os.system(cmd2)
	
	return
def run_BAMsi(BAM, BAMsort, BAMindex):
	cmd1 = 'samtools sort -o %s %s' %(BAMsort, BAM)
	cmd2 = "samtools index -b %s %s" %(BAMsort, BAMindex)
	logfile.write("Sorting BAM file with command:\n"+ cmd1+"\n")
	os.system(cmd1)
	logfile.write("Generating BAM index with command:\n"+ cmd2+"\n")
	os.system(cmd2)
	return

def run_mpileup(Ref,BAMsort,BQ,MQ, MPILEUP):
	cmd1 = "samtools mpileup -d 100000 -q %s -Q %s -f %s %s > %s" % (BQ,MQ,Ref, BAMsort,MPILEUP)
	logfile.write("Generating MPILEUP with command:\n"+ cmd1+"\n")
	os.system(cmd1)
	return

def diversity(MPILEUP, variantDepthCutoff, strandBiasCutoff,refmod):
	shannonDiv = MPILEUP+".shannon.csv"
	cmd1 = "python mpileup2snp.py %s %s %s > %s" %(MPILEUP,variantDepthCutoff, strandBiasCutoff,mpileupSNP)
	cmd2 = "python shannonSum.py %s %s> %s" %(mpileupSNP, refmod, shannonDiv)
	logfile.write("Calculating SNPs from MPILEUP:\n"+ cmd1+"\n")
	os.system(cmd1)
	logfile.write("Summing shannon diversity for each gene:\n"+ cmd2+"\n")
	os.system(cmd2)
	return
#Logfile
if os.path.isfile('Mpileup2DiversityLogFile.txt'):
	logfile = open('Mpileup2DiversityLogFile.txt','a')
	logfile.write("\n")
	logfile.write("Analysis initiated at:"+ TimeStamp+"\n")
else:
	logfile = open('Mpileup2DiversityLogFile.txt','w')
	logfile.write("\n")
	logfile.write("Analysis initiated at:"+ TimeStamp+"\n")

####File parsing
#Update dirPath and suffix ('.fastq') to only include the forward reads
#Currently dirPath is set to look for fastq files in the current directory
#Ensure that the reverse reads are identical to the suffix except for the _R1 flag
filenames = [f for f in os.listdir(work_dir_path) if "_.bam" in f]
if len(filenames) == 0:
	print "Could not locate specified file(s)..."

#Specify filename (including path) of reference file

Ref = '/data2/jessica/programs/ihDiv/shannon_div/PRVrefgen.fasta'
#str(input('Enter absolute path to reference genome:'))
#'/data/jessica/programs/JessQuasDiv/shannon_div/PRVrefgen.fasta'
#refmod = Ref.split(".fasta")[0]+".mod.fasta"
if not os.path.isfile(Ref):
	unwrapFasta(Ref)
if not os.path.isfile(Ref+'.amb'):
	index(Ref)
#Looping and processing
for bam in filenames:
	fileSp = str(bam).split(".bam")
	BAMsort = fileSp[0]+".sort.bam"
	BAMindex = fileSp[0]+".sort.bai"
	MPILEUP = fileSp[0]+"."+str(BQ)+"."+str(MQ)+".mpileup"
	mpileupSNP = MPILEUP+".snp.csv"
	#print bam
	if not os.path.isfile(BAMsort):
		print "Sorting and indexing BAM file..."
		run_BAMsi(bam, BAMsort, BAMindex)
	if not os.path.isfile(MPILEUP):
		print "Generating mpileup..."
		run_mpileup(Ref,BAMsort,BQ,MQ,MPILEUP)
	if not os.path.isfile(mpileupSNP):
		print "Calculating diversity..."
		
		diversity(MPILEUP, variantDepthCutoff, strandBiasCutoff, Ref)









current2 = datetime.datetime.now()
TimeStampFin = current2.strftime("%Y-%m-%d-%H:%M")
logfile.write("Analysis completed at:"+ TimeStampFin+"\n")


