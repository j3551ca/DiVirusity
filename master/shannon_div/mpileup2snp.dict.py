#@author:Vince Montoya
import sys
import re
import string
from decimal import Decimal, getcontext
import decimal
import math

try:
	mpileupFile = sys.argv[1]
except:
	print 'Specify input file'
	sys.exit(1)
try:
	cutoff = sys.argv[2]
except:
	print 'Input a depth threshold for variant calling (ex: 0.1 etc where 0.1 = 10 percent depth)'
	print "If not important, enter '0.0'"
	sys.exit(1)
try:
	biasFlag = sys.argv[3]
except:
	print 'Input a seq threshold if strand bias is desired (100, where at least a 100 for and rev reads are required)'
	print "If not important, enter '0'"
	sys.exit(1)
try:
	refLen = sys.argv[4]
except:
	print 'Input a reference length, if multiple segments are present, use the longest reference length'
	sys.exit(1)



f = open(mpileupFile,"r")
filehandle = f.readlines()

print "RefGenome,Pos,Total Depth,Adep,Afor,Arev,Cdep,Cfor,Crev,Gdep,Gfor,Grev,Tdep,Tfor,Trev,Deletions,Insertions,ShannonDiv"
count = 0
prevcount = 0
pileupSum = {}
def percentage(varDep,totalDepth):

	if float(varDep) > 0:
		nucPer = float(varDep)/float(totalDepth)
		freq = ((nucPer)*(math.log(nucPer)))

		return freq
	else:
		return 0
def shannonDiv(nucDepths,refdepth,depth):
	percentages = []
	shannonPos = []
	for n in nucDepths:
		percentages.append(abs(percentage(n,depth)))
		# percentages.append(abs(percentage(C,depth)))
		# percentages.append(abs(percentage(G,depth)))
		# percentages.append(abs(percentage(T,depth)))
	#percentages.append(abs(percentage(refdepth,depth)))
	#print percentages
	shannonPos = sum(float(i) for i in percentages)
	return shannonPos
def pileupParse(refI,listNucs,depth):
	depthSum = {}
	count = 0
	for r in refI:
		refD =  str(sum(refI[1:]))
		refFor =  str(refI[1])
		refRev =  str(refI[2])
	for l in listNucs:
		count +=1
		if count == 1:
			Adep =  str(sum(l))
			Afor = str(As[0])
			Arev = str(As[1])
			if 'A' in refI[0]:

				depthSum['A'] = [refD,refFor,refRev]
			else:	
				depthSum['A'] = [Adep,Afor,Arev]
		elif count == 2:
			Cdep =  str(sum(Cs))
			Cfor =  str(Cs[0])
			Crev =  str(Cs[1])
			if 'C' in refI[0]:
				depthSum['C'] = [refD,refFor,refRev]
			else:
				depthSum['C'] = [Cdep,Cfor,Crev]
		elif count == 3:
			Gdep =  str(sum(Gs))
			Gfor =  str(Gs[0])
			Grev =  str(Gs[1])
			if 'G' in refI[0]:
				depthSum['G'] = [refD,refFor,refRev]
			else:
				depthSum['G'] = [Gdep,Gfor,Grev]
		elif count == 4:
			Tdep =  str(sum(Ts))
			Tfor =  str(Ts[0])
			Trev =  str(Ts[1])
			if 'T' in refI[0]:
				depthSum['T'] = [refD,refFor,refRev]
			else:
				depthSum['T'] = [Tdep,Tfor,Trev]
	nucs = ('A','C','G','T')
	depths = []
	for n in nucs:
		depths.append(depthSum[n][0])
		shannon = str(shannonDiv(depths,refD,depth))
	depthSum['div'] = shannon
	return depthSum



for line in filehandle:
	line = line.replace("^>","").replace("$","").replace("^[","").replace("^]","")
	linesp = line.split("\t")
	lineR = linesp[4].replace(",-1","").replace(",+1","").replace(".-1","").replace(".+1","")
	refGenome = linesp[0]
	pos = linesp[1]
	refInfo = []
	if not refGenome in pileupSum:
		pileupSum[refGenome] = {}
		positions = {}
		for x in range(1,int(refLen)+1):
			positions[str(x)] = {}
			# print positions
		pileupSum[refGenome] = positions

	pileupSum[refGenome][pos] = [0,{},[],[],[],[],0,0]
	#Depth,RefD,A,C,G,T,Del,Ins
	totalDep = int(linesp[3])
	pileupSum[refGenome][pos][0] = totalDep
	if totalDep >= float(cutoff):
		

		ref=linesp[2]
		Refdepth=[]
		RefdepthFor=[]
		aster=[ (i.start()) for i in re.finditer('[*]', linesp[4])]
		Del=[ (i.start()) for i in re.finditer('[-]', linesp[4])]
		Ins=[ (i.start()) for i in re.finditer('[+]', linesp[4])]
		RefInd=[ (i.start()) for i in re.finditer('[.]', lineR)]
		RefIndR=[ (i.start()) for i in re.finditer('[,]', lineR)]
		RefD = int(len(RefInd))+int(len(RefIndR))
		if RefD==0:
			ref=''
		RefdepthFor=int(len(RefInd))
		RefdepthRev=int(len(RefIndR))
		AInd=[ (i.start()) for i in re.finditer('A', linesp[4])]
		CInd=[ (i.start()) for i in re.finditer('C', linesp[4])]
		GInd=[ (i.start()) for i in re.finditer('G', linesp[4])]
		TInd=[ (i.start()) for i in re.finditer('T', linesp[4])]
		aInd=[ (i.start()) for i in re.finditer('a', linesp[4])]
		cInd=[ (i.start()) for i in re.finditer('c', linesp[4])]
		gInd=[ (i.start()) for i in re.finditer('g', linesp[4])]
		tInd=[ (i.start()) for i in re.finditer('t', linesp[4])]
		depthG = len(GInd)+ len(gInd)
		depthC = len(CInd)+len(cInd)
		depthA = len(AInd)+len(aInd)
		depthT = len(TInd)+len(tInd)
		refInfo = [ref,RefdepthFor,RefdepthRev]
		pileupSum[refGenome][pos][1] = refInfo
		pileupSum[refGenome][pos][2] = [len(AInd),len(aInd)]
		pileupSum[refGenome][pos][3] = [len(CInd),len(cInd)]
		pileupSum[refGenome][pos][4] = [len(GInd),len(gInd)]
		pileupSum[refGenome][pos][5] = [len(TInd),len(tInd)]
		pileupSum[refGenome][pos][6] = len(Del)
		pileupSum[refGenome][pos][7] = len(Ins)

for r in pileupSum:
	for p in pileupSum[r]:
		# print r+'\t'+str(p)+'\t'+str(pileupSum[r][p])
		if len(pileupSum[r][p]) > 0:
			# print r+'\t'+str(p)+'\t'+str(pileupSum[r][p])
			depth = pileupSum[r][p][0]
			if float(depth) >= float(cutoff):
				refI = pileupSum[r][p][1]
				refD = sum(refI[1:])
				As = pileupSum[r][p][2]
				Cs = pileupSum[r][p][3]
				Gs = pileupSum[r][p][4]
				Ts = pileupSum[r][p][5]
				Dels = pileupSum[r][p][6]
				Ins = pileupSum[r][p][7]
				depthInfo = pileupParse(refI,[As,Cs,Gs,Ts],depth)
				#shannonPos = shannonDiv(depthInfo['A'][0],depthInfo['C'][0],depthInfo['G'][0],depthInfo['C'][0],refD,depth)
				#print r+'\t'+str(p)+'\t'+str(shannonPos)
				# print r+'\t'+str(p)
				# print depthInfo['C'][0]
				# print refD
				# print depth
				# print depthInfo['div']
				print r+','+str(p)+','+str(depth)+','+','.join(depthInfo['A'])+','+','.join(depthInfo['C'])+','+','.join(depthInfo['G'])+','+','.join(depthInfo['T'])+','+str(Dels)+','+str(Ins)+','+str(depthInfo['div'])
			else:
				print r+','+str(p)+','+','.join(['0']*16)


		else:
			print r+','+str(p)+','+','.join(['0']*16)
	# flagRef2 = 0
	# # print RefdepthFor
	# # print RefdepthRev
