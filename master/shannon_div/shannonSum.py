#@author:Vince Montoya
import sys
mpileupOut = sys.argv[1]
refseqs = sys.argv[2]
shortName = {
	'NC_036476.1': 'coreShell', 
	'NC_036468.1':'coreTurret',
	'NC_036477.1':'RdRp',
	'NC_036469.1':'NTPase',
	'NC_036470.1': 'outerShell',
	'NC_036471.1':'nonstructuralProt',
	'NC_036472.1':'outerClamp',
	'NC_036473.1':'coreClamp',
	'NC_036474.1':'nonstructuralRNA',
	'NC_036475.1':'outerFiber'
}
seqLengths = {}
with open(refseqs) as f1:
	for line in f1:
		if ">" in line:
			refid = line.split(' ')[0][1:]
		else:
			seqLengths[refid] = len(line.rstrip())
shannonTotals = {}
with open(mpileupOut) as f2:
	for line in f2:
		if not 'RefGenome' in line:
			linesp = line.rstrip().split(',')
			ref = linesp[1]
			shannon = linesp[19]
			shannonTotals.setdefault(ref,[]).append(float(shannon))
#ShannonDivSum is the shannon diversity for a given gene, summed per position
#ShannonDivSumNorm is more comparable for gene comparisons as diversity here is normalized per gene length
print 'Reference,Gene,ShannonDivSum,ShannonDivSumNorm'
for s in shannonTotals:
	try:
		gene = shortName[s]
	except:
		gene = s
	print s+","+gene+','+str(sum(shannonTotals[s]))+','+str(sum(shannonTotals[s])/int(seqLengths[s]))

