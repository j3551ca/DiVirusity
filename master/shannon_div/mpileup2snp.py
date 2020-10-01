#@author: Vince Montoya
import sys
import re
import string
from decimal import Decimal, getcontext
import decimal
import math


def percentage(varDep,totalDepth):

	if float(varDep) > 0:
		nucPer = float(varDep)/float(totalDepth)
		freq = ((nucPer)*(math.log(nucPer)))
		return freq
	else:
		return 0
def shannonDiv(A,C,G,T,refdepth,ref,depth):
	percentages = []
	shannonPos = []
	percentages.append(abs(percentage(A,depth)))
	percentages.append(abs(percentage(C,depth)))
	percentages.append(abs(percentage(G,depth)))
	percentages.append(abs(percentage(T,depth)))
	percentages.append(abs(percentage(refdepth,depth)))
	shannonPos = sum(float(i) for i in percentages)
	return shannonPos
class pileup(object):
	def __init__(self,line):
		line = line.replace("^>","").replace("$","").replace("^[","").replace("^]","")
		linesp = line.split("\t")
		lineR = linesp[4].replace(",-1","").replace(",+1","").replace(".-1","").replace(".+1","")
		
		self.refGenome = linesp[0]
		self.pos = linesp[1]
		self.depth=int(linesp[3])
		self.ref=linesp[2]
		self.Refdepth=[]
		self.RefdepthFor=[]
		self.aster=[ (i.start()) for i in re.finditer('[*]', linesp[4])]
		self.Del=[ (i.start()) for i in re.finditer('[-]', linesp[4])]
		self.Ins=[ (i.start()) for i in re.finditer('[+]', linesp[4])]
		self.RefInd=[ (i.start()) for i in re.finditer('[.]', lineR)]
		self.RefIndR=[ (i.start()) for i in re.finditer('[,]', lineR)]
		self.RefD = int(len(self.RefInd))+int(len(self.RefIndR))
		if self.RefD==0:
			self.ref=''
		self.RefdepthFor=int(len(self.RefInd))
		self.RefdepthRev=int(len(self.RefIndR))
		self.AInd=[ (i.start()) for i in re.finditer('A', linesp[4])]
		self.CInd=[ (i.start()) for i in re.finditer('C', linesp[4])]
		self.GInd=[ (i.start()) for i in re.finditer('G', linesp[4])]
		self.TInd=[ (i.start()) for i in re.finditer('T', linesp[4])]
		self.aInd=[ (i.start()) for i in re.finditer('a', linesp[4])]
		self.cInd=[ (i.start()) for i in re.finditer('c', linesp[4])]
		self.gInd=[ (i.start()) for i in re.finditer('g', linesp[4])]
		self.tInd=[ (i.start()) for i in re.finditer('t', linesp[4])]
		self.depthG = len(self.GInd)+ len(self.gInd)
		self.depthC = len(self.CInd)+len(self.cInd)
		self.depthA = len(self.AInd)+len(self.aInd)
		self.depthT = len(self.TInd)+len(self.tInd)
		self.flagRef2 = 0
		# print self.RefdepthFor
		# print self.RefdepthRev
		if (self.RefdepthFor > 0) and (self.RefdepthRev > 0):
			self.flagRef2 = ' '
		else:
			self.flagRef2 = 'Ref Strand Bias'

	def variantSummary(self):
		
		indels = [self.Del, self.aster, self.Ins]
		variants  = ['A', 'C', 'G', 'T']
		#variants = [['A',self.AInd, self.aInd, self.depthA], ['C', self.CInd, self.cInd, self.depthC], ['G', self.GInd, self.gInd, self.depthG], ['T',self.TInd, self.tInd, self.depthT]]
		for i in indels:
			#print self.depth
			if (float(len(i))/float(self.depth)) >= float(cutoff):
				if i == 'aster':
					self.prevDel = str(len(self.aster))

				elif i == 'Del':
					self.Del = str(len(self.Del))

				elif i == 'Ins':
					self.Ins = str(len(self.Ins))
		if (float(len(self.aster))/float(self.depth)) >= float(cutoff):
			self.aster = str(len(self.aster))
		else:
			self.aster = 0
		if (float(len(self.Del))/float(self.depth)) >= float(cutoff):
			self.Del = str(len(self.Del))
		else:
			self.Del = 0
		if (float(len(self.Ins))/float(self.depth)) >= float(cutoff):
			self.Ins = str(len(self.Ins))
		else:
			self.Ins = 0

		varA = {'Depth': str(self.depthA), 'for': str(len(self.AInd)), 'rev': str(len(self.aInd))}
		varC = {'Depth': str(self.depthC), 'for': str(len(self.CInd)), 'rev': str(len(self.cInd))}
		varG = {'Depth': str(self.depthG), 'for': str(len(self.GInd)), 'rev': str(len(self.gInd))}
		varT = {'Depth': str(self.depthT), 'for': str(len(self.TInd)), 'rev': str(len(self.tInd))}
		if float(self.depthA)/float(self.depth) >= float(cutoff):
			self.varA = varA
		else:
			self.varA = {'Depth': '0', 'for': '0', 'rev': '0'}
		if float(self.depthC)/float(self.depth) >= float(cutoff):
			self.varC = varC
		else:
			self.varC = {'Depth': '0', 'for': '0', 'rev': '0'}
		if float(self.depthG)/float(self.depth) >= float(cutoff):
			self.varG = varG
		else:
			self.varG = {'Depth': '0', 'for': '0', 'rev': '0'}
		if float(self.depthT)/float(self.depth) >= float(cutoff):
			self.varT = varT
		else:
			self.varT = {'Depth': '0', 'for': '0', 'rev': '0'}


		
#This program takes a -q 20 restriction cutoff samtools pileup, and outputs the reference seq and the variants that meet the required depths.  It also flags those variants that have only forward or rev reads (flag=1)
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



f = open(mpileupFile,"r")
filehandle = f.readlines()

print "Pos,RefGenome,Total Depth,Reference,DepthFor,DepthRev,Ref Strand Bias,Total A,For A,Rev A,Total C,For C,Rev C,Total G,For G,Rev G,Total T,For T,Rev T"
count = 0
prevcount = 0
for line in filehandle:
	pileupSum = pileup(line)
	if pileupSum.depth > 0:
		pileupSum.variantSummary()
		print str(pileupSum.pos)+','+pileupSum.refGenome+','+str(pileupSum.depth)+','+pileupSum.ref+','+str(pileupSum.RefdepthFor)+','+str(pileupSum.RefdepthRev)+','+str(pileupSum.flagRef2)+','+pileupSum.varA['Depth']+','+pileupSum.varA['for']+','+pileupSum.varA['rev']+','+pileupSum.varC['Depth']+','+pileupSum.varC['for']+','+pileupSum.varC['rev'] +','+pileupSum.varG['Depth']+','+pileupSum.varG['for']+','+pileupSum.varG['rev']+','+pileupSum.varT['Depth']+','+pileupSum.varT['for']+','+pileupSum.varT['rev']+','+str(shannonDiv(pileupSum.varA['Depth'],pileupSum.varC['Depth'],pileupSum.varG['Depth'],pileupSum.varT['Depth'],pileupSum.RefD, pileupSum.ref, pileupSum.depth))
