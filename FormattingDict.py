import time
import os.path
import sys
from sklearn import svm
import sklearn
import sklearn.externals

call = time.time()

def parse_fasta(filename):
	inputfasta=open(filename,'r')
	stringmaker = [str(line.rstrip('\n')) for line in inputfasta]
	for z in range(0,len(stringmaker)-1,3):
		stringmaker[z]=stringmaker[z].lstrip('>')
	parsedseqs = {stringmaker[x]: stringmaker[x+1] for x in range(0,len(stringmaker)-1,3)}
	inputfasta.close()
	# print(parsedseqs)
	return parsedseqs

def parse_classes(filename):
	inputfasta=open(filename,'r')
	stringmaker = [str(line.rstrip('\n')) for line in inputfasta]
	for z in range(0,len(stringmaker),3):
		stringmaker[z]=stringmaker[z].lstrip('>')
	parsedclasses = {stringmaker[x]: stringmaker[x+2] for x in range(0,len(stringmaker),3)}
	inputfasta.close()
	# print(parsedclasses)
	return parsedclasses

parsedseqdict = parse_fasta('testpeps.txt')
parsedclassdict = parse_classes('testpeps.txt')

aaIndex = {'A':(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
'C':(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
'D':(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
'E':(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
'F':(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
'G':(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
'H':(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
'I':(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
'K':(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
'L':(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
'M':(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
'N':(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
'P':(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
'Q':(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
'R':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
'S':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
'T':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
'V':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
'W':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
'Y':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
'?':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
'X':(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05),
'B':(0,0,0.5,0,0,0,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0),
'Z':(0,0,0,0.5,0,0,0,0,0,0,0,0,0,0.5,0,0,0,0,0,0),
'J':(0,0,0,0,0,0,0,0.5,0,0.5,0,0,0,0,0,0,0,0,0,0)}

SpIndex = {'S':0,'G':1}

def FormatSeqs(filename, windowsize):
	parsedseqs = ParseSeqstoDict(filename)
	#print(parsedseqs)
	#print(parsedclasses)
	prots = []
	features = []
	classifications = []
	for key in parsedseqs:

		#print(parsedseqs[key])
		for x in range(0,len(parsedseqs[key])):
			if x<windowsize//2:
				flist=[]
				flist.extend(aaIndex['?']*(windowsize//2-x))
				z=0
				while z <= x+windowsize//2:
					flist.extend(aaIndex[parsedseqs[key][z]])
					z+=1
				features.append(flist)
			elif len(parsedseqs[key])-windowsize//2 > x >= (windowsize//2-1):
				middlewindow = []
				middlewindow.extend(zz for attribute in parsedseqs[key][x-(windowsize//2):x+(windowsize//2)+1] for zz in aaIndex[attribute])
				features.append(middlewindow)
			elif(x>=(len(parsedseqs[key])-windowsize//2)-1):
				elist=[]
				z = x - (windowsize//2)
				while z <= len(parsedseqs[key]) - 1:
					elist.extend(aaIndex[parsedseqs[key][z]])
					z+=1
				elist.extend(aaIndex['?'] * (x+windowsize//2-(len(parsedseqs[key])-1)))
				features.append(elist)
	return features

################ old dict-only sliding window maker

joblib.dump(FormatAndFit(parsedseqdict,parsedclassdict,5), 'test.pkl')
print ("\nDumped the model in %0.2f seconds.\n" %(time.time()-call))


print(FormatAndFit(parsedseqdict,parsedclassdict,5))
