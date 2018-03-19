import time
from sklearn.externals import joblib
import sys
from sklearn import svm
import sklearn
import sklearn.externals
import pickle
import os.path

call = time.time()

def ParseSeqstoDict(filename):
	inputfasta=open(filename,'r')
	stringmaker = [str(line.rstrip('\n')) for line in inputfasta]
	for z in range(0,len(stringmaker)-1,3):
		stringmaker[z]=stringmaker[z].lstrip('>')
	parsedseqs = {stringmaker[x]: stringmaker[x+1] for x in range(0,len(stringmaker)-1,3)}
	inputfasta.close()
	return parsedseqs

def ParseSeqs(filename):
	inputfasta=open(filename,'r')
	lines = inputfasta.readlines()
	z = str(lines[3][0])
	inputfasta.close()
	inputfasta = open(filename, 'r')
	if(z=='>'):
		stringmaker = [str(line.rstrip('\n').lstrip('>')) for line in inputfasta]
		parsedseqs = [stringmaker[x+1] for x in range(0,len(stringmaker)-1,3)]
	elif(z=='>'):
		stringmaker = [str(line.rstrip('\n').lstrip('>')) for line in inputfasta]
		parsedseqs = [stringmaker[x+1] for x in range(0,len(stringmaker)-1,2)]
	inputfasta.close()
	return parsedseqs

def ParseClasses(filename):
	inputfasta=open(filename,'r')
	stringmaker = [str(line.rstrip('\n').lstrip('>')) for line in inputfasta]
	parsedseqs = [stringmaker[x+2] for x in range(0,len(stringmaker)-1,3)]
	inputfasta.close()
	return parsedseqs

def ParseClassestoDict(filename):
	inputfasta=open(filename,'r')
	stringmaker = [str(line.rstrip('\n')) for line in inputfasta]
	for z in range(0,len(stringmaker),3):
		stringmaker[z]=stringmaker[z].lstrip('>')
	parsedclasses = {stringmaker[x]: stringmaker[x+2] for x in range(0,len(stringmaker),3)}
	inputfasta.close()
	return parsedclasses

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

def Format(filename, windowsize):
	features = []
	parsedseqs = ParseSeqs(filename)
	
	if type(parsedseqs[0]) is str:
		parsed2 = []
		for item in parsedseqs:
			seq = []
			for attribute in item:
				seq.append(aaIndex[attribute])
			parsed2.append(seq)
	
	for item in parsed2:
		prot = item
		count = 0
		for aa in prot:
			if count<windowsize//2:
				flist=[]
				flist.extend(aaIndex['?']*(windowsize//2-count))
				z=0
				while z <= count+windowsize//2:
					flist.extend(prot[z])
					z+=1
				features.append(flist)
			elif len(prot)-windowsize//2 > count >= (windowsize//2)-1:
				middlewindow = []
				for zz in prot[count-(windowsize//2):count+(windowsize//2)+1]:
					middlewindow.extend(zz)
				features.append(middlewindow)
			elif(count>=(len(prot)-windowsize//2)-1):
				elist=[]
				z = count - (windowsize//2)
				while z <= len(prot) - 1:
					elist.extend(prot[z])
					z+=1
				elist.extend(aaIndex['?'] * (count+windowsize//2-(len(prot)-1)))
				features.append(elist)
			count+=1
	with open('../Extrct/%sFeatureExtraction%s' %(os.path.basename(filename).rstrip('.txt'),windowsize), 'wb') as f:
		joblib.dump(features, f, compress = 9)
	pass

def FormatClasses(filename):
	parsedclasses = ParseClassestoDict(filename)
	classifications = []
	for key in parsedclasses:
		classified = [SpIndex[attribute] for attribute in parsedclasses[key]]
		classifications.extend(classified)
	with open('../Extrct/%sClassExtraction' %(os.path.basename(filename).rstrip('.txt')), 'wb') as f:
		joblib.dump(classifications, f, compress = 9)
	pass


#if statement checking to see if files exist
#Features = Format('globular_signal_peptide_2state.3line.txt', int(sys.argv[1]))
# Classifications = FormatClasses('globular_signal_peptide_2state.3line.txt')

#with open('FullSetFeatureExtraction%s' %int(sys.argv[1]), 'wb') as f:
#	pickle.dump(Features, f)
'''with open('FullSetClassExtraction', 'wb') as d:
	pickle.dump(Classifications, d)'''
if('%0.2f' %(time.time()-call) != '0.00'):
	print ("\nFormatting took %0.2f seconds.\n" %(time.time()-call))

