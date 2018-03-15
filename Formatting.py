import time
from sklearn.externals import joblib
import sys
from sklearn import svm
import sklearn
import sklearn.externals
import pickle
from sklearn.model_selection import cross_val_score

call = time.time()

def ParseSeqstoDict(filename):
	inputfasta=open(filename,'r')
	stringmaker = [str(line.rstrip('\n')) for line in inputfasta]
	for z in range(0,len(stringmaker)-1,3):
		stringmaker[z]=stringmaker[z].lstrip('>')
	parsedseqs = {stringmaker[x]: stringmaker[x+1] for x in range(0,len(stringmaker)-1,3)}
	inputfasta.close()
	# print(parsedseqs)
	return parsedseqs

def ParseSeqs(filename):
	inputfasta=open(filename,'r')
	stringmaker = [str(line.rstrip('\n').lstrip('>')) for line in inputfasta]
	parsedseqs = [stringmaker[x+1] for x in range(0,len(stringmaker)-1,3)]
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
	# print(parsedclasses)
	return parsedclasses

# parsedseqdict = parse_fasta('testpeps.txt')
# parsedclassdict = parse_classes('testpeps.txt')

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

#sklearn.feature_extraction.DictVectorizer(parse_fasta('singletest1.fasta'))
#FE = sklearn.preprocessing.OneHotEncoder(n_values=21)
#print(FE)



def Format(filename, windowsize):
	features = []
	classifications = []
	parsedseqs = ParseSeqs(filename)
	
	if type(parsedseqs[0]) is str:
		parsed2 = []
		for item in parsedseqs:
			seq = []
			for attribute in item:
				seq.append(aaIndex[attribute])
			parsed2.append(seq)
	parsedclasses = ParseClasses(filename)
	
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
					# middlewindow.extend(zz for zz in prot[count-(windowsize//2):count+(windowsize//2)+1])
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
	return features

def FormatClasses(filename):
	parsedclasses = ParseClassestoDict(filename)
	classifications = []
	for key in parsedclasses:
		classified = [SpIndex[attribute] for attribute in parsedclasses[key]]
		classifications.extend(classified)
	return classifications

#Features = Format('globular_signal_peptide_2state.3line.txt', int(sys.argv[1]))
# Classifications = FormatClasses('globular_signal_peptide_2state.3line.txt')

#with open('FullSetFeatureExtraction%s' %int(sys.argv[1]), 'wb') as f:
#	pickle.dump(Features, f)
'''with open('FullSetClassExtraction', 'wb') as d:
	pickle.dump(Classifications, d)'''
'''
print(len(Features))
print(len(Classifications))'''
'''with open('testpepsFeatures.pkl', 'wb') as out:
	joblib.dump(Features, out)
with open('testpepsClasses.pkl', 'wb') as ot:
	joblib.dump(Classifications, ot)'''
# clf = joblib.load('FullModel.pkl')
'''clf=svm.LinearSVC(cache_size=1000)
clf.fit(Features, Classifications)'''




print ("\nFormatting took %0.2f seconds.\n" %(time.time()-call))

