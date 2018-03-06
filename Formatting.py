import time
import joblib
import sys
from sklearn import svm
import sklearn
import sklearn.externals

call = time.time()

def ParseSeqs(filename):
	inputfasta=open(filename,'r')
	stringmaker = [str(line.rstrip('\n')) for line in inputfasta]
	for z in range(0,len(stringmaker)-1,3):
		stringmaker[z]=stringmaker[z].lstrip('>')
	parsedseqs = {stringmaker[x]: stringmaker[x+1] for x in range(0,len(stringmaker)-1,3)}
	inputfasta.close()
	# print(parsedseqs)
	return parsedseqs

def ParseClasses(filename):
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

def FormatSeqs(filename, windowsize):
	parsedseqs = ParseSeqs(filename)
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
		# classified = [SpIndex[attribute] for attribute in parsedclasses[key]]
		# print(parsedclasses[key])
		# classifications.extend(classified)
	# clf = svm.SVC()
	# clf.fit(features,classifications)
	return features

def FormatClasses(filename):
	parsedclasses = ParseClasses(filename)
	classifications = []
	for key in parsedclasses:
		classified = [SpIndex[attribute] for attribute in parsedclasses[key]]
		classifications.extend(classified)
	return classifications



	# Will probably need to parse PSSM file somehow - let's do that with another function - and into a list.

# small prediction test! do a list(clf.predict(from testpeps2.txt))
# maybe we do need to do this inside the for loop - with an if statement ? 
Features = FormatSeqs('globular_signal_peptide_2state.3line.txt',5)
Classifications = FormatClasses('globular_signal_peptide_2state.3line.txt')

'''with open('testpepsFeatures.pkl', 'wb') as out:
	joblib.dump(Features, out)
with open('testpepsClasses.pkl', 'wb') as ot:
	joblib.dump(Classifications, ot)'''

clf=svm.SVC()
clf.fit(Features, Classifications)

with open('testpepsModel.pkl', 'wb') as oot:
	joblib.dump(clf, oot)

SpIndexR = {0:'S',1:'G'}

FeaturesPred = FormatSeqs('testpeps2.txt',5)
prediction = list(clf.predict(FeaturesPred))
gg=''
for g in prediction:
	gg += str(g) 
print(gg)
with open('testpepsPred.txt', 'w+') as o:
	o.write(str(ParseSeqs('testpeps2.txt'))+'\n')
	o.write(str(gg))


'''with open('testpepsPred.txt', 'w+') as owt:
	owt.write('Prediction!' + '\n')
	owt.write(prediction)'''


print ("\nProcessing took %0.2f seconds.\n" %(time.time()-call))

