import time
import joblib
import sys
from sklearn import svm
import sklearn
import sklearn.externals
import pickle
from Formatting import Format, FormatClasses, ParseClasses, ParseSeqs, ParseSeqstoDict
import re
import os.path

call = time.time()

model = '../Models/globular_signal_peptide_2state.3lineDTC15.pkl'
ws = int(re.findall(r'\d+', model)[2])
with open('../Models/globular_signal_peptide_2state.3lineDTC15.pkl', "rb") as f:
	clf = pickle.load(f)

SpIndexR = {0:'S',1:'G'}

PepstoPred = sys.argv[1]
print("formatting fasta features")
FeaturesPred = Format(PepstoPred, ws) 
FeaturesPredOnTest = ParseSeqstoDict(PepstoPred)
print("format complete. Predicting peptides")
prediction = list(clf.predict(FeaturesPred))

#StatsComp = FormatClasses(PepstoPred) 
gg = ''
x = 0
FalseS = 0
FalseG = 0
wrong = 0
TotalS = 0
MissedS = 0

#The StatsComp expression was for internal testing - you know the classification of your fasta. 

for g in prediction:
	gg += str(g) 
	'''if(g == 0):
		TotalS += 1
	if(g == 0 and StatsComp[x] == 1):
		FalseS += 1
	if(g == 1 and StatsComp[x] == 0):
		MissedS += 1
	if(g != StatsComp[x]):
		wrong += 1
	x += 1'''

pos = 0
with open('../results/%sLin%sPred.txt' %(os.path.basename(sys.argv[1]).rstrip('.txt'), ws), 'w+') as o: 
	for key in FeaturesPredOnTest:
		o.write('>' + key + '\n')
		o.write(FeaturesPredOnTest[key] + '\n')
		pos = pos + len(FeaturesPredOnTest[key])
		aa = [SpIndexR[prediction[gg]] for gg in range(pos-len(FeaturesPredOnTest[key]),pos)]
		for ggg in aa:
			o.write(ggg)
		o.write('\n')
'''
print("\nTotal S predicted: %i." %int(TotalS))
print("\nFalse S in this prediction: %0.1f" %(FalseS))
print("\nMissed S in this prediction: %0.1f" %(MissedS))
print("\nTotal accuracy of this prediction is: %0.9f percent. Congratulations" %(float(100*(x-wrong)/x)))'''
print("\nPrediction took %0.2f seconds.\n" %(time.time()-call))
