import time
import joblib
import sys
from sklearn import svm
import sklearn
import sklearn.externals
import pickle
from Formatting import Format, FormatClasses, ParseClasses, ParseSeqs, ParseSeqstoDict
import re

call = time.time()

model = '../Models/FullModelRFC5.pkl'
ws = int(re.findall(r'\d+', model)[0])
print(ws)
clf = pickle.load(model)

SpIndexR = {0:'S',1:'G'}



PepstoPred = sys.argv[1]

FeaturesPred = Format(PepstoPred, ws) 
FeaturesPredOut = ParseSeqstoDict(PepstoPred)

prediction = list(clf.predict(FeaturesPred))


StatsComp = FormatClasses(PepstoPred) #change to depend on 3 or 2 line
gg = ''
x = 0
FalseS = 0
FalseG = 0
wrong = 0
TotalS = 0
MissedS = 0
for g in prediction:
	gg += str(g) 
	if(StatsComp[x]==0):
		TotalS += 1
	if(g == 0 and StatsComp[x] == 1):
		FalseS += 1
	if(g == 1 and StatsComp[x] == 0):
		MissedS += 1
	if(g != StatsComp[x]):
		wrong += 1
	x += 1

pos = 0

#change window size to depend on clf
with open('../results/%sRFC%sPred.txt' %(sys.argv[1].rstrip('.txt'),ws), 'w+') as o: 
	for key in FeaturesPredOut:
		o.write('>' + key + '\n')
		o.write(FeaturesPredOut[key] + '\n')
		pos = pos + len(FeaturesPredOut[key])
		aa = [SpIndexR[prediction[gg]] for gg in range(pos-len(FeaturesPredOut[key]),pos)]
		for ggg in aa:
			o.write(ggg)
		o.write('\n')

#need to write an if statement to make StatsComp dependent
print("\nTotal S in predicted dataset: %i. \n" %int(TotalS))
print("\nPercentage of False S in this prediction: %0.9f percent. \n" %(float(100*FalseS/TotalS)))
print("\nPercentage of missed S in this prediction: %0.9f percent. \n" %(float(100*MissedS/TotalS)))
print("\nTotal accuracy of this prediction is: %0.9f percent. Congratulations\n" %(float(100*(x-wrong)/x)))
print("\nPrediction took %0.2f seconds.\n" %(time.time()-call))
