import time
import joblib
import sys
from sklearn import svm
import sklearn
import sklearn.externals
import pickle
from Formatting import Format, FormatClasses, ParseClasses, ParseSeqs, ParseSeqstoDict

call = time.time()

clf = joblib.load('FullModelDTC15.pkl')

SpIndexR = {0:'S',1:'G'}



PepstoPred = sys.argv[1]

FeaturesPred = Format(PepstoPred, 15)
FeaturesPredOut = ParseSeqstoDict(PepstoPred)
StatsComp = FormatClasses(PepstoPred)

prediction = list(clf.predict(FeaturesPred))
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

with open('%sDTCPred15.txt' %sys.argv[1].rstrip('.txt'), 'w+') as o:
	for key in FeaturesPredOut:
		o.write('>' + key + '\n')
		o.write(FeaturesPredOut[key] + '\n')
		pos = pos + len(FeaturesPredOut[key])
		aa = [SpIndexR[prediction[gg]] for gg in range(pos-len(FeaturesPredOut[key]),pos)]
		for ggg in aa:
			o.write(ggg)
		o.write('\n')

print("\nPercentage of False S in this prediction: %0.9f percent. \n" %(float(100*FalseS/TotalS)))
print("\nPercentage of missed S in this prediction: %0.9f percent. \n" %(float(100*MissedS/TotalS)))
print("\nTotal accuracy of this prediction is: %0.9f percent. Congratulations\n" %(float(100*(x-wrong)/x)))
print ("\nPrediction took %0.2f seconds.\n" %(time.time()-call))
