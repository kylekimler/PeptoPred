import time
import joblib
import sys
from sklearn import svm
import sklearn
import sklearn.externals
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree
from Formatting import Format, FormatClasses, ParseClasses, ParseSeqs, ParseSeqstoDict
import os.path

call = time.time()

print("Creating Model! Do not Disturb")
Features = []
Classifications = []
filename = '../data/'+(os.path.basename(sys.argv[1]))

#overarching if statement checking to see if there are argv's. otherwise just run as windowsize 15 and on the full set in the data


if os.path.isfile('../Extrct/%sClassExtraction' %os.path.basename(sys.argv[1]).rstrip('.txt')):
	with open('../Extrct/%sClassExtraction'  %os.path.basename(sys.argv[1]).rstrip('.txt'), 'rb') as d:
		Classifications = joblib.load(d)
else:
	FormatClasses(str(filename))
	with open('../Extrct/%sClassExtraction'  %os.path.basename(sys.argv[1]).rstrip('.txt'), 'rb') as d:
		Classifications = joblib.load(d)

if os.path.isfile('../Extrct/%sFeatureExtraction' %(str(os.path.basename(sys.argv[1]).rstrip('.txt')))):
	with open('../Extrct/%sFeatureExtraction%s' %(str(os.path.basename(sys.argv[1]).rstrip('.txt')), int(sys.argv[2])), 'rb') as f:
		Features = joblib.load(f)
else:
	Format(str(filename), int(sys.argv[2]))
	with open('../Extrct/%sFeatureExtraction%s' %(str(os.path.basename(sys.argv[1]).rstrip('.txt')), int(sys.argv[2])), 'rb') as f:
		Features = joblib.load(f)

clf=svm.LinearSVC()
model = clf.fit(Features, Classifications)

with open('../Models/FullModelLin%s.pkl' %int(sys.argv[2]), 'wb') as oot:
	pickle.dump(model, oot)

modtime = time.time()
print("\nModeling took %0.2f seconds.\n" %(modtime-call))