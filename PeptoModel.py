import time
import joblib
import sys
from sklearn import svm
import sklearn
import sklearn.externals
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree

call = time.time()

print("Creating Model! Do not Disturb")

with open('./FullSetFeatureExtraction%s' %int(sys.argv[1]), 'rb') as f:
	Features = pickle.load(f)

with open('./FullSetClassExtraction', 'rb') as d:
	Classifications = pickle.load(d)

clf=sklearn.ensemble.RandomForestClassifier()
model = clf.fit(Features, Classifications)

with open('FullModelRFC%s.pkl' %int(sys.argv[1]), 'wb') as oot:
	joblib.dump(model, oot)

modtime = time.time()
print("\nModeling took %0.2f seconds.\n" %(modtime-call))