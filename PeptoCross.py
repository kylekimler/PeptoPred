import time
import joblib
import sys
from sklearn import svm
import sklearn
import sklearn.externals
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree
from sklearn.model_selection import KFold

call = time.time()

with open('./FullSetClassExtraction', 'rb') as d:
	Classifications = pickle.load(d)

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import matthews_corrcoef
for q in range(5, 17, 2):
	crosstime = time.time()
	print("cross validation of RFC window size %s currently under way." %int(q))
	with open('./FullSetFeatureExtraction%s' %int(q), 'rb') as f:
		Features = pickle.load(f)	
	clf = joblib.load('FullModelRFC%s.pkl' %int(q))
	x = cross_val_predict(clf, Features, Classifications, cv = 5, verbose = True, n_jobs=-1)
	print(cross_val_score(clf, Features, Classifications, cv = 5, verbose = True, n_jobs=-1))
	print(x)
	print(matthews_corrcoef(Classifications, x))
	print("cross validation of RFC window took %s seconds" %(time.time()-crosstime))