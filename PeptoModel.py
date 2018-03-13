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

clf=tree.DecisionTreeClassifier()
model = clf.fit(Features, Classifications)

with open('FullModelDTC%s.pkl' %int(sys.argv[1]), 'wb') as oot:
	joblib.dump(model, oot)

modtime = time.time()-call
print("\nModeling took %0.2f seconds.\n" %(modtime))

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import matthews_corrcoef

print("Cross Validating, Do Not Distorb!1!!1")

x = cross_val_predict(model, Features, Classifications, cv = 5, verbose = True)

print(cross_val_score(model, Features, Classifications, cv = 5, verbose = True))
print(x)
print(matthews_corrcoef(Classifications, x))

print("\nCross-Validation took %0.2f seconds.\n" %(time.time()-modtime))