import time
import os.path
import sys
from sklearn import svm
import sklearn
import sklearn.externals

def parse_fasta(filename):
    inputfasta=open(filename,'r')
    stringmaker = [str(line[:len(line)-1]) for line in inputfasta]
    for z in range(0,len(stringmaker)-1,3):
        stringmaker[z]=stringmaker[z].lstrip('>')
    parsedseqs = {stringmaker[x]: stringmaker[x+1] for x in range(0,len(stringmaker)-1,3)}
    inputfasta.close()
    #print(parsedseqs)
    return parsedseqs

def parse_classes(filename):
    inputfasta=open(filename,'r')
    stringmaker = [str(line[:len(line)-1]) for line in inputfasta]
    for z in range(0,len(stringmaker)-1,3):
        stringmaker[z]=stringmaker[z].lstrip('>')
    parsedclasses = {stringmaker[x]: stringmaker[x+2] for x in range(0,len(stringmaker)-1,3)}
    inputfasta.close()
    #print(parsedclasses)
    return parsedclasses

parsedseqdict = parse_fasta('testpeps.txt')
parsedclassdict = parse_classes('testpeps.txt')

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

'''aaIndex2 = {'A':1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7,
            'I':8, 'K':9, 'L':10, 'M':11, 'N':12, 'P':13, 'Q':14,
            'R':15, 'S':16, 'T':17, 'V':18, 'W':19, 'Y':20, '?':0}'''


#sklearn.feature_extraction.DictVectorizer(parse_fasta('singletest1.fasta'))
#FE = sklearn.preprocessing.OneHotEncoder(n_values=21)
#print(FE)

def FormatAndFit(parsedseqs, parsedclasses, windowsize):
    protslist = list(parsedseqs.values())
    classeslist = list(parsedclasses.values())
    print(classeslist)
    #indexlist = [aaIndex(z) for z in str(parse_fasta('x').values())]
    #print(indexlist)
    features = []
    classifications = []
    for prot in protslist:
        count = 0
        for x in range(0,len(prot)):
            if x<windowsize//2:
                flist=[]
                flist.extend(aaIndex['?']*(windowsize//2-x))
                z=0
                while z <= x+windowsize//2:
                    flist.extend(aaIndex[prot[z]])
                    z+=1
                features.append(flist)
            elif len(prot)-windowsize//2 > x >= (windowsize//2-1):
                middlewindow = []
                middlewindow.extend(zz for attribute in prot[x-(windowsize//2):x+(windowsize//2)+1] for zz in aaIndex[attribute])
                features.append(middlewindow)
            elif(x>=(len(prot)-windowsize//2)-1):
                elist=[]
                z = x - (windowsize//2)
                while z <= len(prot) - 1:
                    elist.extend(aaIndex[prot[z]])
                    z+=1
                elist.extend(aaIndex['?'] * (x+windowsize//2-(len(prot)-1)))
                features.append(elist)
        #classifications.append(xx for attribute in parsedclasses(parsedseqs.key(prot)) for xx in SpIndex[attribute])
        print(parsedclasses(parsedseqs[prot]))
        #print(classifications)
        count+=1
    clf = svm.SVC()
    clf.fit(features,classifications)
    return clf

print(FormatAndFit(parsedseqdict,parsedclassdict,5))
