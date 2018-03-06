def FormatSansDicts(PSSMarray, PSSMclasses, windowsize):
	features = []
	classifications = []
	for PSSM in PSSMarray:
		count = 0
		for aa in PSSM:
			if aa<windowsize//2:
				flist=[]
				flist.extend(aaIndex['?']*(windowsize//2-aa))
				z=0
				while z <= aa+windowsize//2:
					flist.extend(PSSM[z])
					z+=1
				features.append(flist)
			elif len(parsedseqs[key])-windowsize//2 > aa >= (windowsize//2-1):
				middlewindow = []
				middlewindow.extend(zz for zz in PSSM[aa-(windowsize//2):aa+(windowsize//2)+1])
				features.append(middlewindow)
			elif(aa>=(len(PSSM)-windowsize//2)-1):
				elist=[]
				z = aa - (windowsize//2)
				while z <= len(PSSM) - 1:
					elist.extend(PSSM[z])
					z+=1
				elist.extend(aaIndex['?'] * (aa+windowsize//2-(len(PSSM)-1)))
				features.append(elist)
		if(PSSMclasses[0] is str):
			classified = [SpIndex[attribute] for attribute in PSSMclasses]
			classifications.extend(classified)
		count+=1