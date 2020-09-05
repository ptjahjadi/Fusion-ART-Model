from ARTfunc import *
import copy
import json

FRACTION = 0.00001

class FusionART:
	def __init__(self,numspace=0,lengths=[],beta=[],alpha=[],gamma=[],rho=[],schema={},numarray=True):
		self.codes=[]
		self.icode = {'F2':0, 'weights':[]}
		self.alpha = list(alpha)
		self.gamma = list(gamma)
		self.rho = list(rho)
		self.pmcriteria = [1.0]*numspace

		if type(beta) is float:
			self.beta = beta
		if type(beta) is list:
			self.beta = list(beta)
		
		self.lastChoice = []
		self.lastMatch = []
		self.lastActRho = []
		
		#Default functions for choice, match (per field), weight update (every field), resonance search and match (all fields).
		#The functions can be programmatically updated or replaced at runtime.
		#------------------------------------------------------------------------ 
		self.choiceAct = choiceActFuzzy
		self.compMatch = matchFuncFuzzy
		self.updWeight = []
		self.resonance = resonanceFuzzy
		self.matchVal = matchValFuzzy
		#------------------------------------------------------------------------ 
		
		self.numarray = numarray
		
		
		
		if len(schema) <= 0:
			self.icode = self.initCode(numspace,lengths)
			if len(self.icode['weights']) > 0:
				self.codes.append(self.icode)
			self.schema = {}
		else:
			self.activityF1 = []
			self.schema = copy.deepcopy(schema)
			self.F1FieldsSetup(self.schema)
		
			
		if len(self.updWeight) < len(self.activityF1):
			for k in range(len(self.activityF1)):
				if self.numarray:
					self.updWeight.append(updWeightsFuzzy)
				else:
					self.updWeight.append(updWeightFuzzy)
		
		
		self.prevF2Sel = 0
		self.prevUncommit = True
	
	#def initCode(self,nspace,lengths,icode={'F2':0, 'weights':[]}):
	def initCode(self,nspace,lengths,ivalue=0.0, wvalue=1.0):
			iicode = {'F2':0, 'weights':[]}
			self.activityF1=[[]]*nspace #old
			if len(lengths) >= nspace:
				wght = []
				for k in range(len(self.activityF1)):
					#self.activityF1[k] = [1.0]*lengths[k]
					self.activityF1[k] = [ivalue]*lengths[k]
					wght.append([wvalue]*lengths[k])
				iicode['weights'] = list(wght)
			return copy.deepcopy(iicode)

			
	def F1FieldsSetup(self, fschemas):
		actTmp = []
		schmTmp = []
		if len(self.activityF1) <= 0:
			lengths = []
			for i in range(len(fschemas)):
				if 'compl' in fschemas[i]:
					if fschemas[i]['compl']:
						lengths.append(len(fschemas[i]['attrib'])*2)
					else:
						lengths.append(len(fschemas[i]['attrib']))
				else:
					lengths.append(len(fschemas[i]['attrib']))
			self.icode = self.initCode(len(fschemas),lengths)
			if len(self.icode['weights']) > 0:
				self.codes.append(self.icode)
			if len(self.alpha) < len(fschemas):
				self.alpha = [1.0]*len(fschemas)
			if len(self.beta) < len(fschemas):
				self.beta = [1.0]*len(fschemas)
			if len(self.gamma) < len(fschemas):
				self.gamma = [1.0]*len(fschemas)
			if len(self.rho) < len(fschemas):
				self.rho = [1.0]*len(fschemas)
			self.pmcriteria = [1.0]*len(fschemas)
			if len(self.updWeight) < len(self.activityF1):			
				for k in range(len(self.activityF1)):
					if self.numarray:
						self.updWeight.append(updWeightsFuzzy)
					else:
						self.updWeight.append(updWeightFuzzy)
		for k in range(len(fschemas)):
			fschema = initFieldSchema(fschemas[k])
			schmTmp.append(fschema)
			factivity = getActivityFromField(fschema)
			self.setActivityF1(factivity,kidx=k)
		self.F1Fields = schmTmp
		
	def buttUpAllF1(self):
		if hasattr(self,'F1Fields'):
			for k in range(len(self.F1Fields)):
				self.setActivityF1(getActivityFromField(self.F1Fields[k]),kidx=k)
	

	def updateF1bySchema(self,fschemas,refresh=True):
		for k in range(len(fschemas)):
			if 'name' in fschemas[k]:
				for kf in self.F1Fields:
					if isSchemaWithAtt(kf,'name',fschemas[k]['name']):
						kf.update(fschemas[k])
						if refresh:
							kf.update(refreshComplSchema(kf))
		self.buttUpAllF1()

	def updateF1byAttVal(self,attvals,kidx=-1,name='',refresh=True):
		if kidx >= 0:
			self.F1Fields[kidx].update(setSchemabyAttVal(self.F1Fields[kidx],attvals))
		if len(name) > 0:
			for kf in self.F1Fields:
				if isSchemaWithAtt(kf,'name',name):
					kf.update(setSchemabyAttVal(kf,attvals))
					if refresh:
						kf.update(refreshComplSchema(kf))
		self.buttUpAllF1()
		
	def updateF1byVals(self,vals,kidx=-1,name='',refresh=True):
		if kidx >= 0:
			self.F1Fields[kidx].update(setValFieldSchema(self.F1Fields[kidx],vals))
		if len(name) > 0:
			for kf in self.F1Fields:
				if isSchemaWithAtt(kf,'name',name):
					kf.update(setValFieldSchema(kf,vals))
					if refresh:
						kf.update(refreshComplSchema(kf))
		self.buttUpAllF1()

	
	def buttUpF1(self,fschema,kidx=-1,fname=''):
		if kidx >= 0:
			self.F1Fields[kidx].update(fschema)
		else:
			if len(fname)>0:
				for k in range(len(self.F1Fields)):
					if 'name' in self.F1Fields[k]:
						if self.F1Fields[k]['name'] == fname:
							self.buttUpF1(fschema,k)
							break
			else:
				for k in range(len(fschema)):
					self.buttUpF1(fschema[k],k)
		for k in range(len(self.F1Fields)):
			self.setActivityF1(getActivityFromField(self.F1Fields[k]),kidx=k)

	def TopDownF1(self):
		F1f = []
		if (len(self.activityF1) > 0) and (len(self.F1Fields)>0):
			for k in range(len(self.activityF1)):
				c = False
				if 'compl' in self.F1Fields[k]:
					c = self.F1Fields[k]['compl']
				self.F1Fields[k].update(readOutVectSym(self.activityF1[k], c))
			F1f = copy.deepcopy(self.F1Fields)
		return F1f
	
	def setActivityF1(self,val,kidx=-1,iidx=-1):
		if kidx > -1:
			if iidx > -1:
				self.activityF1[kidx][iidx] = val
			else:
				self.activityF1[kidx] = list(val)
		else:
			self.activityF1 = list(val)
		
	
	def setActivityF2(self,val,jidx=-1):
		if jidx > -1:
			self.codes[jidx]['F2'] = val 
		else:
			assert (len(val) == len(self.codes))
			for j in range(len(val)):
				self.codes[j]['F2'] = val[j]

	def setParam(self,param,value,k=-1):
		if param == "beta":
			if type(value) is float:
				self.beta = value
			if type(value) is list:
				if k>0:
					self.beta[k]=value
				else:
					self.beta=list(value)
		if param == "alpha":
			if k>=0:
				self.alpha[k]=value
			else:
				self.alpha=list(value)
		if param == "gamma":
			if k>=0:
				self.gamma[k]=value
			else:
				self.gamma=list(value)
		if param == "rho":
			if k>=0:
				self.rho[k]=value
			else:
				self.rho=list(value)

	def compChoice(self):
		a = self.activityF1
		w = listAttVal(self.codes,'weights')
		if self.numarray:
			a = np.array(a)
			w = np.array(w)
		self.codes = attValList(list(self.choiceAct(a,w,self.alpha,self.gamma,listAttVal(self.codes,'F2'))),self.codes,'F2')
			
		
	def expandCode(self):
		tw = []
		for k in range(len(self.activityF1)):
			tw.append([1]*len(self.activityF1[k]))
		self.codes.append({'F2':0, 'weights':list(tw)})
		
	def uncommitted(self,idx):
		for k in range(len(self.codes[idx]['weights'])):
			if self.numarray:
				sumw = np.sum(self.codes[idx]['weights'][k])
			else:
				sumw = 0
				for i in range(len(self.codes[idx]['weights'][k])):
					sumw += self.codes[idx]['weights'][k][i]
			#if sumw <= (len(self.codes[idx]['weights'][k])/2):
			if sumw < len(self.codes[idx]['weights'][k]):
				return False
		return True
		
	def codeCompetition(self):
		maxact = -1
		c = -1
		if self.numarray:
			c = np.argmax(listAttVal(self.codes,'F2'))
		else:
			for j in range(len(self.codes)):
				if self.codes[j]['F2'] > maxact:
					maxact = self.codes[j]['F2']
					if maxact > 0:
						c = j
		return c
		
	def doLearn(self,j):
		for k in range(len(self.activityF1)):
			if self.numarray:
				self.codes[j]['weights'][k] = self.updWeight[k](self.beta[k], self.codes[j]['weights'][k], self.activityF1[k])
			else:
				for i in range(len(self.activityF1[k])):
					self.codes[j]['weights'][k][i] = self.updWeight[k](self.beta[k], self.codes[j]['weights'][k][i], self.activityF1[k][i])
				
				
	def doOverwrite(self,j):
		for k in range(len(self.activityF1)):
			self.codes[j]['weights'][k] = list(self.activityF1[k])
			
	def autoLearn(self,j,overwrite=False):
		if self.uncommitted(j):
			overwrite=True
			self.expandCode()
		if overwrite:
			self.doOverwrite(j)
		else:
			self.doLearn(j)
			
	def doReadout(self,j,k,overwrite=True):
		if overwrite:
			self.activityF1[k] = self.codes[j]['weights'][k]
		else:
			if self.numarray:
				self.activityF1[k] = np.amin([self.activityF1[k],self.codes[j]['weights'][k]],axis=0)
			else:
				for i in range(len(self.activityF1[k])):
					self.activityF1[k][i] = min(self.activity[k][i],self.codes[j]['weights'][k][i])
					
	def doReadoutAllFields(self, j, overwrite=True):
		for k in range(len(self.activityF1)):
			self.doReadout(j,k,overwrite)


	def isResonance(self, j, rhos=[]):
		crhos = list(self.rho)
		if len(rhos)>0:
			crhos = list(rhos)
		w = self.codes[j]['weights']
		if self.numarray:
			w = np.array(w)
		if self.resonance(self.activityF1,w,chros):
			return True
		else:
			return False
			
	def rhotracking(m,fraction):
		return min(m+fraction,1)
	
	def resSearch(self,mtrack=[],rhos=[],F2filter=[], duprep=False, prevSel=[]):
		resetcode = True
		J = -1
		crhos = self.rho
		if len(rhos)>0:
			crhos = list(rhos)
		self.lastActRho = list(crhos)
		self.compChoice()
		self.lastChoice = listAttVal(self.codes,'F2')
		
		while resetcode:
			resetcode = False
			J = self.codeCompetition()
			if J >= 0:
				matches = list(self.matchVal(self.activityF1,self.codes[J]['weights']))
				self.lastMatch = list(matches)
				if not duprep and pmismatch(matches,self.pmcriteria):
					return J
				if(not mresonance(matches,crhos)) or (J in F2filter):
					self.codes[J]['F2'] = 0
					resetcode = True
					for m in range(len(mtrack)):
						if crhos[mtrack[m]] < matches[mtrack[m]]:
							crhos[mtrack[m]] = self.rhotracking(matches[mtrack[m]],FRACTION)
				self.lastActRho = list(crhos)
				
				if duprep and not resetcode:
					if J in prevSel:
							self.codes[J]['F2'] = 0
							resetcode = True
		if J >= 0:
			self.prevF2sel = J
			self.prevUncommit = self.uncommitted(J)
		return J
					
	def displayNetwork(self):
		for j in range(len(self.codes)):
			print('Code: ' + str(j) + ' ' + str(self.codes[j]))
		print ('-----------------------------------------')
		print ('F1: ' + str(self.activityF1))
		
	def displayNetParam(self):
		print('alpha: ' + str(self.alpha))
		print('beta: '+ str(self.beta))
		print('gamma: ' + str(self.gamma))
		print('rho: ' + str(self.rho))

	def expandInput(self,idxs=[],quant=1,ivalue=0,wvalue=0,wvalue_uncommit=1):
		for q in range(quant):
			for i in range(len(idxs)):
				self.activityF1[idxs[i]].append(ivalue)
				#for j in range(len(self.weights)):
				for j in range(len(self.codes)):
					if self.uncommitted(j):
						self.codes[j]['weights'][idxs[i]].append(wvalue_uncommit)
					else:
						self.codes[j]['weights'][idxs[i]].append(wvalue)
				

	def expandInputCompl(self, kidx=-1, ivalue=0, wvalue=0, wvalue_uncommit=1):
		if kidx >= 0:
			if kidx < len(self.activityF1):
				if len(self.activityF1[kidx])%2 == 0:
					midx = int(len(self.activityF1[kidx])/2)
					self.activityF1[kidx].insert(midx,ivalue)
					self.activityF1[kidx].append(ivalue)
					for j in range(len(self.codes)):
						if self.uncommitted(j):
							self.codes[j]['weights'][kidx].insert(midx,wvalue_uncommit)
							self.codes[j]['weights'][kidx].append(wvalue_uncommit)
						else:
							self.codes[j]['weights'][kidx].insert(midx,wvalue)
							self.codes[j]['weights'][kidx].append(wvalue)
						
					

	def expandInputwSchema(self,kidx=-1,name='',ivalue=0,wvalue=0,attname=''):
		if kidx >= 0:
			ffield = self.F1Fields[kidx]
		if len(name) > 0:
			for ki in range(len(self.F1Fields)):
				if isSchemaWithAtt(self.F1Fields[ki],'name',name):
					ffield = self.F1Fields[ki]
					kidx = ki
		compl = False
		if 'compl' in ffield:
			compl = ffield['compl']
		
		if 'attrib' in ffield:
			if len(attname) > 0:
				ffield['attrib'].append(attname)
			else:
				att = "ax" + str(len(ffield['attrib'])+1)
				ffield['attrib'].append(att)
		if compl:	
			ffield['val'].append(ivalue)
			ffield['vcompl'].append(ivalue)
			self.expandInputCompl(kidx=kidx,ivalue=ivalue,wvalue=wvalue)
		else:
			ffield['val'].append(ivalue)
			self.expandInput([kidx],ivalue=ivalue,wvalue=wvalue)
		return ffield



					
	def removeInput(self,k=-1,idx=-1):
		if k >= 0:
			if k < len(self.activityF1):
				if idx >= 0:
					self.activityF1[k].remove(self.activityF1[k][idx])
					for j in range(len(self.codes)):
						self.codes[j]['weights'][k].remove(self.codes[j]['weights'][k][idx])
		return idx
						
	
	def removeInputwSchema(self, k=-1, name=''):
		if k >= 0:
			ffield = self.F1Fields[k]
		if len(name) > 0:
			for ki in range(len(self.F1Fields)):
				if isSchemaWithAtt(self.F1Fields[ki],'name',name):
					ffield = self.F1Fields[ki]
					kidx = ki
		compl = False
		
			
			
	def removeCode(self,idx=-1):
		if idx >= 0:
			if idx < len(self.codes):
				retcode = copy.deepcopy(self.codes[idx])
				self.codes.remove(self.codes[idx])
				return retcode
		return {}

	def getActivityF2(self):
		return listAttVal(self.codes,'F2')


	def gradEncActivateF1(self, k=-1, idx=-1, tau=0.1, tresh=0.0):
		if k>=0:
			if k < len(self.activityF1):
				self.setActivityF1(decayVals(self.activityF1[k], tau, tresh),kidx=k)
				if idx >= 0:
					if idx < len(self.activityF1[k]):
						self.setActivityF1(val=1-tresh,kidx=k,iidx=idx)
						
	def gradComplEnvActivateF1(self, k=-1, idx=-1, tau=0.1, tresh=0.0):
		if k>=0:
			return
						
						
	def stackTopFusionART(self, TopfusART):
		self.TopFusionART = TopfusART
		
	def linkF2TopF1BySchema(self, sNameList):
		if hasattr(self,'TopFusionART'):
			self.toplinkedSchemas = []
			for i in range(len(self.TopFusionART.F1Fields)):
				if self.TopFusionART.F1Fields[i]['name'] in sNameList:
					self.toplinkedSchemas.append(self.TopFusionART.F1Fields[i])
					
	def linkF2TopF1(self, F1idxList, cF1idxList):
		if hasattr(self, 'TopFusionART'):
			self.topIdxList = []
			self.topcIdxList = []
			for k in range(len(self.TopFusionART.activityF1)):
				if k in F1idxList:
					self.topIdxList.append(k)
				if k in cF1idxList:
					self.topcIdxList.append(k)
					
	
	def SequentialResSearch(self, mtrack=[],rhos=[],F2filter=[], tau=0.1, tresh=0.0, maxv=0.9, stau=0.01, accdig=10, duprep=True):
		J = self.resSearch(mtrack, rhos, F2filter)
		if hasattr(self, 'TopFusionART'):
			fTop = self.TopFusionART
			for i in range(len(self.topIdxList)):
				if duprep:
					pSel = [x for x in range(len(fTop.activityF1[self.topIdxList[i]])) if fTop.activityF1[self.topIdxList[i]][x] > 0]
					J = self.resSearch(mtrack, rhos, F2filter, duprep=True, prevSel=pSel)				
				if self.uncommitted(J) and (len(self.codes)>1):
					fTop.expandInput(idxs=[self.topIdxList[i]])
					if len(self.topcIdxList) > 0:
						fTop.expandInput(idxs=[self.topcIdxList[i]])
				if len(self.topcIdxList) > 0:
					v = fTop.activityF1[self.topIdxList[i]]
					cv = fTop.activityF1[self.topcIdxList[i]]
					print("v, cv before: ", v, cv)
					v, cv = insertComplDecayVals(v, cv, tau=tau, tresh=tresh, idx=J, maxv=maxv, stau=stau, accdig=accdig)
					fTop.setActivityF1(v,kidx=self.topIdxList[i])
					fTop.setActivityF1(cv,kidx=self.topcIdxList[i])
				else:
					v = fTop.activityF1[self.topIdxList[i]]
					v = insertDecayVals(v, tau=tau, tresh=tresh, idx=J, maxv=maxv, accdig=accdig)
					fTop.setActivityF1(v,kidx=self.topIdxList[i])
		return J						
	
	
					
	def SchemaBasedSequentialResSearch(self,mtrack=[],rhos=[],F2filter=[], tau=0.1, tresh=0.0, maxv=0.9, stau=0.01, accdig=10, duprep=True):
		J = self.resSearch(mtrack, rhos, F2filter)
		if hasattr(self, 'TopFusionART'):
			fTop = self.TopFusionART			
			for scm in self.toplinkedSchemas:
				if duprep:
					pSel = [i for i in range(len(scm['val'])) if scm['val'][i] > 0]
					J = self.resSearch(mtrack, rhos, F2filter, duprep=True, prevSel=pSel)
				if self.uncommitted(J) and (len(self.codes)>1):
					fTop.expandInputwSchema(name=scm['name'])
				if scm['compl']:
					v = scm['val']
					cv = scm['vcompl']
					v, cv = insertComplDecayVals(v, cv, tau=tau, tresh=tresh, idx=J, maxv=maxv, stau=stau, accdig=accdig)
					fTop.updateF1bySchema([{'name':scm['name'], 'val':v, 'vcompl':cv}],refresh=False)
				else:
					v = insertDecayVals(scm['val'], tau=tau, tresh=tresh, idx=J, maxv=maxv, accdig=accdig)
					fTop.updateF1bySchema([{'name':scm['name'], 'val':v}])					
		return J
		
						

	def seqTopReadoutToF1(self, tau=0.1, tresh=0.0, stau=0.01, accdig=10, queue=True, overwrite=True):
		if hasattr(self, 'TopFusionART'):
			fTop = self.TopFusionART
			J=-1
			for i in self.topIdxList:
				if len(self.topcIdxList) > 0:
					v = fTop.activityF1[self.topIdxList[i]]
					cv = fTop.activityF1[self.topcIdxList[i]]
					J, v, cv = maxComplReadoutVals(v, cv, tau=tau, tresh=tresh, stau=stau, accdig=accdig, queue=queue)
					fTop.setActivityF1(v,kidx=self.topIdxList[i])
					fTop.setActivityF1(cv,kidx=self.topcIdxList[i])
				else:
					v = fTop.activityF1[self.topIdxList[i]]
					J, v = maxReadoutVals(v, tau=tau, tresh=tresh, accdig=accdig, queue=queue)
					fTop.setActivityF1(v,kidx=self.topIdxList[i])
				self.doReadoutAllFields(J, overwrite=overwrite)
				self.TopDownF1()
			return J


	def seqTopReadoutToF1Schema(self, tau=0.1, tresh=0.0, stau=0.01, accdig=10, queue=True, overwrite=True):
		if hasattr(self, 'TopFusionART'):
			fTop = self.TopFusionART
			J=-1
			for scm in self.toplinkedSchemas: 
				if scm['compl']:
					v = scm['val']
					cv = scm['vcompl']
					J, v, cv = maxComplReadoutVals(v, cv, tau=tau, tresh=tresh, stau=stau, accdig=accdig, queue=queue)
					fTop.updateF1bySchema([{'name':scm['name'], 'val':v, 'vcompl':cv}],refresh=False)
				else:
					v = scm['val']
					J, v = maxReadoutVals(v, tau=tau, tresh=tresh, accdig=accdig, queue=queue)
					fTop.updateF1bySchema([{'name':scm['name'], 'val':v}])
				self.doReadoutAllFields(J, overwrite=overwrite)
				self.TopDownF1()
			return J



def saveFusionARTNetwork(nnet, name='fart.net'):
	fartnet = {'file_name': name, 
			'codes': nnet.codes,
			'alpha': nnet.alpha,
			'beta': nnet.beta,
			'gamma': nnet.gamma,
			'rho': nnet.rho,
			'lastChoice': nnet.lastChoice,
			'lastMatch': nnet.lastMatch,
			'lastActRho': nnet.lastActRho,
			'activityF1': nnet.activityF1,
			#'F1Fields': nnet.F1Fields
			}
	if hasattr(nnet, 'F1Fields'):
		fartnet['F1Fields'] = nnet.F1Fields
	with open(name, 'w') as outfile:
		json.dump(fartnet, outfile)

def loadFusionARTNetwork(nnet, name='fart.net'):
	with open(name) as json_file:
		fartnet = json.load(json_file)
	nnet.codes = fartnet['codes']
	nnet.alpha = fartnet['alpha']
	nnet.beta = fartnet['beta']
	nnet.gamma = fartnet['gamma']
	nnet.rho = fartnet['rho']
	nnet.lastChoice = fartnet['lastChoice']
	nnet.lastMatch = fartnet['lastMatch']
	nnet.lastActRho = fartnet['lastActRho']
	nnet.activityF1 = fartnet['activityF1']
	if 'F1Fields' in fartnet:
		nnet.F1Fields = fartnet['F1Fields']
			