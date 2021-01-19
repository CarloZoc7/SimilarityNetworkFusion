import pandas as pd

"""
	Class to create a table to show the approch use for the different omnic and
	to show the result obtained.
"""

class ResultTable:
	def __init__(self, name='', cluster='', usedBonf=False, usedPCA=False, usedMMS=False, usedSS=False, usedLT=False,
				usedSF=False, silhouette=0.0, RandIndex=0.0):
		self.name = name
		self.cluser = cluster
		self.usedPCA = usedPCA
		self.usedMMS = usedMMS
		self.usedSS = usedSS 
		self.usedLT = usedLT
		self.usedSF = usedSF
		self.usedBonf = usedBonf
		self.silhouette = silhouette
		self.randIndex = RandIndex
		
		self.row_for_table = {
				'Omnic Name': '',
				'Cluster Algorithm': '',
				'BonferroniTtest': 'No',
				'MaxMinScaler': 'No',
				'StandardScaler': 'No',
				'PCA': 'No',
				'Logarithm Transformation': 'No',
				'Statistical Features': 'No',
				'Silhouette': 0.0,
				'RandIndex': 0.0,
			}
		
		self.yesNo = {True : 'Yes',
					 False: 'No'}
		
		self.list_values = list()
		self.list_tuples= list ()
		
		return 
	
	def setName(self, name=''):
		self.name = name
		return self
	
	def setCluster(self, name=''):
		self.cluser = name
		return self
	
	def setBonf(self, used=False):
		self.usedBonf = used
		return self
	
	def setPca(self, used=False):
		self.usedPCA = used
		return self
	
	def setStandardScaler(self, used=False):
		self.usedSS = used
		return self
	
	def setMaxMinScaler(self, used=False):
		self.usedMMS = used
		return self
	
	def setClusteringAlghorithm(self, name=''):
		self.cluser = name
		return self 
	
	def setStatisticalFeatures(self, used=False):
		self.usedSF = used 
		return self
	
	def setLogarithmTransformation(self, used=False):
		self.usedLT = used
		return self
	
	def setSilhouette(self, value=0.0):
		self.silhouette = value
		return self
	
	def setRandIndex(self, value=0.0):
		self.randIndex = value
		return self
	
	def update(self):
		self.row_for_table = {
				'Omnic Name': self.name,
				'Cluster Algorithm': self.cluser,
				'BonferroniTtest': self.yesNo[self.usedBonf],
				'MaxMinScaler': self.yesNo[self.usedMMS],
				'StandardScaler': self.yesNo[self.usedSS],
				'PCA': self.yesNo[self.usedPCA],
				'Logarithm Transformation': self.yesNo[self.usedLT],
				'Statistical Features': self.yesNo[self.usedSF],
				'Silhouette': self.silhouette,
				'RandIndex': self.randIndex,
			}
		tuples = tuple()
		values = list()
		
		for key in self.row_for_table.keys():
			if key != 'Silhouette' and key != 'RandIndex':
				if len(tuples) == 0:
					tuples = (self.row_for_table[key],)
				else:
					tuples = tuples + (self.row_for_table[key],)
			else:
				values.append(self.row_for_table[key])
	
		self.list_tuples.append(tuples)
		self.list_values.append(values)
		
		return self
	
	def getDict(self, index):

		indexes = self.list_tuples[index]
		values = self.list_values[index]

		tmp = list(indexes) + values

		for k, value in zip(self.row_for_table.keys(), tmp):
			self.row_for_table[k] = value
		return self.row_for_table

	def getDictFromIndex(self, indexes, values):
		tmp = list(indexes) + values

		for value, k in zip(tmp, self.row_for_table.keys()):
			self.row_for_table[k] = value 

		return self.row_for_table

	def getValues(self):
		return self.row_for_table.values
	
	def getDF(self):
		name_col = [k for k in self.row_for_table.keys() if k!= 'Silhouette' and k != 'RandIndex']
		index = pd.MultiIndex.from_tuples(self.list_tuples, names=name_col)
		df = pd.DataFrame(self.list_values, columns=['Silhouette', 'RandIndex'], index=index)
		return df

	def getBoolForPipe(self, mydict):
		boolResult = dict()
		tmp = {
		'Yes' : True,
		'No' : False
		}

		for k in mydict.keys():
			if 'Yes' == mydict[k] or 'No' == mydict[k]:
				boolResult[k] = tmp[mydict[k]]

		return boolResult


	def maxSilhouette(self):
		df = self.getDF()
		index_max = df.idxmax()[0] # 0 to get the Silhouette
		max_score = df.max()[0]

		max_dict = self.getDictFromIndex(index_max, df.loc[index_max].to_list())
		return max_dict

	def maxRandIndex(self):
		df = self.getDF()
		index_max = df.idxmax()[1] # 0 to get the Silhouette
		max_score = df.max()[1]
		max_dict = self.getDictFromIndex(index_max, df.loc[index_max].to_list())
		return max_dict