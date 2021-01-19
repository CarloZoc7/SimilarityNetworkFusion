import pandas as pd
import numpy as np
import pickle
import platform

"""
	Class that allows to clean and merge the tumor and normal dataset.

"""

class Clean_Merge_Dataset:
	def __init__(self, name=''):
		self.path_to_final = './data-ready/final_dataset_'+name+'.pkl'
		if platform.system() == 'Windows':
			self.path_to_final = self.path_to_final.replace('/','\\')
		return

	def fit(self, X, y=None):
		"""
			function created to add the class element to the Pipeline function.
		"""
		return self
		
	def transform(self, data_normal, data_tumor, X=None, y=None):
		"""
			The process trasform is the following:
				1 - removing the label 'TCGA-MESO' and 'CPTAC-3'
				2 - deleting 0-values and Nan-values
				3 - Return the dataset, the labels and the cases_id to the main function
		"""
		print('Data_normal:', data_normal.shape)
		print('Data_tumor:', data_tumor.shape)

		dataset = pd.concat([data_tumor, data_normal], ignore_index=True)
		print('All data:', dataset.shape)
		
		# removing TCGA-MESO e CPTAC-3 label items
		dataset = dataset[dataset['label'] != 'TCGA-MESO']
		dataset = dataset[dataset['label'] != 'CPTAC-3']
		
		print(set(dataset['label']))
		
		# union target with label
		for index, element in dataset.iterrows():
			if element['target'] is False:
				element['label'] = element['target']

			dataset.at[index, 'label'] = element['label']
		dataset.drop(['target'], inplace=True, axis=1)
		
		# Check null and zero value
		# count non zero value

		# delete before the 0-values features
		sum_count = 0
		index_to_delete = list()
		for i, element in enumerate(dataset.isin([0]).sum()):
			if element == dataset.shape[0]:
				sum_count += 1
				index_to_delete.append(i)
			
		dataset.drop(dataset.columns[index_to_delete], inplace=True, axis=1) # delete 0-values features
		print('Features completly 0 values', sum_count, 'removed')

		# counting Nan values
		sum_count = 0
		index_to_delete = list()
		for i, element in enumerate(dataset.isna().sum()):
			if element == dataset.shape[0]:
				sum_count+=1
				index_to_delete.append(i)

		dataset.dropna(inplace=True, axis=1)
		print('Features completely Nan', sum_count, 'removed')
		
		print('Final dataset shape', dataset.shape)
		dataset.to_pickle(self.path_to_final)
		
		del data_normal
		del data_tumor

		y = dataset.loc[:, 'label']
		cases_id = dataset.loc[:, 'case_id']
		X = dataset.drop(columns=['label', 'case_id'])
		#X = dataset.iloc[:, dataset.columns != 'case_id']
		
		return X, y, cases_id
