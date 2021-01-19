import pandas as pd
import numpy as np
import pickle
import platform

from scipy import stats
from itertools import combinations

"""
	Class to apply the Bonferroni T Test.

"""

class Bonferroni_Ttest():

	def __init__(self, alpha=0.05, debug=False, label_case_id_into_X=False):
		self.alpha = alpha
		self.debug = debug
		self.condition = label_case_id_into_X
		self.pvalue_for_tumor = dict()
		self.pvalue_clean = dict()
		np.warnings.filterwarnings('ignore')

	def fit(self, X_, y=None):

		X = X_.copy()
		if self.condition is True:
			X.drop(['label', 'case_id'], axis=1, inplace=True)

		self.threshold = self.alpha / X.shape[1]
		if y is None:
			raise NameError("y must be an array!")

		# dict_to_dirk is a dictionary where the key is the label the value is
		# all the values of the dataset for that label

		if len(y.shape) > 1:
			label_in_dataset = y.loc[:, 'label'].unique()
		else:
			label_in_dataset = y.unique()

		self.dict_to_work = dict()
		for label in label_in_dataset:
			if len(y.shape) > 1:
				idx = y.index[y['label'] == label]
			else:
				idx = y.index[y == label]
			self.dict_to_work[label] = X.loc[idx, :]

		return self

	def transform(self, X, y=None):
		"""	
			The fuction works in the following way:
				1 - Create all the possibile combinations, with k=2, between the tumor types.
				2 - Than for each combination we calculate the ttest.
				3 - at the end each value which is greater than the threshold is selected to be dropped.
		"""
		X_ = X.copy()

		# if in X there are 'label' and 'case_id' features
		if self.condition is True:
			X_.drop(['case_id', 'label'], inplace=True, axis=1)

		if self.debug is True:
			for label in self.dict_to_work.keys():
				print(label,":", self.dict_to_work[label].shape)

		labels_tumor = [l for l in self.dict_to_work.keys() if str(l).__contains__('False') is False] # take just the tumor
		
		# all combinations for the set of labels tumor and normal
		comb_tumor = list(combinations(labels_tumor, r=2))

		# calculate p_value for tumor 
		p_value_tumor = list()
		for i, combo in enumerate(comb_tumor):
			if self.debug is True:
				print("step",i, combo[0], combo[1])

			arr1 = self.dict_to_work[combo[1]]
			arr2 = self.dict_to_work[combo[0]]
			_, p_value = stats.ttest_ind(np.array(arr1), np.array(arr2), equal_var=False)
			p_value_tumor.append(p_value)
			self.pvalue_for_tumor[combo] = p_value
			self.pvalue_clean[combo] = list()

		index_to_delete = set()

		for combo, p_values in zip(comb_tumor, p_value_tumor):
			for i, value in enumerate(p_values):
				if value > self.threshold:
					index_to_delete.add(i)
				else:
					self.pvalue_clean[combo].append(value)
			
		index_to_delete = list(index_to_delete)      
		X_.drop(X_.columns[index_to_delete], axis=1, inplace=True)

		print('Final dataset shape:', X_.shape)
		
		if self.condition is True:
			return pd.concat([X_, X.loc[:, ['label', 'case_id']]], axis=1)
		return X_.values

	def fit_transform(self, X, y):
		return self.fit(X,y).transform(X,y)
	