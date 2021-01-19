import pandas as pd
import numpy as np

from tqdm import tqdm

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score

"""
	Apply the SBS (Sequential backward Selection)

"""


class FeatureSelection:
	def __init__(self, model=KMeans(n_clusters=3), eval_op=silhouette_score):
		self.model = model
		self.eval_op = eval_op
		return

	def fit(self, X, y=None):
		return self

	def transform(self, X, y):
		subset_features = X.copy()
		sub_removed = list()
		tmp_removed = list()
		count = 0
		features_scores = pd.DataFrame(columns=['features', 'score'])

		for j in tqdm(range(0, X.shape[1])):
			local_scores = pd.DataFrame(columns=['features', 'score'])
            
			for col in X.columns:
				tmp_removed.append(col)      
				#print("		tmp removed (1):",tmp_removed)
				subset_features.drop(tmp_removed, inplace=True, axis=1)
				self.model.fit(subset_features)
				#print("Number features:", subset_features.shape[1], "score: ", self.eval_op(X, self.model.labels_))
				
				tmp_dict = {
					'features': tmp_removed.copy(),
					'score': self.eval_op(X, self.model.labels_)
					}
				features_scores = features_scores.append(tmp_dict, ignore_index=True)
				local_scores = local_scores.append(tmp_dict, ignore_index=True)
                
				tmp_removed.remove(col)
				#print("		tmp removed (2):",tmp_removed)

				subset_features = X.copy()
            
			# miglior locale
			tmp_removed = list()
			id_max = local_scores.loc[:, 'score'].idxmax()
			top_removed_features = features_scores.iloc[id_max, 0]
			tmp_removed = tmp_removed + top_removed_features
        
		best_score = features_scores.loc[:, 'score'].max()
		id_max = features_scores.loc[:, 'score'].idxmax()
		top_removed_features = features_scores.iloc[id_max, 0]

		#print(features_scores)
		print("Best_situation:", best_score, top_removed_features)

		X.drop(top_removed_features, axis=1, inplace=True)
		return X

	def fit_transform(self, X, y=None):
		return self.fit(X, y).transform(X, y)