import numpy as np
import pandas as pd

from scipy.spatial.distance import pdist, squareform
from copy import deepcopy

class SimilarityNetworkFusion:
    def __init__(self, dataset, k=3):
        self.dataset = dataset.copy()
        self.cases_id = dataset.loc[:, 'case_id']
        self.k = k
        

    def weights(self):
        self.check_columns()
        df = pd.DataFrame(columns=self.cases_id, data=self.dataset.T.values)
        
        #calculate euclidean distance
        dist = pdist(self.dataset, 'euclidean')
        df_dist = pd.DataFrame(columns=self.cases_id, index=self.cases_id, data=squareform(dist))
        weights = df.cov()**2
                
        for i, patient_i in enumerate(tqdm(self.cases_id)):
            for patient_j in self.cases_id.iloc[i:]:
                    topK_i_mean = np.sort(df_dist.loc[patient_i, :].to_numpy())[:self.k].mean()
                    topK_j_mean = np.sort(df_dist.loc[patient_j, :].to_numpy())[:self.k].mean()
                    
                    mean = (topK_i_mean + topK_j_mean)/2

                    weights.loc[patient_i, patient_j] = np.exp(-(weights.loc[patient_i, patient_j]/mean))
                    weights.loc[patient_j, patient_i] = np.exp(-(weights.loc[patient_j, patient_i]/mean))
                    
        return weights
    
    def similarity(self):
        return
    
    def check_columns(self):
        if 'label' in self.dataset.columns:
            self.dataset.drop(['label'], axis=1, inplace=True)
        if 'case_id' in self.dataset.columns:
            self.dataset.drop(['case_id'], axis=1, inplace=True)
        return


    def find_k_neighbors(self, row, i, k=None): 
        row=deepcopy(row)
        #case of P matrix
        if k==None:
            del row[i]  #delete element of the same column of row index
            return row

        #case of S (find k elements with minimum distance value of W[i][j])
        else:
            k_neighbors_index=[]
            neigh = 0
            for j in range(0, len(row)):
                if j!=i:
                    min_index = row.index(min(row))
                    k_neighbors_index.append(min_index)
                    neigh+=1
                    del row[min_index]
                    if neigh == k:
                        return k_neighbors_index


    def P_matrix(self, W, n_case_id):
        P=[]
        for i in tqdm(range(0, n_case_id)):
            row=[]
            for j in range(0,n_case_id):
                if i==j:
                    row.append(1/2)

                else:
                    k_neighbors = self.find_k_neighbors(W[i], i)
                    denominator = 2*sum(k_neighbors)
                    row.append(W[i][j]/denominator)
            P.append(row)
        return P

    def S_matrix(self, W, n_case_id, k):
        S=[]
        for i in tqdm(range(0, n_case_id)):
            S_row=[]
            neighbors_indeces = self.find_k_neighbors(W[i], i, k)
            for j in range(0,n_case_id):
                if j not in neighbors_indeces:
                    S_row.append(0)

                else:
                    np_row = np.array(W[i])
                    denominator = sum(np_row[neighbors_indeces])
                    S_row.append(W[i][j]/denominator)
            S.append(S_row)
        return S