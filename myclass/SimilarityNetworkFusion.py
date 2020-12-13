import numpy as np
import pandas as pd

from scipy.spatial.distance import pdist, squareform
from copy import deepcopy
from sklearn.preprocessing import StandardScaler

class SimilarityNetworkFusion:
    def __init__(self, df_mirna, df_rna, df_illumina, k=3):
        
        self.cases_id = df_rna.loc[:, 'case_id']
        self.rna = df_rna.copy()
        self.mirna = df_mirna.copy()
        self.illumina = df_illumina.copy()
        
        self.k = k
        self.check_columns()
    
    def calculate_matrix(self):
        self.w_rna = self.__weights__(self.rna, 'RNA')
        self.w_mirna = self.__weights__(self.mirna, 'miRNA')
        self.w_illumina = self.__weights__(self.illumina, 'Illumina')
        
        self.p_rna = self.P_matrix(self.w_rna.to_numpy().tolist(), self.cases_id.shape[0], 'RNA')
        self.p_mirna = self.P_matrix(self.w_mirna.to_numpy().tolist(), self.cases_id.shape[0], 'miRNA')
        self.p_illumina = self.P_matrix(self.w_illumina.to_numpy().tolist(), self.cases_id.shape[0], 'Illumina')
        
        self.s_rna = self.S_matrix(self.w_rna.to_numpy().tolist(), self.cases_id.shape[0], 'RNA')
        self.s_mirna = self.S_matrix(self.w_mirna.to_numpy().tolist(), self.cases_id.shape[0], 'miRNA')
        self.s_illumina = self.S_matrix(self.w_illumina.to_numpy().tolist(), self.cases_id.shape[0], 'Illumina')
        
        return self
        

    def __weights__(self, dataset, name):
        print('Calculating weights for {}...'.format(name))
        df = pd.DataFrame(columns=self.cases_id, data=dataset.T.values)
        
        #calculate euclidean distance
        dist = pdist(dataset, 'euclidean')
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
    
    def check_columns(self):
        scaler = StandardScaler()
        if 'label' in self.mirna.columns:
            self.mirna.drop(['label'], axis=1, inplace=True)
        if 'case_id' in self.mirna.columns:
            self.mirna.drop(['case_id'], axis=1, inplace=True)
            
        if 'label' in self.rna.columns:
            self.rna.drop(['label'], axis=1, inplace=True)
        if 'case_id' in self.rna.columns:
            self.rna.drop(['case_id'], axis=1, inplace=True)
            
        if 'label' in self.illumina.columns:
            self.illumina.drop(['label'], axis=1, inplace=True)
        if 'case_id' in self.illumina.columns:
            self.illumina.drop(['case_id'], axis=1, inplace=True)
            
        self.mirna = pd.DataFrame(scaler.fit_transform(self.mirna))
        self.rna = pd.DataFrame(scaler.fit_transform(self.rna))
        self.illumina = pd.DataFrame(scaler.fit_transform(self.illumina))

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
            max_value = max(row)
            for j in range(0, len(row)):
                if j!=i:
                    min_index = row.index(min(row))
                    k_neighbors_index.append(min_index)
                    neigh+=1
                    row[min_index] = max_value
                    if neigh == k:
                        return k_neighbors_index


    def P_matrix(self, W, n_case_id, name):
        print('Calculating P matrix for {}...'.format(name))
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
        return np.array(P)

    def S_matrix(self, W, n_case_id, name):
        print('Calculating S matrix for {}...'.format(name))
        S=[]
        for i in tqdm(range(0, n_case_id)):
            S_row=[]
            neighbors_indeces = self.find_k_neighbors(W[i], i, self.k)
            for j in range(0,n_case_id):
                if j not in neighbors_indeces:
                    S_row.append(0)

                else:
                    np_row = np.array(W[i])
                    denominator = sum(np_row[neighbors_indeces])
                    S_row.append(W[i][j]/denominator)
            S.append(S_row)
        return np.array(S)
    
    def product_matrix(self, S_matrix, P_matrix):
        result = np.dot(S_matrix, P_matrix)
        result = np.dot(result, S_matrix.T)
        return result
    
    def sum_matrix_P(self, P1, P2):
        return np.add(P1,P2)/2
    
    def fit(self, num_iter=None):
        if num_iter is not None:
            for i in range(0, num_iter):
                self.p_rna_t1 = self.product_matrix(self.s_rna, self.sum_matrix_P(self.p_mirna, self.p_illumina))
                self.p_mirna_t1 = self.product_matrix(self.s_mirna, self.sum_matrix_P(self.p_rna, self.p_illumina))
                self.p_illumina_t1 = self.product_matrix(self.s_illumina, self.sum_matrix_P(self.p_mirna, self.p_rna))

                self.p_rna = self.p_rna_t1
                self.p_mirna = self.p_mirna_t1
                self.p_illumina = self.p_illumina_t1
        else:
            print('ciao')
        return