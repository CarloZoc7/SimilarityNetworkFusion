import numpy as np
import pandas as pd
import os

from scipy.spatial.distance import pdist, squareform, cdist
from copy import deepcopy
from sklearn.preprocessing import StandardScaler

class SimilarityNetworkFusion:
    def __init__(self, df_mirna, df_rna, df_illumina, k=3, mu=0.3):
        
        self.cases_id = df_rna.loc[:, 'case_id']
        self.rna = df_rna.copy()
        self.mirna = df_mirna.copy()
        self.illumina = df_illumina.copy()
        
        self.k = k
        self.mu = mu
        self.check_columns()
    
    def calculate_matrix(self):
        self.dict_dist = self.calculate_sim_matrix()
        if hasattr(self, 'w_rna') is False:
            self.w_rna = self.__weights__(self.rna, 'RNA', save_matrix=True)
            self.w_mirna = self.__weights__(self.mirna, 'miRNA', save_matrix=True)
            self.w_illumina = self.__weights__(self.illumina, 'Illumina', save_matrix=True)
        
        if hasattr(self, 'p_rna') is False:
            self.starting_p_rna = self.P_matrix(self.w_rna.to_numpy().tolist(), self.cases_id.shape[0], 'RNA', save_matrix=True)
            self.starting_p_mirna = self.P_matrix(self.w_mirna.to_numpy().tolist(), self.cases_id.shape[0], 'miRNA', save_matrix=True)
            self.starting_p_illumina = self.P_matrix(self.w_illumina.to_numpy().tolist(), self.cases_id.shape[0], 'Illumina', save_matrix=True)
                
        self.s_rna = self.S_matrix(self.w_rna.to_numpy().tolist(), self.cases_id.shape[0], 'RNA')
        self.s_mirna = self.S_matrix(self.w_mirna.to_numpy().tolist(), self.cases_id.shape[0], 'miRNA')
        self.s_illumina = self.S_matrix(self.w_illumina.to_numpy().tolist(), self.cases_id.shape[0], 'Illumina')
        
        return self
    
    def calculate_sim_matrix(self):
        tot = 0
        distance = 'euclidean'
        tot += pdist(self.rna, distance)
        tot += pdist(self.mirna, distance)
        tot += pdist(self.illumina, distance)
        tot = tot/3 
        df_dist = pd.DataFrame(columns=self.cases_id, index=self.cases_id, data=squareform(tot))
        
        return df_dist

    def __weights__(self, dataset, name, save_matrix=False):
        
        if 'weights_matrix_'+name+'.pkl' in os.listdir('.'):
            print('Read file pickle for weights matrix of {}'.format(name))
            weights = pd.read_pickle('weights_matrix_'+name+'.pkl')
            #dist = pdist(dataset, 'euclidean')
            #df_dist = pd.DataFrame(columns=self.cases_id, index=self.cases_id, data=squareform(dist))
            #self.dict_dist[name] = df_dist.copy()
            return weights
        
        print('Calculating weights for {}...'.format(name))
        df = pd.DataFrame(columns=self.cases_id, data=dataset.T.values)
        
        #calculate euclidean distance
        dist = pdist(dataset, 'euclidean')
        df_dist = pd.DataFrame(columns=self.cases_id, index=self.cases_id, data=squareform(dist))
        weights = pd.DataFrame(columns=self.cases_id, index=self.cases_id, data=[])
                
        for i, patient_i in enumerate(tqdm(self.cases_id)):
            for patient_j in self.cases_id.iloc[i:]:
                    tokK_mean_i = np.sort(df_dist.loc[patient_i, :].to_numpy())[:self.k].mean()
                    topK_mean_j = np.sort(df_dist.loc[patient_j, :].to_numpy())[:self.k].mean()
                    
                    eps = (tokK_mean_i + tokK_mean_i + df_dist.loc[patient_i, patient_j])/3

                    weights.loc[patient_i, patient_j] = np.exp(-(df_dist.loc[patient_i, patient_j]**2/(eps*self.mu)))
                    weights.loc[patient_j, patient_i] = np.exp(-(df_dist.loc[patient_j, patient_i]**2/(eps*self.mu)))
        if save_matrix:
            weights.to_pickle('weights_matrix_'+name+'.pkl')
        #self.dict_dist[name] = df_dist.copy()
        return weights       
    
    def check_columns(self):
        scaler = MinMaxScaler()
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


    def P_matrix(self, W, n_case_id, name, save_matrix=False):
        if 'pStarting_matrix_'+name+'.pkl' in os.listdir('.'):
            print('Reading the file pickle for the p starting matrix {}'.format(name))
            df_p = pd.read_pickle('pStarting_matrix_'+name+'.pkl')
            return df_p.to_numpy()
            
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
        #print(np.array(P))
        
        if save_matrix:
            df_P = pd.DataFrame(np.array(P))
            df_P.to_pickle('pStarting_matrix_'+name+'.pkl')
        return np.array(P)

    def S_matrix(self, W, n_case_id, name):
        print('Calculating S matrix for {}...'.format(name))
        S=[]
        for i in tqdm(range(0, n_case_id)):
            S_row=[]
            neighbors_indeces = self.find_k_neighbors(self.dict_dist.iloc[i,:].to_numpy().tolist(), i, self.k)
            for j in range(0,n_case_id):
                if j not in neighbors_indeces:
                    S_row.append(0)

                else:
                    np_row = np.array(W[i])
                    denominator = sum(np_row[neighbors_indeces])
                    S_row.append(W[i][j]/denominator)
            S.append(S_row)
        print(np.array(S))
        return np.array(S)
    
    def product_matrix(self, S_matrix, P_matrix):
        result = np.dot(S_matrix, P_matrix)
        result = np.dot(result, S_matrix.T)
        return result
    
    def sum_matrix_P(self, P1, P2):
        return np.add(P1,P2)/2
    
    def fit(self, num_iter=None):
        if num_iter is not None:
            self.p_rna = self.starting_p_rna.copy()
            self.p_mirna = self.starting_p_mirna.copy()
            self.p_illumina = self.starting_p_illumina.copy()
            for i in range(0, num_iter):
                self.p_rna_t1 = self.product_matrix(self.s_rna, self.sum_matrix_P(self.p_mirna, self.p_illumina))
                self.p_mirna_t1 = self.product_matrix(self.s_mirna, self.sum_matrix_P(self.p_rna, self.p_illumina))
                self.p_illumina_t1 = self.product_matrix(self.s_illumina, self.sum_matrix_P(self.p_mirna, self.p_rna))
                #print(self.p_rna_t1)
                #print(self.p_mirna_t1)
                #print(self.p_illumina_t1)
                self.p_rna = self.p_rna_t1
                self.p_mirna = self.p_mirna_t1
                self.p_illumina = self.p_illumina_t1
        else:
            print('ciao')

        return self
    
   
    def iterations_fit(self, matrices_diff=None, max_iter=100):
        if matrices_diff is not None:
            self.p_rna = self.starting_p_rna.copy()
            self.p_mirna = self.starting_p_mirna.copy()
            self.p_illumina = self.starting_p_illumina.copy()
            
            for step in range(0, max_iter):
                self.p_rna_t1 = self.product_matrix(self.s_rna, self.sum_matrix_P(self.p_mirna, self.p_illumina))
                self.p_mirna_t1 = self.product_matrix(self.s_mirna, self.sum_matrix_P(self.p_rna, self.p_illumina))
                self.p_illumina_t1 = self.product_matrix(self.s_illumina, self.sum_matrix_P(self.p_mirna, self.p_rna))
               
                self.p_rna = self.p_rna_t1
                self.p_mirna = self.p_mirna_t1
                self.p_illumina = self.p_illumina_t1

                diff_matrix = 0
                for i in range(0, len(self.p_rna)):
                    for j in range(i, len(self.p_rna)):
                        diff_matrix += np.abs(self.p_rna[i][j] - self.p_mirna[i][j])
                        diff_matrix += np.abs(self.p_illumina[i][j] - self.p_mirna[i][j])
                        diff_matrix += np.abs(self.p_illumina[i][j] - self.p_rna[i][j])
                
                diff_matrix = diff_matrix**0.5
                print(step, ':', diff_matrix)
                
                #diff_matrix = np.abs(np.subtract(self.p_rna, self.p_mirna)) + np.abs(np.subtract(self.p_rna, self.p_illumina)) + np.abs(np.subtract(self.p_mirna, self.p_illumina))
                #diff_matrix= np.abs(np.mean(diff_matrix))
                if diff_matrix<=np.abs(matrices_diff):
                    print('number of iterations to reach difference: ', step)
                    break
                    
                if step == max_iter-1: ##impossible to reach matrices difference
                    print('impossible to reach indicated difference, try with a bigger difference value')
        else:
            print('no difference for matrices found')

        return self
    
    def local_minimum_fit(self, iters_to_min=None, max_iter=100):
        if iters_to_min is not None:
            self.p_rna = self.starting_p_rna.copy()
            self.p_mirna = self.starting_p_mirna.copy()
            self.p_illumina = self.starting_p_illumina.copy()
            count=0
            prev_diff=0
            for step in range(0, max_iter):
                self.p_rna_t1 = self.product_matrix(self.s_rna, self.sum_matrix_P(self.p_mirna, self.p_illumina))
                self.p_mirna_t1 = self.product_matrix(self.s_mirna, self.sum_matrix_P(self.p_rna, self.p_illumina))
                self.p_illumina_t1 = self.product_matrix(self.s_illumina, self.sum_matrix_P(self.p_mirna, self.p_rna))
               
                self.p_rna = self.p_rna_t1
                self.p_mirna = self.p_mirna_t1
                self.p_illumina = self.p_illumina_t1

                diff_matrix = 0
                for i in range(0, len(self.p_rna)):
                    for j in range(i, len(self.p_rna)):
                        diff_matrix += np.abs(self.p_rna[i][j] - self.p_mirna[i][j])
                        diff_matrix += np.abs(self.p_illumina[i][j] - self.p_mirna[i][j])
                        diff_matrix += np.abs(self.p_illumina[i][j] - self.p_rna[i][j])
                
                diff_matrix = diff_matrix**0.5
                print(step, ':', diff_matrix)
                
                #diff_matrix = np.abs(np.subtract(self.p_rna, self.p_mirna)) + np.abs(np.subtract(self.p_rna, self.p_illumina)) + np.abs(np.subtract(self.p_mirna, self.p_illumina))
                #diff_matrix= np.abs(np.mean(diff_matrix))
                
                #check if a local minimum is found
                if int(diff_matrix)==prev_diff:
                    count+=1
                    if count>=iters_to_min:
                        print('local minimum reached in ', step, 'iterations')
                        break
                else:
                    count=0
                    
                prev_diff = int(diff_matrix)
                    
                if step == max_iter-1: ##impossible to reach matrices difference
                    print('impossible to reach local minimum, matrices seem to not converge')
        else:
            print('no minimum iterations for matrices found')

        return self
    
    def clean(self):
        del self.p_rna
        del self.p_mirna
        del self.p_illumina
        
        del self.p_rna_t1
        del self.p_mirna_t1
        del self.p_illumina_t1
        
        
        del self.w_rna
        del self.w_mirna
        del self.w_illumina
        
        return self