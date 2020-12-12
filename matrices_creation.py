import numpy as np
from copy import deepcopy

def find_k_neighbors(row, i, k=None): 
    row=deepcopy(row)
    #case of P matrix
    if k==None:
        del row[i]  #delete element of the same column of row index
        return row
    
    #case of S (find k elements with maximimum similarity value of W[i][j])
    else:
        k_neighbors_index=[]
        neigh = 0
        for j in range(0, len(row)):
            if j!=i:
                max_index = row.index(max(row))
                k_neighbors_index.append(max_index)
                neigh+=1
                del row[max_index]
                if neigh == k:
                    return k_neighbors_index


def P_matrix(W, n_case_id):
    P=[]
    for i in range(0, n_case_id):
        row=[]
        for j in range(0,n_case_id):
            if i==j:
                row.append(1/2)
            
            else:
                k_neighbors = find_k_neighbors(W[i], i)
                denominator = 2*sum(k_neighbors)
                row.append(W[i][j]/denominator)
        P.append(row)
    return P

def S_matrix(W, n_case_id, k):
    S=[]
    for i in range(0, n_case_id):
        S_row=[]
        neighbors_indeces = find_k_neighbors(W[i], i, k)
        for j in range(0,n_case_id):
            if j not in neighbors_indeces:
                S_row.append(0)
            
            else:
                np_row = np.array(W[i])
                denominator = sum(np_row[neighbors_indeces])
                S_row.append(W[i][j]/denominator)
        S.append(S_row)
    return S
            