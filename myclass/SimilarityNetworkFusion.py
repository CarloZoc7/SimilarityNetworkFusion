import numpy as np
import pandas as pd

class SimilarityNetworkFusion:
    def __init__(self, dataset):
        self.dataset = dataset
    

    def weights(self):
        self.check_columns()

        for i, case_id in self.dataset.iterrows():
            
        return
    
    def similarity(self):
        return
    
    def check_columns(self):
        if 'label' in self.dataset.columns:
            self.dataset.drop(['label'], axis=1, inplace=True)
        if 'case_id' in self.dataset.columns:
            self.cases_id = self.dataset.loc[:, 'case_id']
            self.dataset.drop(['case_id'], axis=1, inplace=True)

        return