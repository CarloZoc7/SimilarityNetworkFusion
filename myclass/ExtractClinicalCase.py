import os
import json
import pandas as pd

"""
    Class to extract the clinical_case from the json file.
    It returns a dataframe in the following form:
    case_id| ...[clinical case]...| label.
"""


class ExtractClinicalCase:
    def __init__(self, cases_id):
        with open('./data-ready/clinical.cases_selection.2020-11-12.json', 'r') as json_file:
            data = json.load(json_file)

        remove_el = list()
        for el in data:
            if el['case_id'] not in cases_id:
                remove_el.append(el)

        for el in remove_el:
            data.remove(el)
        
        clinical_data = {'case_id': '',
                        'tumor_stage': '',
                        'prior_malignancy': '',
                        'age_at_diagnosis': None,
                        'morphology': '',
                        'label': ''}

        self.df = pd.DataFrame(data=[], columns=clinical_data.keys())

        for i, el in enumerate(data):
            clinical_data['case_id'] = el['case_id']
            clinical_data['tumor_stage'] = el['diagnoses'][0]['tumor_stage']
            clinical_data['prior_malignancy'] = el['diagnoses'][0]['prior_malignancy']
            if el['diagnoses'][0]['age_at_diagnosis'] is not None:
                value = int(el['diagnoses'][0]['age_at_diagnosis'])/365
                clinical_data['age_at_diagnosis'] = self.__truncate__(value)

            clinical_data['morphology'] = el['diagnoses'][0]['morphology']

            self.df = self.df.append(pd.DataFrame(clinical_data, index=[i]), ignore_index=True)
        
    def get_df_clinical_case(self):
        return self.df

    def __truncate__(self, n, decimals=-1):
        """
            Function to take the decade of the age.
        """
        multiplier = 10 ** decimals
        return int(n * multiplier) / multiplier
                

