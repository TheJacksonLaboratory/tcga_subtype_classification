import os
import numpy as np
import pandas as pd

dropdown_headers = ['Tumor type', 'Cancer type', 'Cancer subtype']
dropdown_columns = ['tumor.type.matched', 'cancertype', 'subtype']

def read_data():
    data_file = os.path.join(os.path.dirname(__file__),
                             'data/sample_dr.csv')
    sample_data = pd.read_csv(data_file, header=0, index_col=0)
    sample_data.fillna('None', inplace=True)

    sample_types = sample_data['sample.type'].astype(int).values
    tumor_mask = sample_types < 11
    sample_data.loc[tumor_mask, 'tumor.type.matched'] =\
        sample_data.loc[tumor_mask, 'tumor.type'] + ' T'
    sample_data.loc[~tumor_mask, 'tumor.type.matched'] =\
        sample_data.loc[~tumor_mask, 'tumor.type'] + ' N'

    return sample_data

#def get_tumortype_checkboxes(dataframe):
#    # matched normals are anything above 11
#    subdata = dataframe[['tumor.type', 'sample.type']]
#    checkboxes = {

def get_available_checkboxes(dataframe):

    available_checkboxes = dict((col, sorted(dataframe[col].unique()))
                                for col in dropdown_columns)

    for col, checkboxes in available_checkboxes.items():
        if 'None' not in checkboxes: continue
        new_checkboxes = ['None'] + sorted(set(checkboxes) - set(['None']))
        available_checkboxes[col] = new_checkboxes

    for col, checkboxes in available_checkboxes.items():
        lengths = [sum(dataframe[col] == cb) for cb in checkboxes]
        new_checkboxes = [{'label': f'{cb}  ({length})', 'value': cb}
                          for cb, length in zip(checkboxes, lengths)]
        available_checkboxes[col] = new_checkboxes

    return available_checkboxes
