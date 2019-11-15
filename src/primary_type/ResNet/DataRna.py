import numpy as np
import pandas as pd
from math import sqrt
import random
import copy
import matplotlib.pyplot as plt
import sys
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
import os
import Config
import pickle

from imblearn.over_sampling import SMOTE
#from Model import quantileNormalize


#tf.enable_eager_execution()
    
dataDir = '../DataPool/'
df = pd.read_csv(dataDir+'tcga_training.csv', index_col = 0)
string_labels = df['tumor.type']

# integer encode
values = np.array(string_labels)

label_encoder = LabelEncoder()
integer_encoded = label_encoder.fit_transform(values)


#generate mapping from string label to int
from_int_to_label = {}
for i in range(len(integer_encoded)):
    if integer_encoded[i] in from_int_to_label:
        continue
    else:
        from_int_to_label[integer_encoded[i]] = values[i]

from_label_to_int = dict((from_int_to_label[k], k) for k in from_int_to_label)
print('The string labels are %s', from_int_to_label)

# binary encode
onehot_encoder = OneHotEncoder(sparse=False)
integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
onehot_encoded = onehot_encoder.fit_transform(integer_encoded)
labels = onehot_encoded

df = df.drop(['tumor.type'], axis = 1)

####################preprocessing
#df = df.loc[:, (df != 0).any(axis=0)] #drop all zero columns
#df /= df.max() #scale by column so that range from 0 to 1


num_class = labels.shape[1]

    
######processing pdx data###############
df_pdx = pd.read_csv(dataDir+'ExternalDataPdx.csv', index_col = 0)
df_pdx.dropna(axis='index', inplace=True) #NA values in tumor.type
string_labels_pdx = df_pdx['tumor.type']

# integer encode
# integer encode
#values = np.array(string_labels_pdx)
#print(list(np.unique(values)))
############map current label to the 33 labels

#int_encoded = ['NA' for _ in string_labels]
#for i,l in enumerate(string_labels):
#    if l == 'NotSure':
#        continue
#    else:
#        int_encoded[i] = from_label_to_int[l]
#
#print(int_encoded)
#labels_ex_pdx = []
#for j in int_encoded:
#    ar = [0 for _ in range(num_class)]
#    ar[j] = 1
#    labels_ex_pdx.append(ar)
#
#labels_ex_pdx = np.array(labels_ex_pdx)


df_pdx = df_pdx.drop(['tumor.type'], axis=1)
#df_ex = df_ex/df_ex.max()

######processing Clinical data###############
df_clinical = pd.read_csv(dataDir+'/ExternalDataClinical.csv', index_col = 0)

print(df_clinical.head())
string_labels_clinical = df_clinical['tumor.type']
df_clinical = df_clinical.drop(['tumor.type'], axis=1)

######processing TCGA Meta data###############
df_meta = pd.read_csv(dataDir+'/ExternalDataMeta.csv', index_col = 0)

print(df_meta.head())
string_labels_meta = df_meta['tumor.type']
df_meta = df_meta.drop(['tumor.type'], axis=1)

    
##########extract previous work features
print("Using Selected features.")
selected_cols = pd.read_csv(dataDir+'features_70_median.csv', index_col = 0)
#selected_cols = pd.read_csv(dataDir+'features_1073.csv', index_col = 0)
#selected_cols = pd.read_csv(dataDir+'features_254.csv', index_col = 0)
#selected_cols = pd.read_csv(dataDir+'Cluster_Features_400.csv', index_col = 0)
selected_cols = list(selected_cols.iloc[:,0])#last col name is tumor type
print(selected_cols)
df = df[selected_cols]
#print(df.head())

common_columns = list(set(list(df.columns)) & set(list(df_pdx.columns)) &  set(list(df_clinical.columns)) &  set(list(df_meta.columns)))

print('Common Columns:')
common_columns.sort()

common_columns = common_columns[:1024]

####order the genes by location
geneOrders = pd.read_csv(dataDir+'GeneLocations.csv', index_col = 0)
geneOrders = geneOrders['x.entid']
ll = []
for gene in geneOrders:
    if ('X'+str(gene)) in common_columns:
        ll.append(('X'+str(gene)))
common_columns = ll
###############################
print(common_columns)

df_pdx = df_pdx[common_columns]
df_clinical = df_clinical[common_columns]
df_meta = df_meta[common_columns]
df = df[common_columns]

df = df.loc[:,~df.columns.duplicated()]#pandas wont remove duplicate column when doing submatrix

train_data, test_data, train_labels, test_labels = train_test_split(df.values, labels, test_size = Config.split_rate , random_state=Config.random_seed)

print("Before oversampling:")
print(list(sum(train_labels)))

#####oversampling using SMOTE
#sm = SMOTE(random_state=Config.random_seed, sampling_strategy={4:400, 5:400, 23:400})
#train_data, train_labels = sm.fit_sample(train_data, train_labels)    

print("Train data size:", train_data.shape)
print("Test label size:", test_labels.shape)
print("After oversampling")
print(list(sum(train_labels)))
num_feature = df.shape[1]#actually the number of features for each x
num_class = labels.shape[1]
#df /= df.max()
#df_ex = quantileNormalize(df_input=df_ex, reference_df=df)


#####################training data processing
#df = df.applymap(lambda x: np.log2(x+1))

# demonstrate data standardization with sklearn
from sklearn.preprocessing import StandardScaler
all_data = df.transpose().values

scaler = StandardScaler()
scaler.fit(all_data)
all_data = scaler.transform(all_data)
all_data = np.transpose(all_data)

all_labels = labels
print('Training data:')
print(all_data)

print("TCGA data size")
print(all_data.shape)
print('Whole label size')
print(all_labels.shape)

######################external

pdx_X = df_pdx.transpose().values
scaler.fit(pdx_X)
pdx_X = scaler.transform(pdx_X)
pdx_X = np.transpose(pdx_X)

pdx_Y = string_labels_pdx

print("PDX labels:", string_labels_pdx)
print("PDX data size:", pdx_X.shape)
print("PDX label size:", pdx_Y.shape)
print(pdx_X)


clinical_X = df_clinical.transpose().values
scaler.fit(clinical_X)
clinical_X = scaler.transform(clinical_X)
clinical_X = np.transpose(clinical_X)

clinical_Y = string_labels_clinical
clinical_Y_integer = label_encoder.fit_transform(clinical_Y)###same encoder for tcga data
clinical_Y_integer = clinical_Y_integer.reshape(len(clinical_Y_integer), 1)
clinical_Y_onehot = onehot_encoder.transform(clinical_Y_integer)

print("clinical labels:", string_labels_clinical)
print("clinical data size:", clinical_X.shape)
print("clinical label size:", clinical_Y.shape)
print(clinical_X)
#print(clinical_X)

meta_X = df_meta.transpose().values
scaler.fit(meta_X)
meta_X = scaler.transform(meta_X)
meta_X = np.transpose(meta_X)

meta_Y = string_labels_meta
meta_Y[meta_Y=='COAD'] = 'COADREAD'

print("meta labels:", string_labels_meta)
print("meta data size:", meta_X.shape)
print("meta label size:", meta_Y.shape)
print(meta_X)



#df_test = pd.DataFrame({'C1': {'A': 5, 'B': 2, 'C': 3, 'D': 4},
#                   'C2': {'A': 4, 'B': 1, 'C': 4, 'D': 2},
#                   'C3': {'A': 3, 'B': 4, 'C': 6, 'D': 8}})
#                   
#df_ref = pd.DataFrame({'C1': {'A': 0, 'B': 2, 'C': 0, 'D': 4},
#                   'C2': {'A': 0, 'B': 1, 'C': 0, 'D': 2},
#                   'C3': {'A': 0, 'B': 0, 'C': 6, 'D': 0}})
#                   
#result = quantileNormalize(df_test, df_ref)
#print(result)



