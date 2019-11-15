#ResNet data processing
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
import pickle

split_rate = 0.1
random_seed = 1024

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

train_data_resnet, test_data_resnet, train_labels_resnet, test_labels_resnet = train_test_split(all_data, labels, test_size = split_rate , random_state=random_seed)

#print("Before oversampling:")
#print(list(sum(train_labels)))

#####oversampling using SMOTE
#sm = SMOTE(random_state=Config.random_seed, sampling_strategy={4:400, 5:400, 23:400})
#train_data, train_labels = sm.fit_sample(train_data, train_labels)    

#print("Train data size:", train_data.shape)
#print("Test label size:", test_labels.shape)
#print("After oversampling")
#print(list(sum(train_labels)))
num_feature = df.shape[1]#actually the number of features for each x
num_class = labels.shape[1]


#Inception cnn1d data processing
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

import pickle

#from imblearn.over_sampling import SMOTE
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
#selected_cols = pd.read_csv(dataDir+'features_1073.csv', index_col = 0)
selected_cols = pd.read_csv(dataDir+'features_10_median.csv', index_col = 0)
#selected_cols = pd.read_csv(dataDir+'Cluster_Features_400.csv', index_col = 0)
selected_cols = list(selected_cols.iloc[:,0])#last col name is tumor type
print(selected_cols)
df = df[selected_cols]
#print(df.head())

common_columns = list(set(list(df.columns)) & set(list(df_pdx.columns)) &  set(list(df_clinical.columns)) &  set(list(df_meta.columns)))

print('Common Columns:')
common_columns.sort()
#common_columns = common_columns[:961]

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

train_data, test_data, train_labels, test_labels = train_test_split(df.values, labels, test_size = split_rate , random_state=random_seed)

print("Before oversampling:")
print(list(sum(train_labels)))

#####oversampling using SMOTE
#sm = SMOTE(random_state=Config.random_seed, sampling_strategy={4:400, 5:400, 23:400})
#train_data, train_labels = sm.fit_sample(train_data, train_labels)    

print("Train data size:", train_data.shape)
print("Test label size:", test_labels.shape)
print("After oversampling")
print(list(sum(train_labels)))

####generate history
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve
from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
import pprint
import os, sys
import time as tt
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.metrics import f1_score
import numpy as np
import pandas as pd
import talos as ta
#import Scan #required for pickle to load scan object correctly
import pickle
from keras.optimizers import SGD
import matplotlib.pyplot as plt

start_time = tt.time()

def reset_weights(model):
    from keras import backend as K
    session = K.get_session()
    for layer in model.layers: 
        if hasattr(layer, 'kernel_initializer'):
            layer.kernel.initializer.run(session=session)
            
            
def adjust_label(pred_labels, data_set):
    pred_labels = np.array(pred_labels, dtype=object)###no object option will truncate the string length
    if data_set == 'PDX':
        pred_labels[np.where((pred_labels=='LGG') | (pred_labels=='GBM'))] = 'LGG/GBM'
        pred_labels[np.where((pred_labels=='KIRC') | (pred_labels=='KIRP') | (pred_labels=='KIRH'))] = 'KIRC/KIRP/KIRH'
        pred_labels[np.where((pred_labels=='LUAD') | (pred_labels=='LUSC'))] = 'LUAD/LUSC'
        
    #print(pred_labels)
    return pred_labels


def train_keras(keras_model, best_talos_params, model_name, train_x, train_y, val_x, val_y):
    
    if model_name=='resnet':
        train_x = np.reshape(train_x, (train_x.shape[0], 32, -1))
        train_x = np.expand_dims(train_x, axis=3) 
        val_x = np.reshape(val_x, (val_x.shape[0], 32, -1))
        val_x = np.expand_dims(val_x, axis=3)
        
    elif model_name=='cnn1d' or model_name=='inception':
        
        train_x = np.expand_dims(train_x, axis=2)
        val_x = np.expand_dims(val_x, axis=2)
    else:
        print('Invalid model name')
        sys.exit()


    reset_weights(keras_model)


    sgd = SGD(lr=best_talos_params['learning_rate'], nesterov=best_talos_params['isNadam'])###Nadams
    keras_model.compile(loss='categorical_crossentropy',optimizer=sgd,metrics=['accuracy'])
    history = keras_model.fit(train_x, train_y, nb_epoch=best_talos_params['num_epochs'], validation_data=(val_x, val_y), batch_size=best_talos_params['batch_size'], verbose=0)
    
    return history
 
def savePickle(obj,output_name):
    f = open(output_name+".pickle", "wb+")
    pickle.dump(obj, f)
    f.close()
    
def loadPickle(output_name):
    f = open(output_name+".pickle", "rb")
    obj = pickle.load(f)
    f.close()
    return obj

talos_best_params = {'num_epochs': 500, 'batch_size':32, 'learning_rate':0.001, 'isNadam':False}
tag = '_6' ##791 e500_b32_r0.001_Adam
tag = '_7' ##241 e500_b32_r0.001_Adam

#scan_object = ta.Restore('resnet_tumor_type_t_ResNet_DEG_Feature_1000_smote_orderChromo_Best.zip')
#res_history = train_keras(scan_object.model, talos_best_params, 'resnet', train_data_resnet, train_labels_resnet, test_data_resnet, test_labels_resnet)
#savePickle(res_history,'resnet_history'+tag)

scan_object = ta.Restore('cnn_1d_tumor_type_t_CNN1d_DEG_Feature_10_smote_OrderChromo.zip')
cnn1d_history = train_keras(scan_object.model, talos_best_params, 'cnn1d', train_data, train_labels, test_data, test_labels)
savePickle(cnn1d_history,'cnn1d_history'+tag)

scan_object = ta.Restore('InceptionNet1d_tumor_type_t_InceptionNet1d_diverse_DEG_Feature_10_deep_orderChromo_scan.zip')
inception_history = train_keras(scan_object.model, talos_best_params, 'inception', train_data, train_labels, test_data, test_labels)
savePickle(inception_history,'inception_history'+tag)

#Terminal can't display
#SCALE_FACTOR = 1.25
#
#FIG_WIDTH_1COL = 5.2 * SCALE_FACTOR
#FIG_WIDTH_2COL = 7.5 * SCALE_FACTOR
#FIG_HEIGHT_MAX = 8.75 * SCALE_FACTOR
#
#for key in inception_history.history.keys():
#    fig = plt.figure(figsize=(FIG_WIDTH_1COL, FIG_HEIGHT_MAX))
#    plt.plot(inception_history.history[key])
#    plt.plot(res_history.history[key])
#    plt.plot(cnn1d_history.history[key])
#  
#    plt.title('model '+key)
#    plt.ylabel(key)
#    plt.xlabel('epoch')
#    plt.legend(['Inception', 'ResNet','CNN1d'], loc='upper left')
#
#    plt.show()
#    fig.savefig(key+'Visualization.png', dpi=600)
    
print("--- %s seconds ---" % (tt.time() - start_time))




