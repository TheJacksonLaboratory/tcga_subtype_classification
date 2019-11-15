###Talos scan

from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve
from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
import tensorflow as tf
import pprint
import os, sys
import time as tt
from sklearn.metrics import confusion_matrix
import numpy as np
import pandas as pd

import talos as ta
import wrangle as wr
from keras.optimizers import Nadam, Adam
import pickle
import shutil

sys.path.insert(1, '/projects/compsci/Yue/SubtypeClassifier/Helpers/')
import Config
import Model
import DataRna as Data

start_time = tt.time()
####1d inception best feature_40
p = {
     'batch_size': [32],
     'num_epochs': [200],
     'number_filter':[8],
     'filter_length_11':[1],
     'filter_length_21':[100],
     'filter_length_22':[5],
     'filter_length_31':[1],
     'filter_length_32':[1],
     'filter_length_33':[10],
     'learning_rate':[0.1],
     #'kernel_initializer': ['normal', 'glorot_normal'],
     'max_pooling_size_1':[1],
     'max_pooling_size_2':[10],
     'max_pooling_size_3':[1],
     #'activation_func':['relu', 'elu'],
     #'shapes': ['brick'],          # <<< required
     #'first_neuron': [256,128,64,32],     # <<< required
     #'hidden_layers':[0,1,2,3,4,5],    # <<< required
     'dropout_1': [0.3],          # <<< required
     'dropout_2': [0.5],          # <<< required
     'dropout_3': [0.9],          # <<< required
     'dropout_4': [0.1],          # <<< required
     'isNadam':[False]
     #'activation':['relu', 'elu']
     }
    
#p={
#    'batch_size': [32],
#    'num_epochs': [200],
#    'number_filter':[128,64,32,8],
#    'filter_length_11':[1,5,10,50,100], #only the first node in each module can have 100-filter
#    'filter_length_21':[1,5,10,50,100],
#    'filter_length_22':[1,5,10,50],
#    'filter_length_31':[1,5,10,50,100],
#    'filter_length_32':[1,5,10,50],
#    'filter_length_33':[1,5,10,50],
#    'learning_rate':[0.1],
#    'max_pooling_size_1':[1,5,10], #cannot be too big, will cause negative dim
#    'max_pooling_size_2':[1,5,10],
#    'max_pooling_size_3':[1,5,10],
#    'dropout_1': [0, 0.1,0.3,0.5,0.7,0.9],          
#    'dropout_2': [0, 0.1,0.3,0.5,0.7,0.9],          
#    'dropout_3': [0, 0.1,0.3,0.5,0.7,0.9],          
#    'dropout_4': [0, 0.1,0.3,0.5,0.7,0.9],          
#    'isNadam':[False]
#
#}     
     
# and run the experiment
scan_object = ta.Scan(x=Data.all_data,
                      y=Data.all_labels,
                      #x_val=Data.test_data,#using all data to train
                      #y_val=Data.test_labels,
                      model=Model.inception_1d_model_deeper,
                      params=p,
                      experiment_name='inception_1d_cross_validation',
                      #performance_target=['val_acc', 0.98, False],
                      fraction_limit = 0.5,
                      #time_limit = '2019-09-25 10:00',
                      seed=Config.random_seed)
            
f = open("Talos_scan_object_inception_1d.pickle", "wb+")
pickle.dump(scan_object, f)
f.close()    
 

out = 'InceptionNet1d_tumor_type_t_'+Config.treat
if os.path.isdir(out):
    print('removing previous dir')
    shutil.rmtree(out)
    os.remove(out+'.zip')#has to manually remove existing zip
ta.Deploy(scan_object, out, metric='val_acc')

print("--- %s seconds ---" % (tt.time() - start_time))