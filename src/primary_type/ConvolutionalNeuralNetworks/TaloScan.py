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

p = {
     'batch_size': [64],
     'num_epochs': [200],
     'number_filter':[64],
     'filter_length':[4],
     'learning_rate':[0.02],
     'kernel_initializer': ['glorot_normal'],
     'max_pooling_size':[1],
     'activation_func':['relu'],
     'shapes': ['brick'],          # <<< required
     'first_neuron': [256],     # <<< required
     'hidden_layers':[4],    # <<< required
     'dropout': [0],          # <<< required
     'activation':['elu']
     }

# and run the experiment
scan_object = ta.Scan(x=Data.all_data,
                      y=Data.all_labels,
                      #x_val=Data.test_data, #using all data to train
                      #y_val=Data.test_labels,
                      model=Model.cnn_1d_model,
                      params=p,
                      experiment_name='cnn_cross_validation',
                      #performance_target=['val_acc', 0.96, False],
                      fraction_limit = 1,
                      #time_limit = '2019-09-25 10:00',
                      seed=Config.random_seed)
            
f = open("Talos_scan_object.pickle", "wb+")
pickle.dump(scan_object, f)
f.close()    

#from talos import Evaluate
#
# create the evaluate object
#e = Evaluate(scan_object)
#
## perform the evaluation
#
#tcga_x = np.expand_dims(Data.all_data, axis=2)
#tcga_y = Data.integer_encoded 
#####evaluate change the one hot encode to int encode, and not mentioned in talos doc.....
#print('Cross validation:')
#scan_object.evaluate_models(x=tcga_x,
#           y=tcga_y,
#           folds = 10,
#           metric='val_acc',
#           )
 

out = 'cnn_1d_tumor_type_t_'+Config.treat
if os.path.isdir(out):
    shutil.rmtree(out)
ta.Deploy(scan_object, out, metric='val_acc')

print("--- %s seconds ---" % (tt.time() - start_time))