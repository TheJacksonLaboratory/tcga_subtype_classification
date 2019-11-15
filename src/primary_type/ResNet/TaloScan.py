###Pure Cross Validation
import Config
import DataRna as Data
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

sys.path.insert(1, '../Helpers/')
import Model


start_time = tt.time()

#p = {
#     'batch_size': [256,128,64,32],
#     'num_epochs': [200,300,500,1000],
#     'learning_rate':(0.01,0.1,10),
#     'kernel_initializer': [None],
#     'dropout': (0,1,10),          # <<< required
#     'isNadam':[True, False]
#     }
     
     
p = {
     'batch_size': [64],
     'num_epochs': [500],
     'learning_rate':[0.064],
     'kernel_initializer': [None],
     'dropout': [0],          # <<< required
     'isNadam':[True]
     }     

# and run the experiment
scan_object = ta.Scan(x=Data.all_data,
                      y=Data.all_labels,
                      x_val=Data.test_data,
                      y_val=Data.test_labels,
                      model=Model.resnet_model,
                      params=p,
                      experiment_name='resnet_cross_validation',
                      performance_target=['val_acc', 0.98, False],
                      fraction_limit = 1,
                      seed=Config.random_seed)
            
f = open("Talos_scan_object_resnet.pickle", "wb+")
pickle.dump(scan_object, f)
f.close()
 

out = 'resnet_tumor_type_t_'+Config.treat
if os.path.isdir(out):
    shutil.rmtree(out)
ta.Deploy(scan_object, out, metric='val_acc')

print("--- %s seconds ---" % (tt.time() - start_time))