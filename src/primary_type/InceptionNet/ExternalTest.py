from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve
from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
import tensorflow as tf
import pprint
import os, sys
import time as tt
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.metrics import f1_score
import numpy as np
import pandas as pd
import talos as ta
import pickle
from keras.optimizers import SGD

sys.path.insert(1, '/projects/compsci/Yue/SubtypeClassifier/Helpers/')
import Config
import Model
from Model import ExternalValidate, cross_validate_keras, reshape_data_1d, reset_weights
import DataRna as Data

start_time = tt.time()

scan_object = ta.Restore('InceptionNet1d_tumor_type_t_'+Config.treat+'.zip')
print(scan_object.params)
talos_best_params = {'num_epochs': 500, 'batch_size':32, 'learning_rate':0.001, 'isNadam':False}

scan_object.model.summary()

train_x_all = Data.all_data 
train_y_all = Data.all_labels
train_x_all = np.expand_dims(train_x_all, axis=2)
cross_validate_keras(scan_object.model, talos_best_params, 10, train_x_all, train_y_all)

####retrain with all data
keras_model = scan_object.model
reset_weights(keras_model)

sgd = SGD(lr=talos_best_params['learning_rate'], nesterov=talos_best_params['isNadam'])###Nadams
keras_model.compile(loss='categorical_crossentropy',optimizer=sgd,metrics=['accuracy'])
keras_model.fit(train_x_all, train_y_all, nb_epoch=talos_best_params['num_epochs'], batch_size=talos_best_params['batch_size'], verbose=0)


test_x = reshape_data_1d(Data.pdx_X)
ExternalValidate(test_x, Data.pdx_Y, 'PDX', keras_model)

test_x = reshape_data_1d(Data.clinical_X)
ExternalValidate(test_x, Data.clinical_Y, 'Clinical', keras_model)

test_x = reshape_data_1d(Data.meta_X)
ExternalValidate(test_x, Data.meta_Y, 'Meta', keras_model)



print("--- %s seconds ---" % (tt.time() - start_time))