import Config
import tensorflow as tf
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve
from sklearn.metrics import f1_score
import pandas as pd
import sys
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels
import DataRna as Data
from sklearn.model_selection import KFold

import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from keras.models import Sequential
from keras.layers import Dense, Activation, Flatten, Convolution1D, Dropout, MaxPooling1D
from keras.optimizers import SGD
from keras.utils import np_utils
from talos.utils import hidden_layers
from sklearn.metrics import confusion_matrix, classification_report
from keras.layers import Conv1D, MaxPooling1D
from keras.layers import Flatten, Dense
from keras.layers import Input
from keras.layers import concatenate


from keras.models import load_model
from sklearn.datasets import load_files   
from keras.utils import np_utils
from glob import glob
from keras import applications
from keras.preprocessing.image import ImageDataGenerator 
from keras import optimizers
from keras.models import Sequential,Model,load_model
from keras.layers import Dense, Dropout, Flatten, Conv2D, MaxPool2D,GlobalAveragePooling2D
from keras.callbacks import TensorBoard,ReduceLROnPlateau,ModelCheckpoint
import talos as ta
import pickle


from keras.models import load_model
from sklearn.datasets import load_files   
from keras.utils import np_utils
from keras import applications
from keras.preprocessing.image import ImageDataGenerator 
from keras import optimizers
from keras.models import Sequential,Model,load_model
from keras.layers import Dense, Dropout, Flatten, Conv2D, MaxPool2D,GlobalAveragePooling2D
from keras.callbacks import TensorBoard,ReduceLROnPlateau,ModelCheckpoint
from keras.optimizers import SGD



def reshape_data_2d(data):
    data = np.reshape(data, (data.shape[0], 32, -1))
    data = np.expand_dims(data, axis=3)
    return data

def reshape_data_1d(data):
    data = np.expand_dims(data, axis=2)
    return data

def top_n_accuracy(truths, preds, n):
    #truths: int true label
    #preds: vector prob prediction
    
    best_n = np.argsort(preds, axis=1)[:,-n:]#index sorted
    #ts = np.argmax(truths, axis=1)
    ts = np.array(truths)
    successes = 0

    for i in range(ts.shape[0]):
      if ts[i] in list(best_n[i,:]):
        successes += 1

    return float(successes)/ts.shape[0]


def adjust_label(pred_labels, data_set):
    pred_labels = np.array(pred_labels, dtype=object)###no object option will truncate the string length
    if data_set == 'PDX':
        pred_labels[np.where((pred_labels=='LGG') | (pred_labels=='GBM'))] = 'LGG/GBM'
        pred_labels[np.where((pred_labels=='KIRC') | (pred_labels=='KIRP') | (pred_labels=='KIRH'))] = 'KIRC/KIRP/KIRH'
        pred_labels[np.where((pred_labels=='LUAD') | (pred_labels=='LUSC'))] = 'LUAD/LUSC'
        

    return pred_labels
def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 10)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value


def resnet_model(x_train, y_train, x_val, y_val, params):


    x_train = np.reshape(x_train, (x_train.shape[0], 32, -1))
    x_train = np.expand_dims(x_train, axis=3)
    
    x_val = np.reshape(x_val, (x_val.shape[0], 32, -1))
    x_val = np.expand_dims(x_val, axis=3)
    
    base_model = applications.resnet50.ResNet50(weights=None, 
                                                include_top=False, 
                                                input_shape=(32,32,1))
    x = base_model.output
    x = GlobalAveragePooling2D()(x)
    x = Dropout(params['dropout'])(x)
    predictions = Dense(Data.num_class, activation= 'softmax')(x)
    model = Model(inputs = base_model.input, outputs = predictions)

    sgd = SGD(lr=params['learning_rate'], nesterov=params['isNadam'])###Nadams
    model.compile(loss='categorical_crossentropy',optimizer=sgd,metrics=['accuracy'])

    nb_epoch = params['num_epochs']
    history = model.fit(x_train, y_train, nb_epoch=nb_epoch, validation_data=(x_val, y_val), batch_size=params['batch_size'], verbose=0)
    
    return history, model

    
    
# first we have to make sure to input data and params into the function
def inception_1d_model(x_train, y_train, x_val, y_val, params):

    x_train = np.expand_dims(x_train, axis=2)
    x_val = np.expand_dims(x_val, axis=2)
    
    input = Input(shape = (x_train.shape[1],1))
    tower_1 = Conv1D(params['number_filter'], params['filter_length_11'], activation='relu')(input)
    tower_1 = MaxPooling1D(pool_size=params['max_pooling_size_1'])(tower_1)
    tower_1 = Dropout(params['dropout_1'])(tower_1)

    
    tower_2 = Conv1D(params['number_filter'], params['filter_length_21'], activation='relu')(input)
    tower_2 = Conv1D(params['number_filter'], params['filter_length_22'], activation='relu')(tower_2)
    tower_2 = MaxPooling1D(pool_size=params['max_pooling_size_2'])(tower_2)
    tower_2 = Dropout(params['dropout_2'])(tower_2)

    
    tower_3 = Conv1D(params['number_filter'], params['filter_length_31'] ,activation='relu')(input)
    tower_3 = Conv1D(params['number_filter'], params['filter_length_32'] ,activation='relu')(tower_3)
    tower_3 = Conv1D(params['number_filter'], params['filter_length_33'] ,activation='relu')(tower_3)
    tower_3 = MaxPooling1D(pool_size=params['max_pooling_size_3'])(tower_3)
    tower_3 = Dropout(params['dropout_3'])(tower_3)
    
    output = concatenate([tower_1, tower_2, tower_3], axis = 1)
    output = Flatten()(output)
    output = Dropout(params['dropout_4'])(output)
    predictions = Dense(Data.num_class, activation='softmax')(output)
    
    model = Model(inputs = input, outputs = predictions)

    sgd = SGD(lr=params['learning_rate'], nesterov=params['isNadam'])###Nadams
    model.compile(loss='categorical_crossentropy',optimizer=sgd,metrics=['accuracy'])

    nb_epoch = params['num_epochs']
    history = model.fit(x_train, y_train, nb_epoch=nb_epoch, validation_data=(x_val, y_val), batch_size=params['batch_size'], verbose=0)
    
    return history, model


def inception_1d_model_deeper(x_train, y_train, x_val, y_val, params):

    x_train = np.expand_dims(x_train, axis=2)
    x_val = np.expand_dims(x_val, axis=2)
    
    input = Input(shape = (x_train.shape[1],1))
    tower_1 = Conv1D(params['number_filter'], params['filter_length_11'], activation='relu')(input)
    tower_1 = MaxPooling1D(pool_size=params['max_pooling_size_1'])(tower_1)
    tower_1 = Dropout(params['dropout_1'])(tower_1)

    
    tower_2 = Conv1D(params['number_filter'], params['filter_length_21'], activation='relu')(input)
    tower_2 = Conv1D(params['number_filter'], params['filter_length_22'], activation='relu')(tower_2)
    tower_2 = MaxPooling1D(pool_size=params['max_pooling_size_2'])(tower_2)
    tower_2 = Dropout(params['dropout_2'])(tower_2)

    
    tower_3 = Conv1D(params['number_filter'], params['filter_length_31'] ,activation='relu')(input)
    tower_3 = Conv1D(params['number_filter'], params['filter_length_32'] ,activation='relu')(tower_3)
    tower_3 = Conv1D(params['number_filter'], params['filter_length_33'] ,activation='relu')(tower_3)
    tower_3 = MaxPooling1D(pool_size=params['max_pooling_size_3'])(tower_3)
    tower_3 = Dropout(params['dropout_3'])(tower_3)
    
    output = concatenate([tower_1, tower_2, tower_3], axis = 1)
    output = Flatten()(output)
    output = Dropout(params['dropout_4'])(output)
    output = Dense(100, activation='relu')(output)
    predictions = Dense(Data.num_class, activation='softmax')(output)
    
    model = Model(inputs = input, outputs = predictions)

    sgd = SGD(lr=params['learning_rate'], nesterov=params['isNadam'])###Nadams
    model.compile(loss='categorical_crossentropy',optimizer=sgd,metrics=['accuracy'])

    nb_epoch = params['num_epochs']
    history = model.fit(x_train, y_train, nb_epoch=nb_epoch, validation_data=(x_val, y_val), batch_size=params['batch_size'], verbose=0)
    
    return history, model



def xavier_init(size):
    in_dim = size[0]
    xavier_stddev = 1. / tf.sqrt(in_dim / 2.)
    return tf.random_normal(shape=size, stddev=xavier_stddev)

def auc_multi_class(labels, preds):
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    prc_auc = dict()
    f1_scores = dict()
    
    n_classes = labels.shape[1]
    #####from score vector to binary vector
    pred_label = np.array([np.argmax(x) for x in preds])
    onehot_pred_label = []
    for j in pred_label:
        ar = [0 for _ in range(n_classes)]
        ar[j] = 1
        onehot_pred_label.append(ar)
    onehot_pred_label = np.array(onehot_pred_label)
    

    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(labels[:, i], preds[:, i])
        cancer_type = Data.from_int_to_label[i]
        roc_auc[cancer_type] = auc(fpr[i], tpr[i])
        prc_auc[cancer_type] = average_precision_score(labels[:, i], preds[:, i])
        f1_scores[cancer_type] = f1_score(labels[:, i], onehot_pred_label[:, i])
    return roc_auc, prc_auc, f1_scores


def ExternalValidate(test_x, test_y, data_set, keras_model):

    ###############load the trained NN

    pred_prob = keras_model.predict(test_x)
    pred_labels = np.argmax(pred_prob, axis=1)

    pred_labels = [Data.from_int_to_label[x] for x in pred_labels]

    pred_labels = adjust_label(pred_labels, data_set)
    #print(list(test_y))
    #print(list(pred_labels))
    cm_labels = list(set(pred_labels) | set(test_y))
    cm_labels.sort()

    _confusion_matrix = confusion_matrix(test_y, pred_labels, labels = cm_labels)
    _confusion_matrix_df = pd.DataFrame(_confusion_matrix, index = cm_labels, columns=cm_labels)
    _confusion_matrix_df.to_csv('Confusion_matrix_'+data_set+'_'+Config.treat+'.csv')

    f = open('By_Class_metric_'+data_set+'_'+Config.treat+'.txt', 'w+')
    report = classification_report(test_y, pred_labels, labels = cm_labels)
    f.write(report)
    f.close()

    acc = np.sum(_confusion_matrix.diagonal()) / np.sum(_confusion_matrix)
    print('Overall '+ data_set + ' accuracy: {} %'.format(acc*100))
    
    #top 5 accuracy
    if data_set != 'PDX':
        real_int = [Data.from_label_to_int[x] for x in test_y]
        pred_probs = keras_model.predict(test_x)
        top5_acc = top_n_accuracy(real_int, pred_probs, n=5)
        print('Overall top 5 accuracy  '+data_set+' : {} %'.format(top5_acc*100))


  
        
def reset_weights(model):
    from keras import backend as K
    session = K.get_session()
    for layer in model.layers: 
        if hasattr(layer, 'kernel_initializer'):
            layer.kernel.initializer.run(session=session)

def cross_validate_keras(keras_model, best_talos_params, split_size, train_x_all, train_y_all):
    
    kf = KFold(n_splits=split_size, shuffle = True, random_state = Config.random_seed)
    
    reals = []
    real_bins = []
    preds = []
    pred_probs = []

    
    for train_idx, val_idx in kf.split(train_x_all, train_y_all):
        reset_weights(keras_model)
        
        train_x = train_x_all[train_idx]
        train_y = train_y_all[train_idx]
        val_x = train_x_all[val_idx]
        val_y = train_y_all[val_idx]
        
        sgd = SGD(lr=best_talos_params['learning_rate'], nesterov=best_talos_params['isNadam'])###Nadams
        keras_model.compile(loss='categorical_crossentropy',optimizer=sgd,metrics=['accuracy'])
        keras_model.fit(train_x, train_y, nb_epoch=best_talos_params['num_epochs'], validation_data=(val_x, val_y), batch_size=best_talos_params['batch_size'], verbose=0)
        #pred_val = session.run(pred, feed_dict={x: val_x, y: val_y, is_training_tensor: False})
        pred_prob = keras_model.predict(val_x)
        pred_labels = np.argmax(pred_prob, axis=1)

        pred_labels = [Data.from_int_to_label[x] for x in pred_labels]
        
        true_labels = np.argmax(val_y, axis=1)
        true_labels = [Data.from_int_to_label[x] for x in true_labels]
        reals += list(true_labels)
        preds += list(pred_labels)
        
        real_bins += list(val_y)
        pred_probs += list(pred_prob)
    
        #######auroc calculation using sklearn, for each class
    
    reals = np.array(reals)
    preds = np.array(preds)
    
    cm_labels = list(set(reals) | set(preds))
    cm_labels.sort()
    
    _confusion_matrix = confusion_matrix(reals, preds, labels = cm_labels)
    _confusion_matrix_df = pd.DataFrame(_confusion_matrix, index = cm_labels, columns=cm_labels)
    _confusion_matrix_df.to_csv('Confusion_matrix_CV_'+Config.treat+'.csv')

    #Model.plot_confusion_matrix(tr_labels, pred_labels, list(v for v in Data.from_int_to_label.values()), file="Confusion_matrix_DNN_CV_test.pdf")
    print('Confusion_matrix_CV_'+Config.treat+'.csv'," CV Test confusion matrix created:")
    
    acc = np.sum(_confusion_matrix.diagonal()) / np.sum(_confusion_matrix)
    print('Overall CV accuracy: {} %'.format(acc*100))
    
    ##########top 5 accuracy
    real_int = np.argmax(real_bins, axis=1)
    top5_acc = top_n_accuracy(real_int, pred_probs, n=5)
    print('Overall top 5 accuracy: {} %'.format(top5_acc*100))



# first we have to make sure to input data and params into the function
def cnn_1d_model(x_train, y_train, x_val, y_val, params):

    # next we can build the model exactly like we would normally do it
    x_train = np.expand_dims(x_train, axis=2)
    x_val = np.expand_dims(x_val, axis=2)
    
    model = Sequential()
    model.add(Convolution1D(nb_filter=params['number_filter'], 
                            filter_length=params['filter_length'], 
                            input_shape=(x_train.shape[1],1),
                            kernel_initializer=params['kernel_initializer']
                            ))
    model.add(MaxPooling1D(pool_size=params['max_pooling_size']))
    model.add(Activation(params['activation_func']))
    model.add(Flatten())
    model.add(Dropout(params['dropout']))
    hidden_layers(model, params, Data.num_class)
    model.add(Dense(Data.num_class))
    model.add(Activation('softmax'))

    sgd = SGD(lr=params['learning_rate'], nesterov=True)###Nadams
    model.compile(loss='categorical_crossentropy',optimizer=sgd,metrics=['accuracy'])

    nb_epoch = params['num_epochs']
    history = model.fit(x_train, y_train, nb_epoch=nb_epoch, validation_data=(x_val, y_val), batch_size=params['batch_size'])
    
    return history, model

def cross_validate(keras_model, split_size=5):
    
    kf = KFold(n_splits=split_size, shuffle = True, random_state = Config.random_seed)
    
    reals = []
    preds = []
    
    for train_idx, val_idx in kf.split(train_x_all, train_y_all):
        train_x = train_x_all[train_idx]
        train_y = train_y_all[train_idx]
        val_x = train_x_all[val_idx]
        val_y = train_y_all[val_idx]
        run_train(session, train_x, train_y, num_epochs)
        #pred_val = session.run(pred, feed_dict={x: val_x, y: val_y, is_training_tensor: False})
        pred_val = session.run(pred, feed_dict={x: val_x, y: val_y})
        
        reals += list(val_y)
        preds += list(pred_val)
    
        #######auroc calculation using sklearn, for each class
    
    reals = np.array(reals)
    preds = np.array(preds)
    auroc, auprc, f1Score = Model.auc_multi_class(reals, preds)
    ######one hot labels
    pred_labels = np.array([np.argmax(x) for x in preds])
    onehot_pred_labels = []
    
    for j in pred_labels:
        ar = [0 for _ in range(num_class)]
        ar[j] = 1
        onehot_pred_labels.append(ar)
    onehot_pred_labels = np.array(onehot_pred_labels)
    
    accuracy = accuracy_score(reals, onehot_pred_labels)
    
    tr_labels = np.array([np.argmax(x) for x in reals])
    
    _string_labels = [Data.from_int_to_label[x] for x in range(num_class)]
    _confusion_matrix_df = pd.DataFrame(confusion_matrix(tr_labels, pred_labels, labels = range(num_class)), index = _string_labels, columns=_string_labels)
    _confusion_matrix_df.to_csv('Confusion_matrix_DNN_254_CV_'+Config.treat+'.csv')
    #Model.plot_confusion_matrix(tr_labels, pred_labels, list(v for v in Data.from_int_to_label.values()), file="Confusion_matrix_DNN_CV_test.pdf")
    print("CV Test confusion matrix created:")
    
    return accuracy, auroc, auprc, f1Score


def xavier_init(size):
    in_dim = size[0]
    xavier_stddev = 1. / tf.sqrt(in_dim / 2.)
    return tf.random_normal(shape=size, stddev=xavier_stddev)



def plot_confusion_matrix(y_true, y_pred, classes,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues,
                          file='Confustion_matrix.pdf'):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    #"""
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred, range(len(classes)))
    # Only use the labels that appear in the data
    print(classes)
#    classes = [classes[i] for i in unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')
    ####fixing label tick position
    ax.set_xticks(ax.get_xticks()[::2])
    
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    fig.savefig(file,figsize=(16, 12), dpi=300)
    return ax



