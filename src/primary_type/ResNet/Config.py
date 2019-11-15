#################################config#######################################
import sys
import argparse
##################parsing args from command line
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--data_set', help='data set name')
#parser.add_argument('--model_name', help='ckpt file name in tmp folder')
#parser.add_argument('--alpha', help='L1 loss coefficient')
parser.add_argument('--treat', help='special treatment for the run')
#parser.add_argument('--batch_size', help='set the batch size')
#parser.add_argument('--num_epochs', help='set the epoch number')
#parser.add_argument('--keep_prob', help='set the drop out rate')
#parser.add_argument('--nn_structure', help='set nn structure')



args = parser.parse_args()

if args.data_set:
    dataSet = args.data_set
else:
    print('Please add dataset name as --data_set=')
    #running on laptop
    os.chdir('/Users/c-zhaoy/NetBeansProjects/JaxHelixCompsci/SubtypeClassifier')

    #sys.exit()


if args.treat:
    treat = args.treat
else:
    print('Please add treatment label as --treat=')
    sys.exit()

#batch size has to be larger than 200, since the class least population is 50
#among 10000 training data samples in order all classes are covered in each randomly sampled
#batch

learning_rate = 0.00001 #initial learning rate


split_rate = 0.2 #test data ratio in the whole data

#for i in range(num_layer):
#    nn_structures.append(str(num_feature//(4**i)))

random_seed = 1024
fold_number = 10



