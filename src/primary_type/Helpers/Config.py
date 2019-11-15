#################################config#######################################
import sys
import argparse
##################parsing args from command line
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--data_set', help='data set name')
parser.add_argument('--treat', help='special treatment for the run')




args = parser.parse_args()

if args.data_set:
    dataSet = args.data_set
else:
    print('Please add dataset name as --data_set=')
    #running on laptop


if args.treat:
    treat = args.treat
else:
    print('Please add treatment label as --treat=')
    sys.exit()


learning_rate = 0.00001 #initial learning rate
split_rate = 0.2 #test data ratio in the whole data

random_seed = 1024
fold_number = 10



