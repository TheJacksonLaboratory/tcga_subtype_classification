#!/bin/bash

# set the number of nodes
#SBATCH --ntasks=1

# set the number of nodes
#SBATCH --cpus-per-task=2

#SBATCH --mem=100gb

# set the number of nodes
#SBATCH --gres gpu:1

# set the number of nodes
#SBATCH --gres-flags=enforce-binding

# set max wallclock time
#SBATCH --time=100:00:00

# set name of job
#SBATCH --job-name=1dCNN_10

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=END

# send mail to this address
#SBATCH --mail-user=yue.zhao@jax.org

#SBATCH -o ./output/currentCNN_10

####Hyper parameter
alpha=0.01
dataSet=tumor_type
num_epochs=1000
batch_size=128
random_seed=1024
num_feature=254
nn_structure=200_128_32
dropout_keep_rate=${1}

treat=CNN1d_DEG_Feature_10_smote_OrderChromo

#####scan
python ConvolutionalNeuralNetworks/TaloScan.py --data_set=${dataSet} --treat=${treat}

#####External validation
python ConvolutionalNeuralNetworks/ExternalTest.py --data_set=tumor_type --treat=${treat}










###python3 /projects/compsci/Yue/SubtypeClassifier/DeepNeuralNetworks/CrossValidation.py subtype ${alpha}
