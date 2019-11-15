#!/bin/bash

# set the number of nodes
#SBATCH --ntasks=1

# set the number of nodes
#SBATCH --cpus-per-task=1

#SBATCH --mem=100gb

# set the number of nodes
#SBATCH --gres gpu:1

# set the number of nodes
#SBATCH --gres-flags=enforce-binding

# set max wallclock time
#SBATCH --time=100:00:00

# set name of job
#SBATCH --job-name=Incept_40_deep

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=END

# send mail to this address
#SBATCH --mail-user=yue.zhao@jax.org

#SBATCH -o ./output/currentInception_40_deep

####Hyper parameter
alpha=0.01
dataSet=tumor_type
num_epochs=1000
batch_size=128
random_seed=1024
num_feature=254
nn_structure=200_128_32
dropout_keep_rate=${1}

treat=InceptionNet1d_diverse_DEG_Feature_40_deep_orderChromo_scan

#####scan
python InceptionNet/TaloScan.py --data_set=${dataSet} --treat=${treat}

#####External validation
python InceptionNet/ExternalTest.py --data_set=tumor_type --treat=${treat}










###python3 /projects/compsci/Yue/SubtypeClassifier/DeepNeuralNetworks/CrossValidation.py subtype ${alpha}
