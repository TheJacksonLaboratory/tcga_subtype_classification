####needs to start conda environment before pasting; starting conda env in sh wont work?
i=0
for var in 0.5;
do
echo $var
sbatch /projects/compsci/Yue/SubtypeClassifier/InceptionNet/runInceptionNet.sh ${var};
((i++))
done

