echo
echo 'Predicting fireball radius mean and std for new model parameter samples'
echo
echo 'Hydro mode =' $1
echo
echo 'New parameter samples =' $2
echo

# generate new parameter samples
#
sh sample_model_parameters.sh $2

# make predictions for fireball radius mean and std
cd ../../python
rm -r fireball_size_predictions
mkdir fireball_size_predictions
python3 predict_fireball_radius.py $1 $2 $3

# $1 = hydro mode (vah, vh, vh2)
# $2 = new model parameter samples
# $3 = number of events per job for training grid (training job batch)
#
# e.g. sh predict_fireball_radius.sh vah 100 1
#
# note: use same command arguments as previous scripts (see README)
