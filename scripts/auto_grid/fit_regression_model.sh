echo
echo 'Fitting regression model to optimize grid for 2+1d nonconformal hydrodynamics with Trento initial conditions...'
echo
echo 'Hydro mode =' $1
echo
echo 'Model parameter samples =' $2
echo
echo 'Hydro events per sample =' $3
echo

cd ../..

# copy results from previous script (generate_training_data) to train_data directory
#
#echo
#echo 'Copying training data...'
#echo

#rm -r tests/auto_grid/train_data/$1/model_parameters
#rm -r tests/auto_grid/train_data/$1/fireball_radius
#
#cp -r python/model_parameters tests/auto_grid/train_data/$1
#cp -r output/fireball_radius tests/auto_grid/train_data/$1


# check if any fireball_radius files are missing (python will skip them below)
#
cd tests/auto_grid/train_data/$1/fireball_radius

seq "$2" | while read -r i
do
    [[ -f "fireball_radius_$i.dat" ]] || echo "fireball_radius_$i.dat is missing"
done


# python results are stored in train_data directory
#
echo
echo 'Running python...'
echo
cd ../../../../../python
python3 process_train_data.py $1 $2 $3
python3 fit_regression_model.py $1 $3


# $1 = hydro simulation mode (vah, vh, vh2)
# #2 = number of jobs (or model parameter samples)
# $3 = number of hydro events per job
#
# e.g. sh fit_regression_model.sh vah 1000 1
#
# note: use the same command arguments as generate_training_data/
