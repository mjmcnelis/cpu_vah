echo
echo 'Generating training data to optimize grid for 2+1d nonconformal hydrodynamics with fluctuating Trento initial conditions'
echo
echo 'Hydro mode =' $1
echo
echo 'Model parameter samples =' $2
echo
echo 'Hydro events per job =' $3
echo

cd ..

sh sample_model_parameters.sh $2

cd ../..

cp tests/auto_grid/generate_training_data_fluctuating/parameters/$1/hydro.properties parameters
cp tests/auto_grid/generate_training_data_fluctuating/parameters/$1/lattice.properties parameters
cp tests/auto_grid/generate_training_data_fluctuating/parameters/$1/initial.properties parameters
cp tests/auto_grid/generate_training_data_fluctuating/parameters/$1/Macros.h rhic/include

sh clear_results.sh

make clean
make

cd jobs/auto_grid
rm *.pbs.o*

for((n = 1; n <= $2; n++))
do
    qsub -v JOBID="$n",EVENTS="$3" run_fluctuating_training_events.pbs
done


# $1 = hydro mode                        (vah, vh, vh2)
# $2 = number of model parameter samples (each job runs their own parameter sample)
# $3 = number of events per job          (e.g. 400 fluctuating events)
#
# e.g. sh generate_training_data_fluctuating.sh vah 1000 400
