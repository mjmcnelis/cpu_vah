echo
echo 'Benchmarking auto grid for 2+1d nonconformal hydrodynamics with fluctuating Trento initial conditions'
echo
echo 'Hydro mode =' $1
echo
echo 'Model parameter samples = 200 (fixed)'
echo
echo 'Hydro events per sample =' $2
echo


echo 'Copying simulation parameters...'

cd ../..
cp tests/auto_grid/benchmark_test/auto_grid/$1/hydro.properties parameters
cp tests/auto_grid/benchmark_test/auto_grid/$1/lattice.properties parameters
cp tests/auto_grid/benchmark_test/auto_grid/$1/initial.properties parameters
cp tests/auto_grid/benchmark_test/auto_grid/$1/Macros.h rhic/include


echo 'Copying test model parameters...'

rm -r python/model_parameters
cp -r tests/auto_grid/benchmark_test/model_parameters python


echo 'Running python...'

cd python
rm -r fireball_size_predictions
mkdir fireball_size_predictions
python3 predict_fireball_radius.py $1 200 $3


# test fit success rate ($4 = margin, $5 = sigma_factor)
#		 smooth: 	margin = 2.5    sigma_factor = 0
#	fluctuating:	margin = 1      sigma_factor = 2
#
# see auto_grid formula in python script:
python3 grid_success_rate.py $1 200 $4 $5


while true; do
    read -p "Do you want to submit batch for auto-grid runtimes?" yn
    case $yn in
        [Yy]* ) echo 'Proceeding with job submissions'; break;;
        [Nn]* ) echo 'Aborting job submissions (edit margin and sigma_factor)'; exit;;
        * ) echo "Please answer yes or no.";;
    esac
done


echo
echo 'Clearing benchmarks and fireball radius data in auto_grid...'

cd ..
rm -r tests/auto_grid/benchmark_test/auto_grid/$1/benchmarks
rm -r tests/auto_grid/benchmark_test/auto_grid/$1/fireball_radius
mkdir tests/auto_grid/benchmark_test/auto_grid/$1/benchmarks
mkdir tests/auto_grid/benchmark_test/auto_grid/$1/fireball_radius


echo 'Compiling and submitting jobs...'

sh clear_results.sh

make clean
make

cd jobs/auto_grid
rm *.pbs.o*

for((n = 1; n <= 200; n++))
do
    qsub -v HYDRO="$1",JOBID="$n",EVENTS="$2" run_auto_grid_events.pbs    # multiply walltime by $3
done

# $1 = hydro mode (vah, vh, vh2)
# $2 = number of events per job for auto grid benchmark test (current job batch)
# $3 = number of events per job for training grid (training job batch)
# $4 = margin parameter (see above)
# $5 = sigma factor (see above)
#
# e.g. sh benchmark_auto_grid.sh vah 50 1 2.5 0




