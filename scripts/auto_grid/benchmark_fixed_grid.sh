echo
echo 'Benchmarking fixed grid for 2+1d nonconformal hydrodynamics with fluctuating Trento initial conditions'
echo
echo 'Hydro mode =' $1
echo
echo 'Model parameter samples = 200 (fixed)'
echo
echo 'Hydro events per sample' $2
echo


# copy simulation parameters to parameters/ and rhic/include
#
cd ../..
cp tests/auto_grid/benchmark_test/fixed_grid/$1/hydro.properties parameters
cp tests/auto_grid/benchmark_test/fixed_grid/$1/lattice.properties parameters
cp tests/auto_grid/benchmark_test/fixed_grid/$1/initial.properties parameters
cp tests/auto_grid/benchmark_test/fixed_grid/$1/Macros.h rhic/include


# copy test model_parameters to python/
#
rm -r python/model_parameters
cp -r tests/auto_grid/benchmark_test/model_parameters python


# clear benchmarks and fireball_radius directories in fixed_grid
#
rm -r tests/auto_grid/benchmark_test/fixed_grid/$1/benchmarks
rm -r tests/auto_grid/benchmark_test/fixed_grid/$1/fireball_radius
mkdir tests/auto_grid/benchmark_test/fixed_grid/$1/benchmarks
mkdir tests/auto_grid/benchmark_test/fixed_grid/$1/fireball_radius


# compile and submit jobs
#
sh clear_results.sh

make clean
make

cd jobs/auto_grid
rm *.pbs.o*

for((n = 1; n <= 200; n++))
do
    qsub -v HYDRO="$1",JOBID="$n",EVENTS="$2" run_fixed_grid_events.pbs    # multiply walltime by $3
done

# $1 = hydro mode               (vah, vh, vh2)
# $2 = number of events per job (e.g. 50 fluctuating events)
#
# e.g. sh benchmark_fixed_grid.sh vah 50
