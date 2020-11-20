echo
echo 'Testing simulation runtime for smooth 3+1d TRENTo without OpenMP acceleration...'
echo
echo 'Hydro mode =' $1
echo 'Events =' $2
echo

cd ../..

cp tests/lattice_trento_smooth/initial_profile/e_block.dat tables
cp tests/lattice_trento_smooth/parameters/$1/hydro.properties parameters
cp tests/lattice_trento_smooth/parameters/$1/lattice.properties parameters
cp tests/lattice_trento_smooth/parameters/$1/initial.properties parameters
cp tests/lattice_trento_smooth/parameters/$1/Macros.h rhic/include
cp omp/no_omp/OpenMP.h rhic/include

rm -r output
mkdir output
mkdir output/benchmarks
mkdir output/fireball_radius
mkdir output/adaptive
rm -r semi_analytic
mkdir semi_analytic

make clean
make

cd jobs/omp

for ((n = 1; n <= $2; n++))
do
    qsub run_no_omp_test.pbs
done

# $1 = directory for hydro mode (vah, vh or vh2)
# $2 = number of runs (e.g. 5)


