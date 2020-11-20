echo
echo 'Testing simulation runtime for smooth 3+1d TRENTo with OpenMP acceleration...'
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
cp omp/yes_omp/OpenMP.h rhic/include

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
    qsub run_omp_28_test.pbs
    qsub run_omp_20_test.pbs
    qsub run_omp_10_test.pbs
    qsub run_omp_5_test.pbs
    qsub run_omp_2_test.pbs
done

# $1 = sub-directory for hydro mode (vah, vh or vh2)
# $2 = number of omp runs (e.g. 5)


