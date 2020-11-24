echo
echo 'Testing simulation runtime for smooth 3+1d TRENTo with OpenMP acceleration...'
echo
echo 'Hydro mode =' $1
echo

cd ../..

cp tables/energy_density/lattice_trento/smooth/e_block.dat tables

cp tests/omp_test/parameters/$1/hydro.properties parameters
cp tests/omp_test/parameters/$1/lattice.properties parameters
cp tests/omp_test/parameters/$1/initial.properties parameters
cp tests/omp_test/parameters/$1/Macros.h rhic/include
cp tests/omp_test/parameters/$1/OpenMP.h rhic/include

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

qsub -v HYDRO="$1" run_omp_28_test.pbs
qsub -v HYDRO="$1" run_omp_20_test.pbs
qsub -v HYDRO="$1" run_omp_10_test.pbs
qsub -v HYDRO="$1" run_omp_5_test.pbs
qsub -v HYDRO="$1" run_omp_2_test.pbs
qsub -v HYDRO="$1" run_omp_1_test.pbs


# $1 = hydro mode (vah, vh or vh2)

