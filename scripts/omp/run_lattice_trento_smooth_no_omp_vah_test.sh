echo
echo 'Testing 3+1d nonconformal anisotropic hydrodynamics for smooth TRENTo...'
echo

cd ../..

cp tests/lattice_trento_smooth/initial_profile/e_block.dat tables
cp tests/lattice_trento_smooth/parameters/vah/hydro.properties parameters
cp tests/lattice_trento_smooth/parameters/vah/lattice.properties parameters
cp tests/lattice_trento_smooth/parameters/vah/initial.properties parameters
cp tests/lattice_trento_smooth/parameters/vah/Macros.h rhic/include
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

cd jobs

for ((n = 1; n <= $1; n++))
do
    qsub run_no_omp_test.pbs
done

# $1 = number of omp runs


