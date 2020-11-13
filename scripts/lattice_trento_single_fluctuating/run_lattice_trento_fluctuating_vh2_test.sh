echo
echo 'Testing nonconformal standard viscous hydrodynamics for fluctuating TRENTo...'
echo

cd ../..

rm -r tests/lattice_trento_fluctuating/freezeout/vh2
mkdir tests/lattice_trento_fluctuating/freezeout/vh2

rm -r tests/lattice_trento_fluctuating/fireball/vh2/xy_plane
mkdir tests/lattice_trento_fluctuating/fireball/vh2/xy_plane

cp tests/lattice_trento_fluctuating/initial_profile/e_block.dat tables
cp tests/lattice_trento_fluctuating/parameters/vh2/hydro.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vh2/lattice.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vh2/initial.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vh2/Macros.h rhic/include
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

cd jobs/lattice_trento_fluctuating

qsub run_lattice_trento_fluctuating_vh2_test.pbs
