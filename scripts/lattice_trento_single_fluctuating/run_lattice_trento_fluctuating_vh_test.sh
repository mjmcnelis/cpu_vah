echo
echo 'Testing nonconformal quasiparticle viscous hydrodynamics for fluctuating TRENTo...'
echo

cd ../..

rm -r tests/lattice_trento_fluctuating/freezeout/vh
mkdir tests/lattice_trento_fluctuating/freezeout/vh

rm -r tests/lattice_trento_fluctuating/fireball/vh/xy_plane
mkdir tests/lattice_trento_fluctuating/fireball/vh/xy_plane

cp tests/lattice_trento_fluctuating/initial_profile/e_block.dat tables
cp tests/lattice_trento_fluctuating/parameters/vh/hydro.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vh/lattice.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vh/initial.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vh/Macros.h rhic/include
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

qsub run_lattice_trento_fluctuating_vh_test.pbs
