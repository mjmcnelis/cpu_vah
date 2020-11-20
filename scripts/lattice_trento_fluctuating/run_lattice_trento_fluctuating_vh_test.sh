echo
echo 'Testing nonconformal quasiparticle viscous hydrodynamics for single fluctuating TRENTo event...'
echo

cd ../..

cp tests/lattice_trento_fluctuating/initial_profile/e_block.dat tables
cp tests/lattice_trento_fluctuating/parameters/vh/hydro.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vh/lattice.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vh/initial.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vh/Macros.h rhic/include
cp tests/lattice_trento_fluctuating/parameters/vh/OpenMP.h rhic/include

sh clear_results.sh

make clean
make

cd jobs/lattice_trento_fluctuating

qsub run_lattice_trento_fluctuating_vh_test.pbs
