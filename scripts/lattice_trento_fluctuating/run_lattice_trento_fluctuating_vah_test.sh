echo
echo 'Testing nonconformal anisotropic hydrodynamics for single fluctuating TRENTo event...'
echo

cd ../..

cp tests/lattice_trento_fluctuating/initial_profile/e_block.dat tables
cp tests/lattice_trento_fluctuating/parameters/vah/hydro.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vah/lattice.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vah/initial.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vah/Macros.h rhic/include
cp tests/lattice_trento_fluctuating/parameters/vah/OpenMP.h rhic/include

sh clear_results.sh

make clean
make

cd jobs/lattice_trento_fluctuating

qsub run_lattice_trento_fluctuating_vah_test.pbs
