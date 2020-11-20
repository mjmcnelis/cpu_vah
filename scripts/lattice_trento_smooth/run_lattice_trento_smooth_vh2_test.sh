echo
echo 'Testing 3+1d nonconformal standard viscous hydrodynamics for smooth TRENTo...'
echo

cd ../..

cp tests/lattice_trento_smooth/initial_profile/e_block.dat tables
cp tests/lattice_trento_smooth/parameters/vh2/hydro.properties parameters
cp tests/lattice_trento_smooth/parameters/vh2/lattice.properties parameters
cp tests/lattice_trento_smooth/parameters/vh2/initial.properties parameters
cp tests/lattice_trento_smooth/parameters/vh2/Macros.h rhic/include
cp tests/lattice_trento_smooth/parameters/vh2/OpenMP.h rhic/include

sh clear_results.sh

make clean
make

cd jobs/lattice_trento_smooth

qsub run_lattice_trento_smooth_vh2_test.pbs
