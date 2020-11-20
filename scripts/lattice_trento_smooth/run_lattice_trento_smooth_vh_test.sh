echo
echo 'Testing 3+1d nonconformal quasiparticle viscous hydrodynamics for smooth TRENTo...'
echo

cd ../..

cp tests/lattice_trento_smooth/initial_profile/e_block.dat tables
cp tests/lattice_trento_smooth/parameters/vh/hydro.properties parameters
cp tests/lattice_trento_smooth/parameters/vh/lattice.properties parameters
cp tests/lattice_trento_smooth/parameters/vh/initial.properties parameters
cp tests/lattice_trento_smooth/parameters/vh/Macros.h rhic/include
cp tests/lattice_trento_smooth/parameters/vh/OpenMP.h rhic/include

sh clear_results.sh

make clean
make

cd jobs/lattice_trento_smooth

qsub run_lattice_trento_smooth_vh_test.pbs
