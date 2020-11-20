echo
echo 'Testing 3+1d nonconformal anisotropic hydrodynamics for smooth TRENTo...'
echo

cd ../..

cp tests/lattice_trento_smooth/initial_profile/e_block.dat tables
cp tests/lattice_trento_smooth/parameters/vah/hydro.properties parameters
cp tests/lattice_trento_smooth/parameters/vah/lattice.properties parameters
cp tests/lattice_trento_smooth/parameters/vah/initial.properties parameters
cp tests/lattice_trento_smooth/parameters/vah/Macros.h rhic/include
cp tests/lattice_trento_smooth/parameters/vah/OpenMP.h rhic/include

sh clear_results.sh

make clean
make

cd jobs/lattice_trento_smooth

qsub run_lattice_trento_smooth_vah_test.pbs
