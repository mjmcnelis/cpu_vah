echo
echo 'Testing conformal viscous hydrodynamics for smooth TRENTo...'
echo

cd ../..

cp tests/conformal_trento/initial_profile/e_block.dat tables
cp tests/conformal_trento/parameters/vah/hydro.properties parameters
cp tests/conformal_trento/parameters/vah/lattice.properties parameters
cp tests/conformal_trento/parameters/vah/initial.properties parameters
cp tests/conformal_trento/parameters/vah/Macros.h rhic/include
cp tests/conformal_trento/parameters/vah/OpenMP.h rhic/include

sh clear_results.sh

make clean
make

cd jobs/conformal_trento

qsub run_conformal_trento_vh_test.pbs
