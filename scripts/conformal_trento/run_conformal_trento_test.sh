echo
echo 'Testing conformal viscous hydrodynamics for smooth TRENTo...'
echo
echo 'Hydro mode = ' $1

cd ../..

cp tables/energy_density/conformal_trento/e_block.dat tables
cp tests/conformal_trento/parameters/$1/hydro.properties parameters
cp tests/conformal_trento/parameters/$1/lattice.properties parameters
cp tests/conformal_trento/parameters/$1/initial.properties parameters
cp tests/conformal_trento/parameters/$1/Macros.h rhic/include
cp tests/conformal_trento/parameters/$1/OpenMP.h rhic/include

sh clear_results.sh

make clean
make

cd jobs/conformal_trento

qsub run_conformal_trento_vh_test.pbs


# $1 = hydro mode (vah, vh)
