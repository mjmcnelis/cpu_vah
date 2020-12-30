echo
echo 'Testing conformal hydrodynamics for smooth TRENTo...'
echo
echo 'Hydro mode = ' $1
echo

cd ../..

cp tables/energy_density/conformal_trento/e_block.dat tables

cp tests/conformal_trento/parameters/$1/hydro.properties parameters
cp tests/conformal_trento/parameters/$1/lattice.properties parameters
cp tests/conformal_trento/parameters/$1/initial.properties parameters
cp tests/conformal_trento/parameters/$1/Macros.h rhic/include
cp tests/conformal_trento/parameters/$1/OpenMP.h rhic/include

rm -r tests/conformal_trento/data/$1/adaptive
rm -r tests/conformal_trento/data/$1/x_axis
rm -r tests/conformal_trento/data/$1/z_axis

mkdir tests/conformal_trento/data/$1/x_axis
mkdir tests/conformal_trento/data/$1/z_axis

sh clear_results.sh
make clean
make

cd jobs/conformal_trento

qsub -v HYDRO="$1" run_conformal_trento_test.pbs


# $1 = hydro mode (vah or vh)
