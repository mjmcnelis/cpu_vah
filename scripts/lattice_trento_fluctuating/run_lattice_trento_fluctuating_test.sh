echo
echo 'Testing nonconformal hydrodynamics for single fluctuating TRENTo event...'
echo
echo 'Hydro mode = ' $1
echo

cd ../..

cp tables/energy_density/lattice_trento/fluctuating/e_block.dat tables

cp tests/lattice_trento_fluctuating/parameters/$1/hydro.properties parameters
cp tests/lattice_trento_fluctuating/parameters/$1/lattice.properties parameters
cp tests/lattice_trento_fluctuating/parameters/$1/initial.properties parameters
cp tests/lattice_trento_fluctuating/parameters/$1/Macros.h rhic/include
cp tests/lattice_trento_fluctuating/parameters/$1/OpenMP.h rhic/include

rm -r tests/lattice_trento_fluctuating/fireball/$1/xy_plane
rm -r tests/lattice_trento_fluctuating/freezeout/$1

mkdir tests/lattice_trento_fluctuating/fireball/$1/xy_plane
mkdir tests/lattice_trento_fluctuating/freezeout/$1

sh clear_results.sh

make clean
make

cd jobs/lattice_trento_fluctuating

qsub -v HYDRO="$1" run_lattice_trento_fluctuating_test.pbs


# $1 = hydro mode (vah, vh, vh2)