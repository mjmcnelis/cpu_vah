echo
echo 'Testing 3+1d nonconformal hydrodynamics for smooth TRENTo...'
echo
echo 'Hydro mode = ' $1
echo

cd ../..

cp tables/energy_density/lattice_trento/smooth/e_block.dat tables

cp tests/lattice_trento_smooth/parameters/$1/hydro.properties parameters
cp tests/lattice_trento_smooth/parameters/$1/lattice.properties parameters
cp tests/lattice_trento_smooth/parameters/$1/initial.properties parameters
cp tests/lattice_trento_smooth/parameters/$1/Macros.h rhic/include

rm -r tests/lattice_trento_smooth/data/$1/adaptive
rm -r tests/lattice_trento_smooth/data/$1/x_axis
rm -r tests/lattice_trento_smooth/data/$1/z_axis

mkdir tests/lattice_trento_smooth/data/$1/x_axis
mkdir tests/lattice_trento_smooth/data/$1/z_axis

sh clear_results.sh
make clean
make
./cpu_vah

cp -r output/adaptive tests/lattice_trento_smooth/data/$1

cd tests/lattice_trento_smooth/data/$1
make


# $1 = hydro mode (vah, vh, vh2)
