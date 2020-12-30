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
./cpu_vah

cp output/surface_slice_x.dat tests/lattice_trento_fluctuating/freezeout/$1
cp output/surface_slice_z.dat tests/lattice_trento_fluctuating/freezeout/$1

cp output/Rinv_shear_x.dat tests/lattice_trento_fluctuating/freezeout/$1
cp output/Rinv_shear_z.dat tests/lattice_trento_fluctuating/freezeout/$1

cp output/Rinv_bulk_x.dat tests/lattice_trento_fluctuating/freezeout/$1
cp output/Rinv_bulk_z.dat tests/lattice_trento_fluctuating/freezeout/$1

cp output/Rinv_dB_x.dat tests/lattice_trento_fluctuating/freezeout/$1
cp output/Rinv_dB_z.dat tests/lattice_trento_fluctuating/freezeout/$1

cp output/regulate_viscous_x.dat tests/lattice_trento_fluctuating/freezeout/$1
cp output/regulate_viscous_z.dat tests/lattice_trento_fluctuating/freezeout/$1

cp output/regulate_dB_x.dat tests/lattice_trento_fluctuating/freezeout/$1
cp output/regulate_dB_z.dat tests/lattice_trento_fluctuating/freezeout/$1

cd tests/lattice_trento_fluctuating/fireball/$1
make

# $1 = hydro mode (vah, vh, vh2)
