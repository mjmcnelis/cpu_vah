echo
echo 'Testing vh lattice trento with fluctuating energy density profile...'
echo

cp lattice_trento_fluctuating/initial_profile/e_block.dat ../tables

cp lattice_trento_fluctuating/parameters/vh/hydro.properties ../parameters
cp lattice_trento_fluctuating/parameters/vh/lattice.properties ../parameters
cp lattice_trento_fluctuating/parameters/vh/initial.properties ../parameters
cp lattice_trento_fluctuating/parameters/vh/Macros.h ../rhic/include
cp lattice_trento_fluctuating/parameters/vh/OpenMP.h ../rhic/include

cd ..

sh hydro.sh 1

cp output/surface_slice_x.dat tests/lattice_trento_fluctuating/freezeout/vh
cp output/surface_slice_z.dat tests/lattice_trento_fluctuating/freezeout/vh
cp output/Rinv_shear_x.dat tests/lattice_trento_fluctuating/freezeout/vh
cp output/Rinv_shear_z.dat tests/lattice_trento_fluctuating/freezeout/vh
cp output/Rinv_bulk_x.dat tests/lattice_trento_fluctuating/freezeout/vh
cp output/Rinv_bulk_z.dat tests/lattice_trento_fluctuating/freezeout/vh

cd tests/lattice_trento_fluctuating/fireball/vh
make
