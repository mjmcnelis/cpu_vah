echo
echo 'Generating freezeout surface and inverse Reynolds number for 2+1d standard viscous hydro...'
echo
echo 'Collision =' $1
echo 'Bulk corrections =' $2
echo

cd ../..

cp tables/energy_density/modified_distribution/$1/e_block.dat tables

# first run to output surface.dat
cp tests/modified_distribution/$1/$2_bulk/parameters/hydro.properties parameters
cp tests/modified_distribution/$1/$2_bulk/parameters/lattice.properties parameters
cp tests/modified_distribution/$1/$2_bulk/parameters/initial.properties parameters
cp tests/modified_distribution/$1/$2_bulk/parameters/Macros.h rhic/include
cp tests/modified_distribution/$1/$2_bulk/parameters/OpenMP.h rhic/include

sh hydro.sh 1

cp output/surface.dat tests/modified_distribution/$1/$2_bulk/data

# rerun to output inverse Reynolds number
cp tests/modified_distribution/$1/$2_bulk/parameters/inverse_reynolds/hydro.properties parameters

sh clear_results.sh
./cpu_vah

cp output/surface_slice_x.dat tests/modified_distribution/$1/$2_bulk/data
cp output/regulate_viscous_x.dat tests/modified_distribution/$1/$2_bulk/data
cp output/Rinv_shear_x.dat tests/modified_distribution/$1/$2_bulk/data
cp output/Rinv_bulk_x.dat tests/modified_distribution/$1/$2_bulk/data

# $1 = (central, noncentral)
# $2 = (small, large)
