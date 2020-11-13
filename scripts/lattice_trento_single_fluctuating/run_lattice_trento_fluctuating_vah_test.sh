echo
echo 'Testing nonconformal anisotropic hydrodynamics for fluctuating TRENTo...'
echo

cd ../..

rm -r tests/lattice_trento_fluctuating/freezeout/vah
mkdir tests/lattice_trento_fluctuating/freezeout/vah

rm -r tests/lattice_trento_fluctuating/fireball/vah/xy_plane
mkdir tests/lattice_trento_fluctuating/fireball/vah/xy_plane

cp tests/lattice_trento_fluctuating/initial_profile/e_block.dat tables
cp tests/lattice_trento_fluctuating/parameters/vah/hydro.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vah/lattice.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vah/initial.properties parameters
cp tests/lattice_trento_fluctuating/parameters/vah/Macros.h rhic/include
cp omp/yes_omp/OpenMP.h rhic/include

rm -r output
mkdir output
mkdir output/benchmarks
mkdir output/fireball_radius
mkdir output/adaptive
rm -r semi_analytic
mkdir semi_analytic

make clean
make

cd jobs/lattice_trento_fluctuating

qsub run_lattice_trento_fluctuating_vah_test.pbs
