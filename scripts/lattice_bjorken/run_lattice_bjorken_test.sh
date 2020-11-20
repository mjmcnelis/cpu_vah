echo
echo 'Testing nonconformal hydrodynamics for Bjorken flow...'
echo
echo 'Hydro mode =' $1
echo

cd ../..

cp tests/lattice_bjorken/parameters/$1/hydro.properties parameters
cp tests/lattice_bjorken/parameters/$1/lattice.properties parameters
cp tests/lattice_bjorken/parameters/$1/initial.properties parameters
cp tests/lattice_bjorken/parameters/$1/Macros.h rhic/include
cp tests/lattice_bjorken/parameters/$1/OpenMP.h rhic/include

sh hydro.sh 1
cp -r semi_analytic tests/lattice_bjorken/data/$1
cp -r output tests/lattice_bjorken/data/$1

#$1 = hydro mode (vah, vh or vh2)
