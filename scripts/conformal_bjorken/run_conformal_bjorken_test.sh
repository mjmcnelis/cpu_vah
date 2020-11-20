echo
echo 'Testing conformal anisotropic hydrodynamics for Bjorken flow...'
echo

cd ../..

cp tests/conformal_bjorken/parameters/hydro.properties parameters
cp tests/conformal_bjorken/parameters/lattice.properties parameters
cp tests/conformal_bjorken/parameters/initial.properties parameters
cp tests/conformal_bjorken/parameters/Macros.h rhic/include
cp tests/conformal_bjorken/parameters/OpenMP.h rhic/include

rm -r tests/conformal_bjorken/data
mkdir tests/conformal_bjorken/data

sh hydro.sh 1
cp -r semi_analytic tests/conformal_bjorken/data
cp -r output tests/conformal_bjorken/data/vah
