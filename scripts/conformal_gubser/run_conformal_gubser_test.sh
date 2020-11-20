echo
echo 'Testing conformal anisotropic hydrodynamics for Gubser flow...'
echo

cd ../..

cp tests/conformal_gubser/parameters/hydro.properties parameters
cp tests/conformal_gubser/parameters/lattice.properties parameters
cp tests/conformal_gubser/parameters/initial.properties parameters
cp tests/conformal_gubser/parameters/Macros.h rhic/include
cp tests/conformal_bjorken/parameters/OpenMP.h rhic/include

rm -r tests/conformal_gubser/data
mkdir tests/conformal_gubser/data

sh hydro.sh 1
cp -r semi_analytic tests/conformal_gubser/data
cp -r output tests/conformal_gubser/data/vah
