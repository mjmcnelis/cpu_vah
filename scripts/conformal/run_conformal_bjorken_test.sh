echo
echo 'Testing conformal anisotropic hydrodynamics for Bjorken flow...'
echo

cd ../..

cp tests/conformal_bjorken/vah_files/hydro.properties parameters
cp tests/conformal_bjorken/vah_files/lattice.properties parameters
cp tests/conformal_bjorken/vah_files/initial.properties parameters
cp tests/conformal_bjorken/vah_files/Macros.h rhic/include

rm -r tests/conformal_bjorken/results
mkdir tests/conformal_bjorken/results

sh hydro.sh 1
cp -r semi_analytic tests/conformal_bjorken/results
cp -r output tests/conformal_bjorken/results/vah
