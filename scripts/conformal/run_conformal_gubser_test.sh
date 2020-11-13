echo
echo 'Testing conformal anisotropic hydrodynamics for Gubser flow...'
echo

cd ../..

cp tests/conformal_gubser/vah_files/hydro.properties parameters
cp tests/conformal_gubser/vah_files/lattice.properties parameters
cp tests/conformal_gubser/vah_files/initial.properties parameters
cp tests/conformal_gubser/vah_files/Macros.h rhic/include

rm -r tests/conformal_gubser/results
mkdir tests/conformal_gubser/results

sh hydro.sh 1
cp -r semi_analytic tests/conformal_gubser/results
cp -r output tests/conformal_gubser/results/vah
