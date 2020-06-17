cp conformal_gubser/vah_files/hydro.properties ../parameters
cp conformal_gubser/vah_files/lattice.properties ../parameters
cp conformal_gubser/vah_files/initial.properties ../parameters
cp conformal_gubser/vah_files/Macros.h ../rhic/include

rm -r conformal_gubser/results
mkdir conformal_gubser/results

cd ..
sh hydro.sh 1
cp -r semi_analytic tests/conformal_gubser/results
cp -r output tests/conformal_gubser/results/vah
