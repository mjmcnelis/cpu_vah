cp conformal_bjorken/vah_files/hydro.properties ../parameters
cp conformal_bjorken/vah_files/lattice.properties ../parameters
cp conformal_bjorken/vah_files/initial.properties ../parameters
cp conformal_bjorken/vah_files/Macros.h ../rhic/include

rm -r conformal_bjorken/results
mkdir conformal_bjorken/results

cd ..
sh hydro.sh 1
cp -r semi_analytic tests/conformal_bjorken/results
cp -r output tests/conformal_bjorken/results/vah
