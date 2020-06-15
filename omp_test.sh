cd pbs
rm *pbs.o*
cd ..
rm -r output
mkdir output
mkdir output/benchmarks
mkdir output/fireball_radius
rm -r semi_analytic
mkdir semi_analytic

make clean
make

cd pbs
qsub hydro_omp_1.pbs
qsub hydro_omp_2.pbs
qsub hydro_omp_5.pbs
qsub hydro_omp_10.pbs
qsub hydro_omp_20.pbs
