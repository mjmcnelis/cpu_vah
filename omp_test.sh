rm *pbs.o*
rm -r output
mkdir output
mkdir output/benchmarks
mkdir output/fireball_radius
rm -r semi_analytic
mkdir semi_analytic

make clean
make

qsub hydro_omp.pbs
