rm -r output
mkdir output
mkdir output/benchmarks
mkdir output/fireball_radius
mkdir output/adaptive
rm -r semi_analytic
mkdir semi_analytic

cp omp/yes_omp/OpenMP.h rhic/include

make clean
make

qsub hydro_omp.pbs
