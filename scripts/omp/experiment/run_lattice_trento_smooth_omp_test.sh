

if[ $1 == "vah" ]
then
    echo 'Testing 3+1d nonconformal anisotropic hydrodynamics with smooth TRENTo and OpenMP...'
elif[ $1 == "vh" ]
then
    echo 'Testing 3+1d nonconformal quasiparticle viscous hydrodynamics with smooth TRENTo and OpenMP...']
elif[ $1 == "vh2" ]
then
    echo 'Testing 3+1d nonconformal standard viscous hydrodynamics with smooth TRENTo and OpenMP...']
else
    echo 'Error: select vah, vh or vh2 as first command argument'
fi
echo



#
#cd ../..
#
#cp tests/lattice_trento_smooth/initial_profile/e_block.dat tables
#cp tests/lattice_trento_smooth/parameters/$hydro_mode/hydro.properties parameters
#cp tests/lattice_trento_smooth/parameters/$hydro_mode/lattice.properties parameters
#cp tests/lattice_trento_smooth/parameters/$hydro_mode/initial.properties parameters
#cp tests/lattice_trento_smooth/parameters/$hydro_mode/Macros.h rhic/include
#cp omp/yes_omp/OpenMP.h rhic/include
#
#rm -r output
#mkdir output
#mkdir output/benchmarks
#mkdir output/fireball_radius
#mkdir output/adaptive
#rm -r semi_analytic
#mkdir semi_analytic
#
#make clean
#make
#
#cd jobs
#
#for ((n = 1; n <= $events; n++))
#do
#    qsub run_omp_28_test.pbs
#    qsub run_omp_20_test.pbs
#    qsub run_omp_10_test.pbs
#    qsub run_omp_5_test.pbs
#    qsub run_omp_2_test.pbs
#done
#
## $1 = hydro mode (vah, vh, vh2)
## $2 = number of omp events
#
#
