echo
echo 'Testing 3+1d nonconformal hydrodynamics for multiple fluctuating TRENTo events...'
echo
echo 'Hydro mode =' $1
echo 'Events =' $2
echo

cd ../..

cp tests/lattice_trento_multiple_events/parameters/$1/hydro.properties parameters
cp tests/lattice_trento_multiple_events/parameters/$1/lattice.properties parameters
cp tests/lattice_trento_multiple_events/parameters/$1/initial.properties parameters
cp tests/lattice_trento_multiple_events/parameters/$1/Macros.h rhic/include
cp tests/lattice_trento_multiple_events/parameters/$1/OpenMP.h rhic/include

sh clear_results.sh

make clean
make

cd jobs/lattice_trento_multiple_events
rm *.pbs.o*

for ((n = 1; n <= $2; n++))
do
    qsub run_lattice_trento_multiple_events_test.pbs
done

# $1 = sub-directory for hydro mode (vah, vh or vh2)
# $2 = number of fluctuating events (e.g. 10)
