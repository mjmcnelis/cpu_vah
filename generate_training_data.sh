rm *pbs.o*
rm -r output
mkdir output
mkdir output/benchmarks
mkdir output/fireball_radius
rm -r semi_analytic
mkdir semi_analytic

sh makefile.sh icpc

sh generate_parameter_samples.sh $1

make clean
make

for ((n = 1; n <= $1; n++))
do
    qsub -v JOBID="$n",EVENTS="$2" hydro.pbs
done

# $1 = number of jobs (each job runs their own parameter sample)
# $2 = number of hydro events per job
