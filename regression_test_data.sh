rm *pbs.o*
rm -r output
mkdir output
rm -r semi_analytic
mkdir semi_analytic

sh load_modules.sh
sh makefile.sh icpc

cd python
rm -r random_model_parameters
cp -r test_model_parameters random_model_parameters
cd ..

make clean
make

# fixed number of test parameters samples = 250

for((n = 1; n <= 250; n++))
do
   qsub -v JOBID="$n",EVENTS="$1" hydro.pbs
done

# $1 = number of hydro events per job
#
#
# benchmark info for walltime needed in hydro.pbs:
#
# 151 x 151 x 1 hydro grid takes ~ 90s maximum per event
# 301 x 301 x 1 hydro grid takes ~ 8x longer