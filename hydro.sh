rm -r output
mkdir output
mkdir output/benchmarks
mkdir output/fireball_radius
rm -r semi_analytic
mkdir semi_analytic

make clean
make

for ((i = 1; i <= $1; i++))
do
    echo
    echo "Running event $i"
    ./cpu_vah $2
    echo
done

# $1 = number of hydro events to run successively ($1 >= 1)
# $2 = file index in python/model_parameters/ (1 <= $2 <= samples or leave blank)

# e.g. sh hydro.sh 1 6 runs a single event with sample model parameters model_parameters_6.dat
#      sh hydro.sh 1   runs a single event with fixed model parameters given in parameters/

# note 1: parameters/ contains default values for the model parameters
# note 2: model parameters are a subset of all parameters in parameters/
# note 3: also need to define RANDOM_MODEL_PARAMETERS in rhic/include/Macros.h
#         if you want to replace model parameters in parameters/ with a random sample
