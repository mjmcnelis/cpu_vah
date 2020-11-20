cd ../../python

rm -r model_parameters
mkdir model_parameters

python3 sample_model_parameters.py $1

# $1 = number of samples you want to generate
# default samples = 10 if no argument is given
