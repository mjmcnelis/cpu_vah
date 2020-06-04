cd python
rm -r model_parameters
mkdir model_parameters
python3 sample_model_parameters.py $1

# $1 = number of samples you want to generate
# note: if no argument is given, default number of samples = 10
