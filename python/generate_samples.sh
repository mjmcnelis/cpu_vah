rm -r random_model_parameters
mkdir random_model_parameters
python3 sample_model_parameters.py $1

# $1 = number of samples you want to generate
# note: if no arguments given, default number of samples = 10
