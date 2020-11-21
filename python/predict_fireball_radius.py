import numpy as np
import os
import sys
from sklearn.externals import joblib
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import Lasso
from sklearn.preprocessing import StandardScaler


if len(sys.argv) > 3:
    hydro_mode = str(sys.argv[1])
    samples = int(sys.argv[2])
    events = int(sys.argv[3])
else:
    print('predict_fireball_radius error: need 3 command arguments (e.g. vah 200 1)\n')
    quit()


print('\n--------------------------------------------------\n\n')
print('Predicting fireball radius mean and std...\n\n')
print('Hydro simulation mode    =', hydro_mode)
print('Model parameter samples  =', samples)
print('Number of events per job =', events)
print()


print('Loading regression model(s) and RMSE...\n')

input_dir = '../tests/auto_grid/regression_model/' + hydro_mode

regression_mean = joblib.load(input_dir + '/regression_fireball_radius_mean.pkl')
regression_std = 0

include_fireball_radius_std = (events > 1)

if include_fireball_radius_std:
    regression_std = joblib.load(input_dir + '/regression_fireball_radius_std.pkl')

[mean_radius_RMSE, mean_radius_alpha] = np.loadtxt(input_dir + '/mean_radius_alpha_RMSE.dat')
[std_radius_RMSE,  std_radius_alpha]  = np.loadtxt(input_dir + '/std_radius_alpha_RMSE.dat')


for i in range (0, samples):
    parameters = np.loadtxt('model_parameters/model_parameters_' + str(i + 1) + '.dat').reshape(1,-1)

    # model predictions
    mean = regression_mean.predict(parameters)[0]  +  mean_radius_RMSE
    std  = 0

    if include_fireball_radius_std:
        std  = regression_std.predict(parameters)[0]  +  std_radius_RMSE

    np.savetxt('fireball_size_predictions/fireball_size_' + str(i + 1) + '.dat', np.array([[mean, std]]))


print('Finished fireball radius predictions\n')






