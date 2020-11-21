import numpy as np
import math
import os
import sys
from sklearn.externals import joblib
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import Lasso
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold


# root mean squared error
def RMSE(y, y_predict):
    return math.sqrt(np.sum((y - y_predict)**2)/len(y))


# command arguments
if len(sys.argv) > 2:
    hydro_mode = str(sys.argv[1])
    events = int(sys.argv[2])
else:
    print('fit_regression_model error: need 2 command arguments (e.g. vah 1)\n')
    quit()


print('\n\n--------------------------------------------------\n\n')
print('Fitting regression model...\n\n')
print('Hydrodynamic simulation mode   =', hydro_mode)
print('Number of hydro events per job =', events)
print()

include_fireball_radius_std = (events > 1)

# smooth Trento event has fireball_radius_std = 0

if include_fireball_radius_std:
    print('X = [b, P_B]                                        (feature matrix)')
    print('Y = [fireball_radius_mean, fireball_radius_std]     (target)')
else:
    print('X = [b, P_B]                   (feature matrix)')
    print('Y = [fireball_radius_mean]     (target)')


print('\nLoading training data...\n')
train_data = np.loadtxt('../tests/auto_grid/train_data/' + hydro_mode + '/train_data.dat')

print('Training data shape =', train_data.shape)




# hard coded number of feature and target variables (adjust for future 3+1d auto grid algorithms)
features = 15
targets  = 2

if(features + targets != train_data.shape[1]):
    print("\nfit_regression_model error: check number of selected features and targets in training data\n")




print("\nOptimizing regression model...\n")

# Lasso regression polynomial order
order = 3

# Lasso regression regularization alpha (or lambda)
points = 11
alpha = np.linspace(0.001, 0.021, points)

# number of folds for cross validation
k = 5
kfold = KFold(n_splits = k, shuffle = True, random_state = 440)

RMSE_test_min = 1000
alpha_min  = 0



# Optimize Lasso regression model for fireball_radius_mean (loop over alpha values)
for a in alpha:
    X = train_data[:,0:features].reshape(-1, features)
    y = train_data[:,-2].reshape(-1, 1)     # assumes second-last column is mean transverse fireball radius

    RMSE_test_sum = 0

    # cross validation
    for train_index, test_index in kfold.split(X):
        X_train = X[train_index]
        y_train = y[train_index]

        X_test  = X[test_index]
        y_test  = y[test_index]

        regression_mean = Pipeline([('scale', StandardScaler()),
                                    ('poly',  PolynomialFeatures(order)),
                                    ('reg',   Lasso(alpha = a))])

        regression_mean.fit(X_train, y_train)

        y_test_predict = regression_mean.predict(X_test).reshape(-1, 1)

        RMSE_test_sum += RMSE(y_test, y_test_predict)

    RMSE_test = RMSE_test_sum / k

    if RMSE_test < RMSE_test_min:
        RMSE_test_min = RMSE_test
        alpha_min = a



mean_radius_RMSE = RMSE_test_min
mean_radius_alpha = alpha_min

print("\nfireball_radius_mean RMSE = %.4f     (alpha = %.4f)" % (mean_radius_RMSE, mean_radius_alpha))



std_radius_RMSE = 0
std_radius_alpha = 0

# only train regression model for fireball_radius_std if used fluctuating Trento events
if include_fireball_radius_std:
    RMSE_test_min = 1000                        # reset value for alpha search


    # Optimize Lasso regression model for fireball_radius_std (loop over alpha values)
    for a in alpha:
        X = train_data[:,0:features].reshape(-1, features)
        y = train_data[:,-1].reshape(-1, 1)     # assumes last column is std transverse fireball radius

        RMSE_test_sum = 0

        # cross validation
        for train_index, test_index in kfold.split(X):
            X_train = X[train_index]
            y_train = y[train_index]

            X_test  = X[test_index]
            y_test  = y[test_index]

            regression_std = Pipeline([('scale', StandardScaler()),
                                      ('poly', PolynomialFeatures(order)),
                                      ('reg', Lasso(alpha = a))])

            regression_std.fit(X_train, y_train)

            y_test_predict = regression_std.predict(X_test).reshape(-1,1)

            RMSE_test_sum += RMSE(y_test, y_test_predict)

        RMSE_test = RMSE_test_sum / k

        if RMSE_test < RMSE_test_min:
            RMSE_test_min = RMSE_test
            alpha_min = a


    std_radius_RMSE = RMSE_test_min
    std_radius_alpha = alpha_min


print("fireball_radius_std  RMSE = %.4f     (alpha = %.4f)" % (std_radius_RMSE, std_radius_alpha))



# train regression model for fireball_radius_mean (using all the training data)
alpha = mean_radius_alpha

X_train = train_data[:,0:features].reshape(-1, features)
y_train = train_data[:,-2].reshape(-1, 1)       # assumes second-last column is mean transverse fireball radius

regression_mean = Pipeline([('scale', StandardScaler()),
                            ('poly', PolynomialFeatures(order)),
                            ('reg', Lasso(alpha = alpha))])

regression_mean.fit(X_train, y_train)



# only train regression model for fireball_radius_std if used fluctuating Trento events
regression_std = Pipeline([('scale', StandardScaler()),
                            ('poly', PolynomialFeatures(order)),
                            ('reg', Lasso(alpha = alpha))])

# train regression model for fireball_radius_std (using all the training data)
if include_fireball_radius_std:
    alpha = std_radius_alpha

    X_train = train_data[:,0:features].reshape(-1, features)
    y_train = train_data[:,-1].reshape(-1, 1)       # assumes last column is std transverse fireball radius

    regression_std.fit(X_train, y_train)


print('\nSaving RMSE and regression model(s)...\n')

target_dir = '../tests/auto_grid/regression_model/' + hydro_mode

np.savetxt(target_dir + '/mean_radius_alpha_RMSE.dat',  [[mean_radius_RMSE, mean_radius_alpha]])
np.savetxt(target_dir + '/std_radius_alpha_RMSE.dat',   [[std_radius_RMSE,  std_radius_alpha]])

joblib.dump(regression_mean, target_dir + '/regression_fireball_radius_mean.pkl')

if include_fireball_radius_std:
    joblib.dump(regression_std,  target_dir + '/regression_fireball_radius_std.pkl')






