import numpy as np
import math
import os
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import Lasso
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold


# root mean squared error
def RMSE(y, y_predict):
    return math.sqrt(np.sum((y - y_predict)**2)/len(y))


# auto grid formula in cpu_vah
def auto_grid_size(mean, std, number_sigmas, margin):
    return 2 * (mean  +  number_sigmas * std  + margin)


# auto grid example parameters
number_sigmas = 2    # number of sigmas
margin        = 1    # additional margin [fm] on each side of fireball


print("\nSearching for optimized cubic polynomial Lasso regression models for fireball radius mean and std...\n")


# load training data
train_data = np.loadtxt('train/train_data.dat')
print("training samples =", train_data.shape[0])


# polynomial order
order = 3


# Lasso regression regularization parameter values
points = 101
alpha = np.linspace(0.001, 0.011, points)

RMSE_test_min = 1000.
alpha_min  = 0


# number of folds
k = 5
kfold = KFold(n_splits = k, shuffle = True, random_state = 440)


# radius mean (loop over alpha values)
for a in alpha:

    X = train_data[:,0:15].reshape(-1,15)
    y = train_data[:,-2].reshape(-1,1)

    RMSE_test_sum = 0

    # k-fold validation
    for train_index, test_index in kfold.split(X):
        X_train = X[train_index]
        y_train = y[train_index]

        X_test  = X[test_index]
        y_test  = y[test_index]

        regression_mean = Pipeline([('scale', StandardScaler()),
                                    ('poly',  PolynomialFeatures(order)),
                                    ('reg',   Lasso(alpha = a))])

        regression_mean.fit(X_train, y_train)

        y_test_predict = regression_mean.predict(X_test).reshape(-1,1)

        RMSE_test_sum += RMSE(y_test, y_test_predict)

    RMSE_test = RMSE_test_sum / k

    if RMSE_test < RMSE_test_min:
        RMSE_test_min = RMSE_test
        alpha_min = a


mean_radius_RMSE = RMSE_test_min
mean_radius_alpha = alpha_min

print("\nmean_test_RMSE = %.4f for alpha = %.4f" % (mean_radius_RMSE, mean_radius_alpha))

RMSE_test_min = 1000.


# radius std (loop over alpha values)
for a in alpha:

    X = train_data[:,0:15].reshape(-1,15)
    y = train_data[:,-1].reshape(-1,1)

    RMSE_test_sum = 0

    # k-fold validation
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


print("\nstd_test_RMSE  = %.4f for alpha = %.4f" % (std_radius_RMSE, std_radius_alpha))


# train regression model for mean fireball radius
alpha = mean_radius_alpha

X_train = train_data[:,0:15].reshape(-1,15)
y_train = train_data[:,-2].reshape(-1,1)

regression_mean = Pipeline([('scale', StandardScaler()),
                            ('poly', PolynomialFeatures(order)),
                            ('reg', Lasso(alpha = alpha))])

regression_mean.fit(X_train, y_train)


# train regression model for std fireball radius
alpha = std_radius_alpha

X_train = train_data[:,0:15].reshape(-1,15)
y_train = train_data[:,-1].reshape(-1,1)

regression_std = Pipeline([('scale', StandardScaler()),
                           ('poly', PolynomialFeatures(order)),
                           ('reg', Lasso(alpha = alpha))])

regression_std.fit(X_train, y_train)


# test auto grid on fireball radius data from launch fixed grid
number_files = len(os.listdir('launch_fixed_grid/model_parameters'))

total   = 0
success = 0
average_area = 0
average_margin = 0

for n in range (0, number_files):
    parameters = np.loadtxt('launch_fixed_grid/model_parameters/model_parameters_' + str(n + 1) + '.dat').reshape(1,-1)

    mean = regression_mean.predict(parameters)[0] + mean_radius_RMSE
    std  = regression_std.predict(parameters)[0]  + std_radius_RMSE

    L = auto_grid_size(mean, std, number_sigmas, margin)
    average_area += L * L

    fireball_radius = np.loadtxt('launch_fixed_grid/fireball_radius/fireball_radius_' + str(n + 1) + '.dat')[:,-1]

    for radius in fireball_radius:
        total += 1

        if L > 2*radius:
            success += 1
            average_margin += (L - 2*radius) / 2


average_area /= number_files
average_margin /= success
success_rate = 100. * success / total

print("\nAuto grid with %.1f sigmas and %.1f fm extra margin was able to fit %ld / %ld fireball samples\n"
      % (number_sigmas, margin, success, total))
print("Success rate   = %.2f"    % success_rate, "%")
print("Average margin = %.2f fm" % average_margin)
print()
print("Average auto grid area = %.0f fm^2" % average_area)
print("Large fixed grid area  = 900 fm^2")


np.savetxt('train/mean_radius_alpha_RMSE.dat', [[mean_radius_RMSE, mean_radius_alpha]])
np.savetxt('train/std_radius_alpha_RMSE.dat',  [[std_radius_RMSE,  std_radius_alpha]])




