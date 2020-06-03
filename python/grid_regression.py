import numpy as np
import os
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import Lasso
from sklearn.preprocessing import StandardScaler

# auto grid formula in cpu_vah
def auto_grid_size(mean, std, number_sigmas, margin):
    return 2 * (mean  +  number_sigmas * std  + margin)


# auto grid example parameters
number_sigmas = 2    # number of sigmas
margin        = 1    # additional margin [fm] on each side of fireball


print("\nTraining optimized Lasso regression models for fireball radius mean and std...\n")


# load training data
train_data = np.loadtxt('train/train_data.dat')
print("training samples =", train_data.shape[0])

[mean_radius_RMSE, mean_radius_alpha] = np.loadtxt('train/mean_radius_alpha_RMSE.dat')
[std_radius_RMSE,  std_radius_alpha]  = np.loadtxt('train/std_radius_alpha_RMSE.dat')


# train regression model for mean fireball radius
order = 3
alpha = mean_radius_alpha

X_train = train_data[:,0:15].reshape(-1,15)
y_train = train_data[:,-2].reshape(-1,1)

regression_mean = Pipeline([('scale', StandardScaler()),
                            ('poly', PolynomialFeatures(order)),
                            ('reg', Lasso(alpha = alpha))])

regression_mean.fit(X_train, y_train)


# train regression model for std fireball radius
order = 3
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
print("Average margin = %.2f fm\n" % average_margin)
print("Average auto grid area = %.0f fm^2" % average_area)
print("Large fixed grid area  = 900 fm^2")



# save model predictions for radius mean, std given random model parameters
for n in range (0,len(os.listdir('model_parameters'))):
    parameters = np.loadtxt('model_parameters/model_parameters_' + str(n + 1) + '.dat').reshape(1,-1)

    # model predictions
    mean = regression_mean.predict(parameters)[0] + mean_radius_RMSE
    std  = regression_std.predict(parameters)[0]  + std_radius_RMSE
    np.savetxt('fireball_size_predictions/fireball_size_' + str(n + 1) + '.dat', np.array([[mean, std]]))






