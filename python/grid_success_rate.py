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
#
def RMSE(y, y_predict):
    return math.sqrt(np.sum((y - y_predict)**2)/len(y))


# transverse auto grid length
#
def auto_grid_size(mean, std, sigma_factor, margin):
    return 2 * (mean  +  sigma_factor * std  +  margin)


# command arguments
#
if len(sys.argv) > 4:
    hydro_mode = str(sys.argv[1])
    samples = int(sys.argv[2])
    margin = float(sys.argv[3])
    sigma_factor = float(sys.argv[4])
else:
    print('grid_success_rate error: need 4 command arguments (e.g. vah 200 2.5 0)\n')
    quit()


print('\n--------------------------------------------------\n\n')
print('Testing auto grid\'s fireball fit success rate...\n\n')
print('Hydrodynamic simulation mode =', hydro_mode)
print('Model parameter samples      =', samples)
print()
print('Margin parameter =', margin, 'fm')
print('Sigma factor     =', sigma_factor)
print()



# test auto grid on fireball radius data from fixed grid
#
total = 0
success = 0
average_area = 0
average_margin = 0


# compute fireball fit success rate and grid area reduction
#
for i in range (0, samples):
    # fireball size predictions for automated grid
    #
    predicted_fireball_size = np.loadtxt('fireball_size_predictions/fireball_size_' + str(i + 1) + '.dat').reshape(1, -1)

    mean = predicted_fireball_size[0, 0]
    std  = predicted_fireball_size[0, 1]

    L = auto_grid_size(mean, std, sigma_factor, margin)

    average_area += L * L

    input_dir = '../tests/auto_grid/benchmark_test/fixed_grid/' + hydro_mode + '/fireball_radius/fireball_radius_'+ str(i + 1) + '.dat'

    if os.path.exists(input_dir):
        # fireball radius in fixed grid events
        #
        fixed_grid_fireball_radius = np.loadtxt(input_dir).reshape(-1, 5)[:, -1]

        for radius in fixed_grid_fireball_radius:
            total += 1

            if L > 2*radius:
                success += 1
                average_margin += (L - 2*radius) / 2


average_area /= samples
average_margin /= success
success_rate = 100. * success / total


print("\nAuto grid with %.1f sigmas and %.1f fm extra margin was able to fit %ld / %ld fireball samples\n"
      % (sigma_factor, margin, success, total))
print("Success rate   = %.2f"    % success_rate, "%")
print("Average margin = %.2f fm" % average_margin)
print()
print("Average auto grid area = %.0f fm^2" % average_area)
print("Large fixed grid area  = 900 fm^2")




