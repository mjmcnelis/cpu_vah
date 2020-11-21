import sys
import numpy as np
import os

print('\n\n--------------------------------------------------\n\n')
print('Processing training data...\n')

if len(sys.argv) > 3:
	hydro_mode = str(sys.argv[1])
	samples = int(sys.argv[2])
	events = int(sys.argv[3])
else:
	print('process_train_data error: need 3 command arguments (e.g. vah 1000 1)\n')
	quit()

# it's possible that fireball_radius_1.dat doesn't exist (from failed single-event) but unlikely

parameters = np.loadtxt('../tests/auto_grid/train_data/' + hydro_mode + '/model_parameters/model_parameters_1.dat').reshape(1, -1)
fireball   = np.loadtxt('../tests/auto_grid/train_data/' + hydro_mode + '/fireball_radius/fireball_radius_1.dat').reshape(-1, 5)   # 5 columns by default
train_data = np.append(parameters, [fireball[:,4].mean(), fireball[:,4].std()]).reshape(1, -1)

total_successful_events = fireball.shape[0]


for i in range(1, samples):
	parameters = np.loadtxt('../tests/auto_grid/train_data/' + hydro_mode + '/model_parameters/model_parameters_' + str(i + 1) + '.dat').reshape(1,-1)

	fireball_dir = '../tests/auto_grid/train_data/' + hydro_mode + '/fireball_radius/fireball_radius_' + str(i + 1) + '.dat'

	if os.path.exists(fireball_dir):
		fireball = np.loadtxt(fireball_dir).reshape(-1, 5)

		total_successful_events += fireball.shape[0]

		train_sample = np.append(parameters, [fireball[:,4].mean(), fireball[:,4].std()]).reshape(1, -1)
		train_data   = np.concatenate((train_data, train_sample))


print('Hydrodynamic simulation mode = ', hydro_mode)
print('Generated parameter samples  = ', samples)
print()
print('Total number of events ran  = ', samples * events)
print('Number of successful events = ', total_successful_events)
print()
print('Training data shape =', train_data.shape)

np.savetxt('../tests/auto_grid/train_data/' + hydro_mode + '/train_data.dat', train_data)




