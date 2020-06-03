import numpy as np
import os

print('\nProcessing train data...')

samples = len(os.listdir('train/model_parameters'))

parameters = np.loadtxt('train/model_parameters/model_parameters_1.dat')
fireball   = np.loadtxt('train/fireball_radius/fireball_radius_1.dat')

train_data = np.append(parameters, [fireball[:,4].mean(), fireball[:,4].std()]).reshape(1,-1)

for n in range(1, samples):
    parameters = np.loadtxt('train/model_parameters/model_parameters_' + str(n + 1) + '.dat')
    fireball   = np.loadtxt('train/fireball_radius/fireball_radius_' + str(n + 1) + '.dat')

    train_sample = np.append(parameters, [fireball[:,4].mean(), fireball[:,4].std()]).reshape(1,-1)
    train_data   = np.concatenate((train_data, train_sample))

print('\ntrain data shape =', train_data.shape)

# save train data
np.savetxt('train/train_data.dat', train_data)