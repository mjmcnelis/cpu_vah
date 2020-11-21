#!/usr/bin/env python3
import sys
import random as r
import numpy as np
import time
from datetime import datetime

# Uniform model parameter sampling
# Linear impact parameter sampling

def sample_impact_parameter(R):
    while True:
        b = r.uniform(0, 2*R)
        w = b / (2*R)

        if(r.random() < w):
            return b

samples = 10                       # default number of parameter samples
R       = 7.0                      # assumed radius of Pb nucleus [fm]

if len(sys.argv) > 1:
    samples = int(sys.argv[1])     # overwrite with command argument

for i in range(0, samples):
    b = sample_impact_parameter(R) # impact parameter [fm]
    N = r.uniform(10, 20)          # normalization Pb+Pb [GeV]
    p = r.uniform(-0.7, 0.7)       # generalized mean
    w = r.uniform(0.5, 1.5)        # nucleon width [fm]
    dmin = r.uniform(0, 1.7)       # minimum nucleon - nucleon distance [fm]
    sigmak = r.uniform(0.3, 2.0)   # multiplicity fluctuation (gamma distribution std)
    Tsw = r.uniform(0.135, 0.165)  # switching temperature [GeV]
    Tk = r.uniform(0.13, 0.3)      # etas kink temperature [GeV]
    etask = r.uniform(0.01, 0.2)   # etas value at Tk
    aL = r.uniform(-2, 1)          # left etas slope [GeV^-1]
    aH = r.uniform(-1, 2)          # right etas slope [GeV^-1]
    zetasN = r.uniform(0.01, 0.25) # max value of zetas
    Tp = r.uniform(0.12, 0.3)      # zetas peak temperature [GeV]
    wz = r.uniform(0.025, 0.15)    # zetas width [GeV]
    lambdaz = r.uniform(-0.8, 0.8) # zetas skewness

    parameters = np.array([[b, N, p, w, dmin, sigmak, Tsw, Tk, etask, aL, aH, zetasN, Tp, wz, lambdaz]])

    np.savetxt('model_parameters/model_parameters_' + str(i + 1) + '.dat', parameters)

