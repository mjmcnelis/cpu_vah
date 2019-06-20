#!/usr/bin/env python3

from scipy import integrate

import numpy
from matplotlib.pylab import *

import equation_of_state as eos
import specific_bulk_viscosity as zetas
#from plot_setup import plt
import StringIO

def load_var(dir, t, var, nx, ny, nz):
    x, y, n, dataRaw = np.loadtxt(dir + '/' + var + '_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    data = np.reshape(dataRaw, (nx, ny, nz))
    return squeeze(data)

def Rtilde(a):
    a2 = a*a
    a3 = a2*a
    a4 = a3*a
    a5 = a4*a
    a6 = a5*a
    a7 = a6*a
    a8 = a7*a
    a9 = a8*a
    a10 = a9*a
    a11 = a10*a
    a12 = a11*a
    a13 = a12*a
    a14 = a13*a
    a15 = a14*a
    return (-6.674731906076046e-6 + 0.004617789933500251*a + 0.7207562721999754*a2 + 9.097427250602184*a3 - 4.475814747302824*a4 - 36.37501529319408*a5 +
     46.868405146729316*a6 - 15.833867583743228*a7)/(
        0.06856675185266 + 2.9181587012768597*a + 11.951184087839218*a2 - 29.708257843442173*a3 - 2.618233802059826*a4 + 34.646239784689065*a5 -
     19.62596366454439*a6 + 2.374808442453899*a7)

def temp(e):
    fac = 13.8997
    return math.pow(e/fac, 0.25)

def eiso(T):
    fac = 13.8997
    return fac*math.pow(T,4.0)

def equationsAH(t, y):
    delta_pipi = 1.33333
    tau_pipi = 1.42857
    delta_PiPi = 0.666667
    lambda_piPi = 1.2
    etabar = 0.2

    # Assign some variables for convenience of notation
    e = y[0]
    pL = y[1]

    a = pL/e

    T = temp(e)
    taupiInv = T / 5 / etabar
    p = e/3

    # Output from ODE function must be a COLUMN vector, with n rows
    n = len(y)  # 2: implies we have two ODEs
    f = np.zeros((n, 1))
    f[0] = -(e + pL) / t
    f[1] = -(pL - e/3) * taupiInv - (3*pL-Rtilde(a)*e) / t
    return f

def equationsVH(t, y):
    delta_pipi = 1.33333
    tau_pipi = 1.42857
    delta_PiPi = 0.666667
    lambda_piPi = 1.2
    etabar = 0.2

    # Assign some variables for convenience of notation
    e = y[0]
    pizz = y[1]

    T = temp(e)
    taupiInv = T / 5 / etabar
    p = e/3
    beta_pi = (e + p) / 5

    # Output from ODE function must be a COLUMN vector, with n rows
    n = len(y)  # 2: implies we have two ODEs
    f = np.zeros((n, 1))
    f[0] = -(e + p - pizz) / t
    f[1] = -pizz * taupiInv + beta_pi * 4. / 3. / t - (tau_pipi / 3 + delta_pipi) * pizz / t
    return f

if __name__ == '__main__':
    rAH = integrate.ode(equationsAH).set_integrator('vode', method='bdf')
    rVH = integrate.ode(equationsVH).set_integrator('vode', method='bdf')

    t0 = 0.5
    t_final = 10.0
    delta_t = 0.05
    num_steps = int(np.floor((t_final - t0) / delta_t) + 1)

    # Set initial condition
    T0 = 3.05
    ed = 93.21038818333332
    pL0 = ed/3
    print ed
    print pL0
    pizz0 = 0
    rAH.set_initial_value([ed, pL0], t0)
    rVH.set_initial_value([ed, pizz0], t0)

    tAH = np.zeros((num_steps, 1))
    tVH = np.zeros((num_steps, 1))
    eAH= np.zeros((num_steps, 1))
    pLAH = np.zeros((num_steps, 1))
    eVH= np.zeros((num_steps, 1))
    pizzVH = np.zeros((num_steps, 1))
    tAH[0] = t0
    eAH[0] = ed
    pLAH[0] = pL0
    tVH[0] = t0
    eVH[0] = ed
    pizzVH[0] = pizz0

    # VH semi-analytic solution
    k = 1
    while rVH.successful() and k < num_steps:
        rVH.integrate(rVH.t + delta_t)
        tVH[k] = rVH.t
        eVH[k] = rVH.y[0]
        pizzVH[k] = rVH.y[1]
        k += 1

    # AH semi-analytic solution
    k = 1
    while rAH.successful() and k < num_steps:
        rAH.integrate(rAH.t + delta_t)
        tAH[k] = rAH.t
        eAH[k] = rAH.y[0]
        pLAH[k] = rAH.y[1]
        k += 1

    ##############################################
    # load numerical solution
    ##############################################
    dataDirAH = '/media/bazow/Data/fluid_dynamic_output_for_thesis/cpu-vah/conformal_bjorken_test_Glb'

    tAH_est = np.zeros((num_steps, 1))
    eAH_est = np.zeros((num_steps, 1))
    plAH_est = np.zeros((num_steps, 1))

    k = 0
    while k < num_steps:
        ti = t0 + k * delta_t
        tAH_est[k] = ti
        eAH_est[k] = load_var(dataDirAH, ti, 'e', 1, 1, 1)
        plAH_est[k] = load_var(dataDirAH, ti, 'pl', 1, 1, 1)
        k += 1
    '''
    dataDir = 'tests/output/bjorken/nonconformal'

    tEst, eEst = np.loadtxt(dataDir+'/e.dat', unpack=True)
    tEst, plEst = np.loadtxt(dataDir+'/pl.dat', unpack=True)
    tEst, piEst = np.loadtxt(dataDir+'/Pi.dat', unpack=True)

    plotDir = 'tests/figs/bjorken/nonconformal'
    '''
    hbarc = 0.197327

    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    labelSize = 16

    #####################################################################################################
    # Plots
    #####################################################################################################
    plt.style.use('fivethirtyeight')
#    mpl.rcParams['font.family'] = 'Ubuntu'
#    plt.rcParams['text.color'] = 'k'
#    plt.rcParams['xtick.color'] = 'k'
#    plt.rcParams['ytick.color'] = 'k'
#    plt.rcParams['axes.labelcolor'] = 'k'
#    plt.rcParams['axes.facecolor'] = 'white'
#    plt.rcParams['axes.grid'] = 'False'

    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.major.size'] = 5.5
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['xtick.minor.size'] = 3.5
    plt.rcParams['xtick.minor.width'] = 1.0
    plt.rcParams['ytick.major.size'] = 5.5
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['font.size'] = 16

    minorLocator = MultipleLocator(1)

    pad=0.5
    h_pad = None
    w_pad = None
    rect = [0, 0, 1, 1]

    fig, ax = plt.subplots()
    ax.plot(tVH, np.divide(eVH/3-pizzVH,eVH/3+pizzVH/2), color='blue', linewidth=3.5, linestyle='--', label='Ideal')
    ax.plot(tAH, np.divide(pLAH,(eAH-pLAH)/2), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(tAH_est, np.divide(plAH_est,(eAH_est-plAH_est)/2), color='red', linewidth=3.5, linestyle='--', label='Ideal')
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathcal{E}\,\mathrm{[GeV/fm^3]}$')
    ax.set_xscale('log')
    ax.set_xticks([0.5, 1, 5, 10])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    plt.show()
