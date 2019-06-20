#!/usr/bin/env python3
import brewer2mpl
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.ndimage
import matplotlib.pyplot as plt
import pylab
from scipy import interpolate
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import colors as mcolors

def Power(a, b):
    return pow(a, b)

def Rbar(e, p, pl):
    a = pl/e
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
    return (0.015267955823446243 + 7.725572805021035*a + 421.0063884634789*a2 + 3422.877939650926*a3 - 5785.670846299543*a4 - 12261.66452089229*a5 +
     31491.409484673808*a6 - 22737.05146992673*a7 + 5441.373392185447*a8)/ (
        0.05470696094814806 + 14.505878005231883*a + 522.6643024173569*a2 + 2731.7776413939037*a3 - 6161.1991042880445*a4 -
     3989.4375208972588*a5 + 15008.260526258282*a6 - 10243.036679405379*a7 + 2116.74060159494*a8)

def load_var(dir, t, var, nx, ny, nz):
    x, y, n, dataRaw = np.loadtxt(dir + '/' + var + '_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    data = np.reshape(dataRaw, (nx, ny, nz))
    return squeeze(data)

def load_var_int(dir, t, var, nx, ny, nz, xx):
    x, y, n, dataRaw = np.loadtxt(dir + '/' + var + '_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    data = np.reshape(dataRaw, (nx, ny, nz))
    dataInt = interpolate.interp1d(x, squeeze(data), kind='cubic')
    return dataInt(xx)

def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]

if __name__ == '__main__':
    root = '/media/bazow/Data/fluid_dynamic_output_for_thesis/'
#    '''
    ### mc glauber
    icTag = 'mcGlb'
    # Ideal
    type = 'gpu-vh/1d/' + icTag + '/conformalEOS/ideal/'
    dataDirIdeal = root + type + 'theta1p1'
    # VH
    type = 'gpu-vh/1d/' + icTag + '/conformalEOS/shear/'
    dataDirVH = root + type + 'reg1-10_theta1p1'
    # AH
    type = 'cpu-vah/1d/' + icTag + '/conformalEOS/'
    dataDirAH = root + type + 'LO_reg1-10_theta1p1'
    # VAH
    type = 'cpu-vah/1d/' + icTag + '/conformalEOS/'
    dataDirVAH = root + type + 'reg1-10_theta1p1'

    plotDir = 'tests/figs/qgp/' + icTag + ''
#    '''
    '''
    ### continous glauber
    icTag = 'Glb'
    # Ideal
    type = 'gpu-vh/1d/' + icTag + '/conformalEOS/ideal/'
    dataDirIdeal = '/media/bazow/Data/fluid_dynamic_output_for_thesis/gpu-vh/Glb_1d_ideal_isotropicIC_conformalEOS'
    # VH
    type = 'gpu-vh/1d/' + icTag + '/conformalEOS/shear/'
    dataDirVH = '/media/bazow/Data/fluid_dynamic_output_for_thesis/gpu-vh/Glb_1d_shear_isotropicIC_conformalEOS_etaOverS0p2_reg1-10'
    # AH
    type = 'cpu-vah/1d/' + icTag + '/conformalEOS/'
    dataDirAH = '/media/bazow/Data/fluid_dynamic_output_for_thesis/cpu-vah/Glb_1d_shear_LOAhydro_isotropicIC_conformalEOS_etaOverS0p2'
    # VAH
    type = 'cpu-vah/1d/' + icTag + '/conformalEOS/'
    dataDirVAH = '/media/bazow/Data/fluid_dynamic_output_for_thesis/cpu-vah/Glb_1d_shear_withPimunu_WTz_isotropicIC_conformalEOS_etaOverS0p2_reg1-10'

    plotDir = 'tests/figs/qgp/' + icTag + ''
    '''

    nx = 201
    ny = 1
    nz = 1
    dx = 0.1
    dy = 0.1
    dz = 0.1
    xmax = (nx - 1) / 2 * dx
    x = np.linspace(-xmax, xmax, num=nx)

    nxx = 10*nx
    xx = np.linspace(-xmax, xmax, num=nxx)

    t = 1.0
    t2 = t*t

    hbarc = 0.197327

    # For t=1 -- Glb
    xF01 = -4.6723
    xF02 = 4.6723
    # For t=3 -- Glb
    xF01 = -4.32059
    xF02 = 4.32059

    # For t=1 -- mcGlb
    xF01 = -5.57862
    xF02 = 5.4526
    # For dynamical Kn freezeout
    xKn1 = [-4.74794, 4.60536]
    xKn0p7 = [-4.37841, 4.24168]
    xKn0p5 = [-3.1386, 2.84833]
    '''
    # For t=3 -- mcGlb
    xF01 = -6.34775
    xF02 = 6.04202
    # For dynamical Kn freezeout
    xKn1 = [-5.55553, 5.27722]
    xKn0p7 = [-4.46538, 4.27531]
    xKn0p5 = [-3.51163, 3.5734]
    #'''
    # ed
    eIdeal = load_var_int(dataDirIdeal, t, 'e', nx, ny, nz, xx)
    eVH = load_var_int(dataDirVH, t, 'e', nx, ny, nz, xx)
    eAH = load_var_int(dataDirAH, t, 'e', nx, ny, nz, xx)
    eVAH = load_var_int(dataDirVAH, t, 'e', nx, ny, nz, xx)
    # ux
    uxIdeal = load_var_int(dataDirIdeal, t, 'ux', nx, ny, nz, xx)
    uxVH = load_var_int(dataDirVH, t, 'ux', nx, ny, nz, xx)
    uxAH = load_var_int(dataDirAH, t, 'ux', nx, ny, nz, xx)
    uxVAH = load_var_int(dataDirVAH, t, 'ux', nx, ny, nz, xx)
    # p
    pIdeal = load_var_int(dataDirIdeal, t, 'p', nx, ny, nz, xx)
    pVH = load_var_int(dataDirVH, t, 'p', nx, ny, nz, xx)
    pAH = load_var_int(dataDirAH, t, 'p', nx, ny, nz, xx)
    pVAH = load_var_int(dataDirVAH, t, 'p', nx, ny, nz, xx)

    # theta
    thetaVH = load_var_int(dataDirVH, t, 'theta', nx, ny, nz, xx)
    thetaAH = load_var_int(dataDirAH, t, 'theta', nx, ny, nz, xx)
    thetaVAH = load_var_int(dataDirVAH, t, 'theta', nx, ny, nz, xx)
    # taupi
    taupiVH = load_var_int(dataDirVH, t, 'taupi', nx, ny, nz, xx)
    taupiAH = load_var_int(dataDirAH, t, 'taupi', nx, ny, nz, xx)
    taupiVAH = load_var_int(dataDirVAH, t, 'taupi', nx, ny, nz, xx)

    KnTaupiVAH = load_var_int(dataDirVAH, t, 'knTaupi', nx, ny, nz, xx)

    inds = indices(KnTaupiVAH, lambda x: x > 1.0)
    print(inds)

    #####################################################################################################
    # Plots
    #####################################################################################################
    plt.style.use('fivethirtyeight')
    plt.style.use('seaborn-whitegrid')
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

    print(plt.style.available)

    minorLocator = MultipleLocator(1)

    pad=0.5
    h_pad = None
    w_pad = None
    rect = [0, 0, 1, 1]

    colorFO = '#8da0cb'
    colorFO = 'dodgerblue'

    # Valid names are: ['Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3']
    bmap = brewer2mpl.get_map('Set2', 'qualitative', 4)
    colors = bmap.mpl_colors

    # KnTaupi
    fig, ax = plt.subplots()
#    ax.plot(xx, np.ones((nxx,1)), color='black', linewidth=1.5, linestyle='-')
    ax.plot(xx, np.multiply(taupiVH,np.abs(thetaVH)), color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, np.multiply(taupiAH,np.abs(thetaAH)), color='blue', linewidth=3.5, linestyle='--', label='AH')
    ax.plot(xx, np.abs(KnTaupiVAH), color='red', linewidth=3.5, linestyle='--', label='VAH')
    plt.legend(loc='best', frameon=False)
    ax.axvspan(xF02, xmax, alpha=0.25, color=colors[2], label=r'$\mathrm{T}\geq T_f$')
    ax.axvspan(xF01, -xmax, alpha=0.25, color=colors[2])
#    ax.axvspan(xF02, xKn1[1], alpha=0.25, color='#8da0cb')
    ax.axvspan(xF01, xKn1[0], alpha=0.25, color=colors[0], label=r'$\mathrm{Kn}\geq 1$')
    ax.axvspan(xKn1[1], xF02, alpha=0.25, color=colors[0])
    ax.axvspan(xKn1[0], xKn0p7[0], alpha=0.25, color=colors[1], label=r'$\mathrm{Kn}\geq 0.7$')
    ax.axvspan(xKn0p7[1], xKn1[1], alpha=0.25, color=colors[1])
    ax.axvspan(xKn0p7[0], xKn0p5[0], alpha=0.25, color=colors[3], label=r'$\mathrm{Kn}\geq 0.5$')
    ax.axvspan(xKn0p5[1], xKn0p7[1], alpha=0.25, color=colors[3])
#    pylab.xlim([-xF01,xF02])
    pylab.ylim([0,2])
    pylab.xlim([-xmax,xmax])
    plt.legend(loc='best', frameon=False)
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathrm{Kn}\equiv\theta\tau_\pi$')
    plt.tight_layout(pad=0.5, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/Kn_dynamicalFO_t-' + '{:.1f}'.format(t) + '.pdf', pad=0.5, h_pad=h_pad, w_pad=w_pad, rect=rect)

    plt.show()
