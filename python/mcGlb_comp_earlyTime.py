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
from matplotlib.offsetbox import AnchoredText

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

if __name__ == '__main__':
    root = '/media/bazow/Data/fluid_dynamic_output_for_thesis/'
#    '''
    ### mc glauber
    # Ideal
    type = 'gpu-vh/1d/mcGlb/conformalEOS/ideal/'
    dataDirIdeal = root + type + 't0-0.25'
    # VH
    type = 'gpu-vh/1d/mcGlb/conformalEOS/shear/'
    dataDirVH = root + type + 't0-0p25_xi0-10/reg1-1_theta1p1'
    #
    type = 'gpu-vh/1d/mcGlb/conformalEOS/shear/'
    dataDirVH_noReg = root + type + 't0-0p25_xi0-10/reg100-100_theta1p1'
    # VAH
    type = 'cpu-vah/1d/mcGlb/conformalEOS/'
    dataDirVAH = root + type + 't0-0.25_xi0-100_reg1-10_theta1p1'

    plotDir = 'tests/figs/qgp/mcGlb_earlyTime'

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

    t = 0.5
    t=3.0
    t2 = t*t

    hbarc = 0.197327

    # For t=0.5 -- mcGlb
    xF01 = -5.54508
    xF02 = 5.42785

    e = np.zeros((1, nx))
    pl = np.zeros((1, nx))

    # ed
    eIdeal = load_var_int(dataDirIdeal, t, 'e', nx, ny, nz, xx)
    eVH = load_var_int(dataDirVH, t, 'e', nx, ny, nz, xx)
    eVH_noReg = load_var_int(dataDirVH_noReg, t, 'e', nx, ny, nz, xx)
    eVAH = load_var_int(dataDirVAH, t, 'e', nx, ny, nz, xx)
    # ux
    uxIdeal = load_var_int(dataDirIdeal, t, 'ux', nx, ny, nz, xx)
    uxVH = load_var_int(dataDirVH, t, 'ux', nx, ny, nz, xx)
    uxVH_noReg = load_var_int(dataDirVH_noReg, t, 'ux', nx, ny, nz, xx)
    uxVAH = load_var_int(dataDirVAH, t, 'ux', nx, ny, nz, xx)
    # p
    pIdeal = load_var_int(dataDirIdeal, t, 'p', nx, ny, nz, xx)
    pVH = load_var_int(dataDirVH, t, 'p', nx, ny, nz, xx)
    pVH_noReg = load_var_int(dataDirVH_noReg, t, 'p', nx, ny, nz, xx)
    pVAH = load_var_int(dataDirVAH, t, 'p', nx, ny, nz, xx)
    # pl
    plVAH = load_var_int(dataDirVAH, t, 'pl', nx, ny, nz, xx)
    # pinn
    pinnVH = load_var_int(dataDirVH, t, 'pinn', nx, ny, nz, xx)*t2
    pinnVH_noReg = load_var_int(dataDirVH_noReg, t, 'pinn', nx, ny, nz, xx)*t2
    pinnVAH = load_var_int(dataDirVAH, t, 'pinn', nx, ny, nz, xx)*t2
    # reg
    regVH = load_var(dataDirVH, t, 'regulations', nx, ny, nz)
    regVH_noReg = load_var(dataDirVH_noReg, t, 'regulations', nx, ny, nz)
    regVAH = load_var(dataDirVAH, t, 'regulations', nx, ny, nz)

    # ptHat
    ptHatVAH = (eVAH-plVAH-Rbar(eVAH,pVAH,plVAH)*(eVAH-3*pVAH))/2
    # pratio
    pratioVH = np.divide(pVH+pinnVH,pVH-pinnVH/2)
    pratioVH_noReg = np.divide(pVH_noReg+pinnVH_noReg,pVH_noReg-pinnVH_noReg/2)
    pratioVAH = np.divide(plVAH,ptHatVAH)

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

    fig, ax = plt.subplots()
    ax.plot(xx, eIdeal*hbarc, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(xx, eVH*hbarc, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, eVAH*hbarc, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax.plot(xx, eVH_noReg*hbarc, color='blue', linewidth=3.5, linestyle=':', label='VH-(100,100)')
    ax.axvspan(xF02, xmax, alpha=0.25, color=colorFO)
    ax.axvspan(xF01, -xmax, alpha=0.25, color=colorFO)
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathcal{E}\,\mathrm{[GeV/fm^3]}$')
    text(0.1, 0.9, '(a)', ha='center', va='center', transform=ax.transAxes,fontsize=16)
    text(0.8, 0.8,
         r'$\tau=$' + '{:.1f}'.format(t) + r'$\,\mathrm{fm/c}$' + '\n'
         + r'$\tau_0=0.25\,\mathrm{fm/c}$' + '\n' r'$T_0=0.6\,\mathrm{GeV}$' + '\n' + r'$\eta/s=0.2$',
         ha='center', va='center', transform=ax.transAxes,
         fontsize=16)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/ed_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(xx, uxIdeal, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(xx, uxVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, uxVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax.plot(xx, uxVH_noReg, color='blue', linewidth=3.5, linestyle=':', label='VH-(10,100)')
    plt.legend(loc='best', frameon=False)
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$u^x$')
    text(0.1, 0.9, '(b)', ha='center', va='center', transform=ax.transAxes,fontsize=16)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/ux_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(xx, np.ones((nxx, 1)), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(xx, pratioVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, pratioVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax.plot(xx, pratioVH_noReg, color='blue', linewidth=3.5, linestyle=':', label='VH-(100,100)')
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathcal{P}_L/\mathcal{P}_\perp$')
    text(0.08, 0.9, '(c)', ha='center', va='center', transform=ax.transAxes,fontsize=16)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/pl_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(x, np.ones((nx, 1)), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(x, regVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(x, regVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax.plot(x, regVH_noReg, color='blue', linewidth=3.5, linestyle=':', label='VH-(100,100)')
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathrm{tanh}\rho/\rho$')
    text(0.1, 0.1, '(d)', ha='center', va='center', transform=ax.transAxes,fontsize=16)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/reg_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    plt.show()
