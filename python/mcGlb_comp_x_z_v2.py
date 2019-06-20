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
    return squeeze(data[xp,0,:])

def load_var_int(dir, t, var, nx, ny, nz, xp):
    x, y, n, dataRaw = np.loadtxt(dir + '/' + var + '_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    data = squeeze(np.reshape(dataRaw, (nx, ny, nz)))
    return squeeze(data[:,xp])

if __name__ == '__main__':
    root = '/media/bazow/Data/fluid_dynamic_output_for_thesis/'

    ### mc glauber
    icTag = 'mcGlb'
    # Ideal
    type = 'gpu-vh/3d/mcGlb/conformalEOS/ideal/'
    dataDirIdeal = root + type + 'theta1p1'
    # VH
    type = 'gpu-vh/3d/mcGlb/conformalEOS/shear/'
    dataDirVH = root + type + 'reg1-10_theta1p1'
    # VAH
    type = 'cpu-vah/3d/mcGlb/conformalEOS/shear/'
    dataDirVAH = root + type + 'reg1-10_theta1p1'

    plotDir = 'tests/figs/qgp/mcGlb_x-z'

    nx = 121
    ny = 1
    nz = 121
    dx = 0.1
    dy = 0.1
    dz = 0.1
    xmax = (nx - 1) / 2 * dx
    zmax = (nz - 1) / 2 * dz
    x = np.linspace(-xmax, xmax, num=nx)
    z = np.linspace(-zmax, zmax, num=nz)

    nxx = 10*nx
    xx = np.linspace(-xmax, xmax, num=nxx)

#    xp = 60 # x=0
#    xp = 90 # x=3

    xpList = [60,90]

    t = 1.0
    t2 = t*t

    hbarc = 0.197327

    ############################################
    ### t=1 ####################################
    ############################################
    zF01 = [-4.51637, -4.42572]
    zF02 = [4.5164, 4.42574]

    ############################################
    ### t=3 ####################################
    ############################################
#    zF01 = [-4.69657, -4.51751]
#    zF02 = [4.69532, 4.51661]


    nt=2

    eIdeal = np.zeros((nt, nz))
    eVH = np.zeros((nt, nz))
    eVAH = np.zeros((nt, nz))
    # un
    unIdeal = np.zeros((nt, nz))
    unVH = np.zeros((nt, nz))
    unVAH = np.zeros((nt, nz))
    # u
    utVH = np.zeros((nt, nz))
    uxVH = np.zeros((nt, nz))
    uyVH = np.zeros((nt, nz))
    utVAH = np.zeros((nt, nz))
    uxVAH = np.zeros((nt, nz))
    uyVAH = np.zeros((nt, nz))
    # p
    pIdeal = np.zeros((nt, nz))
    pVH = np.zeros((nt, nz))
    pVAH = np.zeros((nt, nz))
    # pl
    plVAH = np.zeros((nt, nz))
    # pinn
    pinnVH = np.zeros((nt, nz))
    pinnVAH = np.zeros((nt, nz))
    # reg
    regVH = np.zeros((nt, nz))
    regVAH = np.zeros((nt, nz))

    # ptHat
    ptHatVAH = np.zeros((nt, nz))
    # pratio
    pratioVH = np.zeros((nt, nz))
    pratioVAH = np.zeros((nt, nz))

    for i in range(0,2):
        xp = xpList[i]
        # ed
        eIdeal[i,:] = load_var_int(dataDirIdeal, t, 'e', nx, ny, nz, xp)
        eVH[i,:] = load_var_int(dataDirVH, t, 'e', nx, ny, nz, xp)
        eVAH[i,:] = load_var_int(dataDirVAH, t, 'e', nx, ny, nz, xp)
        # un
        unIdeal[i,:] = load_var_int(dataDirIdeal, t, 'un', nx, ny, nz, xp)
        unVH[i,:] = load_var_int(dataDirVH, t, 'un', nx, ny, nz, xp)
        unVAH[i,:] = load_var_int(dataDirVAH, t, 'un', nx, ny, nz, xp)
        # u
        utVH[i,:] = load_var_int(dataDirVH, t, 'ut', nx, ny, nz, xp)
        uxVH[i,:] = load_var_int(dataDirVH, t, 'ux', nx, ny, nz, xp)
        uyVH[i,:] = load_var_int(dataDirVH, t, 'uy', nx, ny, nz, xp)
        utVAH[i,:] = load_var_int(dataDirVAH, t, 'ut', nx, ny, nz, xp)
        uxVAH[i,:] = load_var_int(dataDirVAH, t, 'ux', nx, ny, nz, xp)
        uyVAH[i,:] = load_var_int(dataDirVAH, t, 'uy', nx, ny, nz, xp)
        # p
        pIdeal[i,:] = load_var_int(dataDirIdeal, t, 'p', nx, ny, nz, xp)
        pVH[i,:] = load_var_int(dataDirVH, t, 'p', nx, ny, nz, xp)
        pVAH[i,:] = load_var_int(dataDirVAH, t, 'p', nx, ny, nz, xp)
        # pl
        plVAH[i,:] = load_var_int(dataDirVAH, t, 'pl', nx, ny, nz, xp)
        # pinn
        pinnVH[i,:] = load_var_int(dataDirVH, t, 'pinn', nx, ny, nz, xp)*t2
        pinnVAH[i,:] = load_var_int(dataDirVAH, t, 'pinn', nx, ny, nz, xp)*t2
        # reg
        regVH[i,:] = load_var_int(dataDirVH, t, 'regulations', nx, ny, nz, xp)
        regVAH[i,:] = load_var_int(dataDirVAH, t, 'regulations', nx, ny, nz, xp)

    # ptHat
    ptHatVAH = (eVAH-plVAH-Rbar(eVAH,pVAH,plVAH)*(eVAH-3*pVAH))/2
    # pratio
    pratioVH = np.divide(pVH+pinnVH,pVH-pinnVH/2)
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

    ##########################################################
    ### Energy density #######################################
    ##########################################################
    i=0
    xp = xpList[i]
    fig1, ax1 = plt.subplots()
    ax1.plot(z, squeeze(eIdeal[i,:])*hbarc, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax1.plot(z, squeeze(eVH[i,:])*hbarc, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax1.plot(z, squeeze(eVAH[i,:])*hbarc, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax1.axvspan(zF02[i], zmax, alpha=0.25, color=colorFO)
    ax1.axvspan(zF01[i], -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax1.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\mathcal{E}\,\mathrm{[GeV/fm^3]}$')
    plt.title(r'$\vec{r}=($' + '{:.1f}'.format(x[xp]) + ', 0.0) fm',fontsize=16)

    i=1
    xp = xpList[i]
    fig2, ax2 = plt.subplots()
    ax2.plot(z, squeeze(eIdeal[i,:])*hbarc, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax2.plot(z, squeeze(eVH[i,:])*hbarc, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax2.plot(z, squeeze(eVAH[i,:])*hbarc, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax2.axvspan(zF02[i], zmax, alpha=0.25, color=colorFO)
    ax2.axvspan(zF01[i], -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax2.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\mathcal{E}\,\mathrm{[GeV/fm^3]}$')
    plt.title(r'$\vec{r}=($' + '{:.1f}'.format(x[xp]) + ', 0.0) fm',fontsize=16)

    y1 = ax1.get_ylim()
    y2 = ax2.get_ylim()
    minVal = np.min([y1[0],y2[0]])
    maxVal = np.max([y1[1],y2[1]])
    ax1.set_ylim([minVal,maxVal])
    ax2.set_ylim([minVal,maxVal])

    ax1.text(0.5, 0.35,
         r'$\tau=$' + '{:.1f}'.format(t) + r'$\,\mathrm{fm/c}$' + '\n' + r'$\tau_0=0.5\,\mathrm{fm/c}$' + '\n' r'$T_0=0.6\,\mathrm{GeV}$' + '\n' + r'$\eta/s=0.2$',
         ha='center', va='center', transform=ax1.transAxes,
         fontsize=16)
    ax2.text(0.5, 0.35,
         r'$\tau=$' + '{:.1f}'.format(t) + r'$\,\mathrm{fm/c}$' + '\n' + r'$\tau_0=0.5\,\mathrm{fm/c}$' + '\n' r'$T_0=0.6\,\mathrm{GeV}$' + '\n' + r'$\eta/s=0.2$',
         ha='center', va='center', transform=ax2.transAxes,
         fontsize=16)

    fig1.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    fig2.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig1.savefig(plotDir+'/ed_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xpList[0]])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    fig2.savefig(plotDir+'/ed_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xpList[1]])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    ##########################################################
    ### un #######################################
    ##########################################################
    i=0
    xp = xpList[i]
    fig1, ax1 = plt.subplots()
    ax1.plot(z, squeeze(unIdeal[i,:]), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax1.plot(z, squeeze(unVH[i,:]), color='green', linewidth=3.5, linestyle='--', label='VH')
    ax1.plot(z, squeeze(unVAH[i,:]), color='red', linewidth=3.5, linestyle='--', label='VAH')
    plt.legend(loc='best', frameon=False)
    ax1.axvspan(zF02[i], zmax, alpha=0.25, color=colorFO)
    ax1.axvspan(zF01[i], -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax1.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$u^\eta$')

    i=1
    xp = xpList[i]
    fig2, ax2 = plt.subplots()
    ax2.plot(z, squeeze(unIdeal[i,:]), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax2.plot(z, squeeze(unVH[i,:]), color='green', linewidth=3.5, linestyle='--', label='VH')
    ax2.plot(z, squeeze(unVAH[i,:]), color='red', linewidth=3.5, linestyle='--', label='VAH')
    plt.legend(loc='best', frameon=False)
    ax2.axvspan(zF02[i], zmax, alpha=0.25, color=colorFO)
    ax2.axvspan(zF01[i], -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax2.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$u^\eta$')

    y1 = ax1.get_ylim()
    y2 = ax2.get_ylim()
    minVal = np.min([y1[0],y2[0]])
    maxVal = np.max([y1[1],y2[1]])
    ax1.set_ylim([minVal,maxVal])
    ax2.set_ylim([minVal,maxVal])

    fig1.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    fig1.savefig(plotDir+'/un_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xpList[0]])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig2.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    fig2.savefig(plotDir+'/un_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xpList[1]])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    ##########################################################
    ### pl/pt #######################################
    ##########################################################
    i=0
    xp = xpList[i]
    fig1, ax1 = plt.subplots()
    ax1.plot(z, np.ones((nz, 1)), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax1.plot(z, squeeze(pratioVH[i,:]), color='green', linewidth=3.5, linestyle='--', label='VH')
    ax1.plot(z, squeeze(pratioVAH[i,:]), color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax1.axvspan(zF02[i], zmax, alpha=0.25, color=colorFO)
    ax1.axvspan(zF01[i], -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax1.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\mathcal{P}_L/\mathcal{P}_\perp$')

    i=1
    xp = xpList[i]
    fig2, ax2 = plt.subplots()
    ax2.plot(z, np.ones((nz, 1)), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax2.plot(z, squeeze(pratioVH[i,:]), color='green', linewidth=3.5, linestyle='--', label='VH')
    ax2.plot(z, squeeze(pratioVAH[i,:]), color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax2.axvspan(zF02[i], zmax, alpha=0.25, color=colorFO)
    ax2.axvspan(zF01[i], -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax2.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\mathcal{P}_L/\mathcal{P}_\perp$')

    y1 = ax1.get_ylim()
    y2 = ax2.get_ylim()
    minVal = np.min([y1[0],y2[0]])
    maxVal = np.max([y1[1],y2[1]])
    ax1.set_ylim([minVal,maxVal])
    ax2.set_ylim([minVal,maxVal])

    fig1.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    fig1.savefig(plotDir+'/pratio_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xpList[0]])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig2.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    fig2.savefig(plotDir+'/pratio_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xpList[1]])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    ##########################################################
    ### reg #######################################
    ##########################################################
    i=0
    xp = xpList[i]
    fig1, ax1 = plt.subplots()
    ax1.plot(z, np.ones((nz, 1)), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax1.plot(z, squeeze(regVH[i,:]), color='green', linewidth=3.5, linestyle='--', label='VH')
    ax1.plot(z, squeeze(regVAH[i,:]), color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax1.axvspan(zF02[i], zmax, alpha=0.25, color=colorFO)
    ax1.axvspan(zF01[i], -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax1.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\mathrm{tanh}\rho/\rho$')

    i=1
    xp = xpList[i]
    fig2, ax2 = plt.subplots()
    ax2.plot(z, np.ones((nz, 1)), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax2.plot(z, squeeze(regVH[i,:]), color='green', linewidth=3.5, linestyle='--', label='VH')
    ax2.plot(z, squeeze(regVAH[i,:]), color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax2.axvspan(zF02[i], zmax, alpha=0.25, color=colorFO)
    ax2.axvspan(zF01[i], -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax2.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\mathrm{tanh}\rho/\rho$')

    y1 = ax1.get_ylim()
    y2 = ax2.get_ylim()
    minVal = np.min([y1[0],y2[0]])
    maxVal = np.max([y1[1],y2[1]])
    ax1.set_ylim([minVal,maxVal])
    ax2.set_ylim([minVal,maxVal])

    fig1.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    fig1.savefig(plotDir+'/reg_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xpList[0]])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig2.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    fig2.savefig(plotDir+'/reg_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xpList[1]])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    plt.show()
