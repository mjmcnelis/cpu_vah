#!/usr/bin/env python3
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.ndimage
import matplotlib.pyplot as plt
import pylab
from scipy import interpolate
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

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

    t = 3.0
    t2 = t*t

    hbarc = 0.197327

    # For t=1 -- Glb
    xF01 = -4.6723
    xF02 = 4.6723
    # For t=3 -- Glb
    xF01 = -4.32059
    xF02 = 4.32059

    e = np.zeros((1, nx))
    pl = np.zeros((1, nx))

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
    # pl
    plAH = load_var_int(dataDirAH, t, 'pl', nx, ny, nz, xx)
    plVAH = load_var_int(dataDirVAH, t, 'pl', nx, ny, nz, xx)
    # pinn
    pinnVH = load_var_int(dataDirVH, t, 'pinn', nx, ny, nz, xx)*t2
    pinnVAH = load_var_int(dataDirVAH, t, 'pinn', nx, ny, nz, xx)*t2
    # reg
    regVH = load_var(dataDirVH, t, 'regulations', nx, ny, nz)
    regVAH = load_var(dataDirVAH, t, 'regulations', nx, ny, nz)

    # theta
    thetaIdeal = load_var_int(dataDirIdeal, t, 'theta', nx, ny, nz, xx)
    thetaVH = load_var_int(dataDirVH, t, 'theta', nx, ny, nz, xx)
    thetaAH = load_var_int(dataDirAH, t, 'theta', nx, ny, nz, xx)
    thetaVAH = load_var_int(dataDirVAH, t, 'theta', nx, ny, nz, xx)
    # taupi
    taupiIdeal = load_var_int(dataDirIdeal, t, 'taupi', nx, ny, nz, xx)
    taupiVH = load_var_int(dataDirVH, t, 'taupi', nx, ny, nz, xx)
    taupiAH = load_var_int(dataDirAH, t, 'taupi', nx, ny, nz, xx)
    taupiVAH = load_var_int(dataDirVAH, t, 'taupi', nx, ny, nz, xx)

    # pimnunu VH
    pittVH = load_var_int(dataDirVH, t, 'pitt', nx, ny, nz, xx)
    pitxVH = load_var_int(dataDirVH, t, 'pitx', nx, ny, nz, xx)
    pityVH = load_var_int(dataDirVH, t, 'pity', nx, ny, nz, xx)
    pixxVH = load_var_int(dataDirVH, t, 'pixx', nx, ny, nz, xx)
    pixyVH = load_var_int(dataDirVH, t, 'pixy', nx, ny, nz, xx)
    piyyVH = load_var_int(dataDirVH, t, 'piyy', nx, ny, nz, xx)
    # pimnunu VAH
    pittVAH = load_var_int(dataDirVAH, t, 'pitt', nx, ny, nz, xx)
    pitxVAH = load_var_int(dataDirVAH, t, 'pitx', nx, ny, nz, xx)
    pityVAH = load_var_int(dataDirVAH, t, 'pity', nx, ny, nz, xx)
    pixxVAH = load_var_int(dataDirVAH, t, 'pixx', nx, ny, nz, xx)
    pixyVAH = load_var_int(dataDirVAH, t, 'pixy', nx, ny, nz, xx)
    piyyVAH = load_var_int(dataDirVAH, t, 'piyy', nx, ny, nz, xx)

    # ptHat
    ptHatAH = (eAH-plAH-Rbar(eAH,pAH,plAH)*(eAH-3*plAH))/2
    ptHatVAH = (eVAH-plVAH-Rbar(eVAH,pVAH,plVAH)*(eVAH-3*plVAH))/2
    # pratio
    pratioVH = np.divide(pVH+pinnVH,pVH-pinnVH/2)
    pratioAH = np.divide(plAH,ptHatAH)
    pratioVAH = np.divide(plVAH,ptHatVAH)

    # enthalpy normalized pimunu
    # VH
    pittEnthalpyNormVH = np.divide(pittVH,eVH+pVH)
    pitxEnthalpyNormVH = np.divide(pitxVH,eVH+pVH)
    pityEnthalpyNormVH = np.divide(pityVH,eVH+pVH)
    pixxEnthalpyNormVH = np.divide(pixxVH,eVH+pVH)
    pixyEnthalpyNormVH = np.divide(pixyVH,eVH+pVH)
    piyyEnthalpyNormVH = np.divide(piyyVH,eVH+pVH)
    pinnEnthalpyNormVH = np.divide(pinnVH,eVH+pVH)
    # VAH
    pittEnthalpyNormVAH = np.divide(pittVAH,eVAH+pVAH)
    pitxEnthalpyNormVAH = np.divide(pitxVAH,eVAH+pVAH)
    pityEnthalpyNormVAH = np.divide(pityVAH,eVAH+pVAH)
    pixxEnthalpyNormVAH = np.divide(pixxVAH,eVAH+pVAH)
    pixyEnthalpyNormVAH = np.divide(pixyVAH,eVAH+pVAH)
    piyyEnthalpyNormVAH = np.divide(piyyVAH,eVAH+pVAH)
    pinnEnthalpyNormVAH = np.divide(pinnVAH,eVAH+pVAH)

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

    print(plt.style.available)

    minorLocator = MultipleLocator(1)

    pad=0.5
    h_pad = None
    w_pad = None
    rect = [0, 0, 1, 1]

    fig, ax = plt.subplots()
    ax.plot(xx, eIdeal*hbarc, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(xx, eVH*hbarc, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, eAH*hbarc, color='blue', linewidth=3.5, linestyle='--', label='AH')
    ax.plot(xx, eVAH*hbarc, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathcal{E}\,\mathrm{[GeV/fm^3]}$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/ed_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(xx, uxIdeal, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(xx, uxVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, uxAH, color='blue', linewidth=3.5, linestyle='--', label='AH')
    ax.plot(xx, uxVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    plt.legend(loc='best', frameon=False)
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$u^x$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/ux_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(xx, np.ones((nxx, 1)), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(xx, pratioVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, pratioAH, color='blue', linewidth=3.5, linestyle='--', label='AH')
    ax.plot(xx, pratioVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathcal{P}_L/\mathcal{P}_\perp$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/pl_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(x, np.ones((nx, 1)), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(x, regVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(x, regVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathrm{tanh}\rho/\rho$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/reg_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    ####################################
    ### plots for Kn ###################
    ####################################

    # T
    fig, ax = plt.subplots()
#    ax.plot(xx, hbarc/taupiIdeal, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(xx, hbarc/taupiVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, hbarc/taupiAH, color='blue', linewidth=3.5, linestyle='--', label='AH')
    ax.plot(xx, hbarc/taupiVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    plt.legend(loc='best', frameon=False)
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel('T [Gev]')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/T_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    # taupi
    fig, ax = plt.subplots()
#    ax.plot(xx, taupiIdeal, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(xx, taupiVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, taupiAH, color='blue', linewidth=3.5, linestyle='--', label='AH')
    ax.plot(xx, taupiVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    plt.legend(loc='best', frameon=False)
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\tau_\pi\,[\mathrm{fm}]$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/taupi_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    # theta
    fig, ax = plt.subplots()
#    ax.plot(xx, thetaIdeal, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(xx, thetaVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, thetaAH, color='blue', linewidth=3.5, linestyle='--', label='AH')
    ax.plot(xx, thetaVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    plt.legend(loc='best', frameon=False)
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\theta\,[\mathrm{fm}^{-1}]$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/theta_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    # KnTaupi
    fig, ax = plt.subplots()
    ax.plot(xx, np.ones((nxx,1)), color='black', linewidth=3.5, linestyle='-')
    ax.plot(xx, np.multiply(taupiVH,np.abs(thetaVH)), color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(xx, np.multiply(taupiAH,np.abs(thetaAH)), color='blue', linewidth=3.5, linestyle='--', label='AH')
    ax.plot(xx, np.multiply(taupiVAH,np.abs(thetaVAH)), color='red', linewidth=3.5, linestyle='--', label='VAH')
    plt.legend(loc='best', frameon=False)
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathrm{Kn}\equiv\theta\tau_\pi$')
    plt.tight_layout(pad=0.5, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/Kn_t-' + '{:.1f}'.format(t) + '.pdf', pad=0.5, h_pad=h_pad, w_pad=w_pad, rect=rect)

    ####################################
    ### plots for \pi^{\mu\nu} #########
    ####################################

    fig, ax = plt.subplots()
    ax.plot(xx, pittVH, color='black', linewidth=3.5, linestyle='-', label=r'$\pi^{\tau\tau}$')
    ax.plot(xx, pitxVH, color='green', linewidth=3.5, linestyle='--', label=r'$\pi^{\tau x}$')
    ax.plot(xx, pityVH, color='blue', linewidth=3.5, linestyle='--', label=r'$\pi^{\tau y}$')
    ax.plot(xx, pixxVH, color='red', linewidth=3.5, linestyle='--', label=r'$\pi^{xx}$')
    first_legend = plt.legend(loc='lower left', frameon=False)
    plt.gca().add_artist(first_legend)
    line5, = ax.plot(xx, pixyVH, color='cyan', linewidth=3.5, linestyle='--', label=r'$\pi^{xy}$')
    line6, = ax.plot(xx, piyyVH, color='orange', linewidth=3.5, linestyle='--', label=r'$\pi^{yy}$')
    line7, = ax.plot(xx, pinnVH, color='purple', linewidth=3.5, linestyle='--', label=r'$\tau^2\pi^{nn}$')
    plt.legend(handles=[line5,line6,line7], loc='lower right', frameon=False)
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'${}^{}$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/pimunuVH_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(xx, pittVAH, color='black', linewidth=3.5, linestyle='-', label=r'$\tilde{\pi}_{\perp}^{\tau\tau}$')
    ax.plot(xx, pitxVAH, color='green', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{\tau x}$')
    ax.plot(xx, pityVAH, color='blue', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{\tau y}$')
    ax.plot(xx, pixxVAH, color='red', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{xx}$')
    first_legend = plt.legend(loc='lower left', frameon=False)
    plt.gca().add_artist(first_legend)
    line5, = ax.plot(xx, pixyVAH, color='cyan', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{xy}$')
    line6, = ax.plot(xx, piyyVAH, color='orange', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{yy}$')
    line7, = ax.plot(xx, pinnVAH, color='purple', linewidth=3.5, linestyle='--', label=r'$\tau^2\tilde{\pi}_{\perp}^{nn}$')
    plt.legend(handles=[line5,line6,line7], loc='lower right', frameon=False)
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'${}^{}$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/pimunuVAH_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    # enthalpyn norm

    fig1, ax1 = plt.subplots()
    ax1.plot(xx, pittEnthalpyNormVH, color='black', linewidth=3.5, linestyle='-', label=r'$\pi^{\tau\tau}$')
    ax1.plot(xx, pitxEnthalpyNormVH, color='green', linewidth=3.5, linestyle='--', label=r'$\pi^{\tau x}$')
    ax1.plot(xx, pityEnthalpyNormVH, color='blue', linewidth=3.5, linestyle='--', label=r'$\pi^{\tau y}$')
    ax1.plot(xx, pixxEnthalpyNormVH, color='red', linewidth=3.5, linestyle='--', label=r'$\pi^{xx}$')
    first_legend = plt.legend(loc='lower left', frameon=False)
    plt.gca().add_artist(first_legend)
    line1_5, = ax1.plot(xx, pixyEnthalpyNormVH, color='cyan', linewidth=3.5, linestyle='--', label=r'$\pi^{xy}$')
    line1_6, = ax1.plot(xx, piyyEnthalpyNormVH, color='orange', linewidth=3.5, linestyle='--', label=r'$\pi^{yy}$')
    line1_7, = ax1.plot(xx, pinnEnthalpyNormVH, color='purple', linewidth=3.5, linestyle='--', label=r'$\tau^2\pi^{nn}$')
    plt.legend(handles=[line1_5,line1_6,line1_7], loc='lower right', frameon=False)
    ax1.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax1.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax1.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\pi^{\mu\nu}/(\mathcal{E}+\mathcal{P}_\mathrm{0})$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig2, ax2 = plt.subplots()
    ax2.plot(xx, pittEnthalpyNormVAH, color='black', linewidth=3.5, linestyle='-', label=r'$\tilde{\pi}_{\perp}^{\tau\tau}$')
    ax2.plot(xx, pitxEnthalpyNormVAH, color='green', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{\tau x}$')
    ax2.plot(xx, pityEnthalpyNormVAH, color='blue', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{\tau y}$')
    ax2.plot(xx, pixxEnthalpyNormVAH, color='red', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{xx}$')
    first_legend = plt.legend(loc='lower left', frameon=False)
    plt.gca().add_artist(first_legend)
    line2_5, = ax2.plot(xx, pixyEnthalpyNormVAH, color='cyan', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{xy}$')
    line2_6, = ax2.plot(xx, piyyEnthalpyNormVAH, color='orange', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{yy}$')
    line2_7, = ax2.plot(xx, pinnEnthalpyNormVAH, color='purple', linewidth=3.5, linestyle='--', label=r'$\tau^2\tilde{\pi}_{\perp}^{nn}$')
    plt.legend(handles=[line2_5,line2_6,line2_7], loc='lower right', frameon=False)
    ax2.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax2.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax2.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\tilde{\pi}_{\perp}^{\mu\nu}/(\mathcal{E}+\mathcal{P}_\mathrm{0})$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    y1 = ax1.get_ylim()
    y2 = ax2.get_ylim()
    minVal = np.min([y1[0],y2[0]])
    maxVal = np.max([y1[1],y2[1]])
    ax1.set_ylim([minVal,maxVal])
    ax2.set_ylim([minVal,maxVal])

    fig1.savefig(plotDir+'/pimunuEnthalpyNormVH_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    fig2.savefig(plotDir+'/pimunuEnthalpyNormVAH_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    plt.show()
