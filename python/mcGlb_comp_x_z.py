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

    xp = 60 # x=0
#    xp = 90 # x=3

    t = 3.0
    t2 = t*t

    hbarc = 0.197327

    ############################################
    ### t=1 ####################################
    ############################################
    # For x=0
    zF01 = -4.51637
    zF02 = 4.5164
    # For x=3
    zF01 = -4.42572
    zF02 = 4.42574

    ############################################
    ### t=3 ####################################
    ############################################
    #'''
    # For x=0
    zF01 = -4.69657
    zF02 = 4.69532
    # For x=3
    #zF01 = -4.51751
    #zF02 = 4.51661
    #'''
    e = np.zeros((1, nx))
    pl = np.zeros((1, nx))

    # ed
    eIdeal = load_var_int(dataDirIdeal, t, 'e', nx, ny, nz, xp)
    eVH = load_var_int(dataDirVH, t, 'e', nx, ny, nz, xp)
    eVAH = load_var_int(dataDirVAH, t, 'e', nx, ny, nz, xp)
    # un
    unIdeal = load_var_int(dataDirIdeal, t, 'un', nx, ny, nz, xp)
    unVH = load_var_int(dataDirVH, t, 'un', nx, ny, nz, xp)
    unVAH = load_var_int(dataDirVAH, t, 'un', nx, ny, nz, xp)
    # u
    utVH = load_var_int(dataDirVH, t, 'ut', nx, ny, nz, xp)
    uxVH = load_var_int(dataDirVH, t, 'ux', nx, ny, nz, xp)
    uyVH = load_var_int(dataDirVH, t, 'uy', nx, ny, nz, xp)
    utVAH = load_var_int(dataDirVAH, t, 'ut', nx, ny, nz, xp)
    uxVAH = load_var_int(dataDirVAH, t, 'ux', nx, ny, nz, xp)
    uyVAH = load_var_int(dataDirVAH, t, 'uy', nx, ny, nz, xp)
    # p
    pIdeal = load_var_int(dataDirIdeal, t, 'p', nx, ny, nz, xp)
    pVH = load_var_int(dataDirVH, t, 'p', nx, ny, nz, xp)
    pVAH = load_var_int(dataDirVAH, t, 'p', nx, ny, nz, xp)
    # pl
    plVAH = load_var_int(dataDirVAH, t, 'pl', nx, ny, nz, xp)
    # pinn
    pinnVH = load_var_int(dataDirVH, t, 'pinn', nx, ny, nz, xp)*t2
    pinnVAH = load_var_int(dataDirVAH, t, 'pinn', nx, ny, nz, xp)*t2
    # reg
    regVH = load_var_int(dataDirVH, t, 'regulations', nx, ny, nz, xp)
    regVAH = load_var_int(dataDirVAH, t, 'regulations', nx, ny, nz, xp)

    # theta
    thetaIdeal = load_var_int(dataDirIdeal, t, 'theta', nx, ny, nz, xp)
    thetaVH = load_var_int(dataDirVH, t, 'theta', nx, ny, nz, xp)
    thetaVAH = load_var_int(dataDirVAH, t, 'theta', nx, ny, nz, xp)
    # taupi
    taupiIdeal = load_var_int(dataDirIdeal, t, 'taupi', nx, ny, nz, xp)
    taupiVH = load_var_int(dataDirVH, t, 'taupi', nx, ny, nz, xp)
    taupiVAH = load_var_int(dataDirVAH, t, 'taupi', nx, ny, nz, xp)

    # pimnunu VH
    pittVH = load_var_int(dataDirVH, t, 'pitt', nx, ny, nz, xp)
    pitxVH = load_var_int(dataDirVH, t, 'pitx', nx, ny, nz, xp)
    pityVH = load_var_int(dataDirVH, t, 'pity', nx, ny, nz, xp)
    pitnVH = load_var_int(dataDirVH, t, 'pitn', nx, ny, nz, xp)
    pixxVH = load_var_int(dataDirVH, t, 'pixx', nx, ny, nz, xp)
    pixyVH = load_var_int(dataDirVH, t, 'pixy', nx, ny, nz, xp)
    pixnVH = load_var_int(dataDirVH, t, 'pixn', nx, ny, nz, xp)
    piyyVH = load_var_int(dataDirVH, t, 'piyy', nx, ny, nz, xp)
    piynVH = load_var_int(dataDirVH, t, 'piyn', nx, ny, nz, xp)
    # pimnunu VAH
    pittVAH = load_var_int(dataDirVAH, t, 'pitt', nx, ny, nz, xp)
    pitxVAH = load_var_int(dataDirVAH, t, 'pitx', nx, ny, nz, xp)
    pityVAH = load_var_int(dataDirVAH, t, 'pity', nx, ny, nz, xp)
    pitnVAH = load_var_int(dataDirVAH, t, 'pitn', nx, ny, nz, xp)
    pixxVAH = load_var_int(dataDirVAH, t, 'pixx', nx, ny, nz, xp)
    pixyVAH = load_var_int(dataDirVAH, t, 'pixy', nx, ny, nz, xp)
    pixnVAH = load_var_int(dataDirVAH, t, 'pixn', nx, ny, nz, xp)
    piyyVAH = load_var_int(dataDirVAH, t, 'piyy', nx, ny, nz, xp)
    piynVAH = load_var_int(dataDirVAH, t, 'piyn', nx, ny, nz, xp)

    # ptHat
    ptHatVAH = (eVAH-plVAH-Rbar(eVAH,pVAH,plVAH)*(eVAH-3*pVAH))/2
    # pratio
    pratioVH = np.divide(pVH+pinnVH,pVH-pinnVH/2)
    pratioVAH = np.divide(plVAH,ptHatVAH)

    # W
    WtVAH = load_var_int(dataDirVAH, t, 'WtTz', nx, ny, nz, xp)
    WxVAH = load_var_int(dataDirVAH, t, 'WxTz', nx, ny, nz, xp)
    WyVAH = load_var_int(dataDirVAH, t, 'WyTz', nx, ny, nz, xp)
    WnVAH = load_var_int(dataDirVAH, t, 'WnTz', nx, ny, nz, xp)

    #
    z0 = np.divide(t*unVAH,np.sqrt(1+np.multiply(uxVAH,uxVAH)+np.multiply(uyVAH,uyVAH)))
    z3 = np.divide(utVAH/t,np.sqrt(1+np.multiply(uxVAH,uxVAH)+np.multiply(uyVAH,uyVAH)))
    pitt2 = pittVAH + 2*np.multiply(WtVAH,z0) + np.multiply((3*plVAH-eVAH)/2,2*np.multiply(z0,z0)+1-np.multiply(utVAH,utVAH)+np.multiply(z0,z0))/3
    pitx2 = pitxVAH + np.multiply(WxVAH,z0) + np.multiply((3*plVAH-eVAH)/2,-np.multiply(utVAH,uxVAH))/3
    pity2 = pityVAH + np.multiply(WyVAH,z0) + np.multiply((3*plVAH-eVAH)/2,-np.multiply(utVAH,uyVAH))/3
    pitn2 = pitnVAH + np.multiply(WtVAH,z3) + np.multiply(WnVAH,z0) + np.multiply((3*plVAH-eVAH)/2,2*np.multiply(z0,z3)-np.multiply(utVAH,unVAH)+np.multiply(z0,z3))/3
    pixx2 = pixxVAH + np.multiply((3*plVAH-eVAH)/2,-1-np.multiply(uxVAH,uxVAH))/3
    pixy2 = pixyVAH + np.multiply((3*plVAH-eVAH)/2,-np.multiply(uxVAH,uyVAH))/3
    pixn2 = pixnVAH + np.multiply(WxVAH,z3) + np.multiply((3*plVAH-eVAH)/2,-np.multiply(uxVAH,unVAH))/3
    piyy2 = piyyVAH + np.multiply((3*plVAH-eVAH)/2,-1-np.multiply(uyVAH,uyVAH))/3
    piyn2 = piynVAH + np.multiply(WyVAH,z3) + np.multiply((3*plVAH-eVAH)/2,-np.multiply(uyVAH,unVAH))/3
    pinn2 = pinnVAH + 2*np.multiply(WnVAH,z3) + np.multiply((3*plVAH-eVAH)/2,2*np.multiply(z3,z3)-1/t/t-np.multiply(unVAH,unVAH)+np.multiply(z3,z3))/3

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
    ax.plot(z, eIdeal*hbarc, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(z, eVH*hbarc, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(z, eVAH*hbarc, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax.axvspan(zF02, zmax, alpha=0.25, color=colorFO)
    ax.axvspan(zF01, -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\mathcal{E}\,\mathrm{[GeV/fm^3]}$')
    plt.title(r'$\vec{r}=($' + '{:.1f}'.format(x[xp]) + ', 0.0) fm',fontsize=16)
    #text(0.1, 0.9, '(a)', ha='center', va='center', transform=ax.transAxes,fontsize=16)
    text(0.5, 0.35,
         r'$\tau=$' + '{:.1f}'.format(t) + '\n' + r'$\tau_0=0.5\,\mathrm{fm/c}$' + '\n' r'$T_0=0.6\,\mathrm{GeV}$' + '\n' + r'$\eta/s=0.2$',
         ha='center', va='center', transform=ax.transAxes,
         fontsize=16)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/ed_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xp])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(z, unIdeal, color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(z, unVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(z, unVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    plt.legend(loc='best', frameon=False)
    ax.axvspan(zF02, zmax, alpha=0.25, color=colorFO)
    ax.axvspan(zF01, -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$u^\eta$')
    #text(0.1, 0.9, '(b)', ha='center', va='center', transform=ax.transAxes,fontsize=16)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/un_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xp])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(z, np.ones((nz, 1)), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(z, pratioVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(z, pratioVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax.axvspan(zF02, zmax, alpha=0.25, color=colorFO)
    ax.axvspan(zF01, -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\mathcal{P}_L/\mathcal{P}_\perp$')
    #text(0.08, 0.9, '(c)', ha='center', va='center', transform=ax.transAxes,fontsize=16)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/pratio_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xp])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(z, np.ones((nz, 1)), color='black', linewidth=3.5, linestyle='-', label='Ideal')
    ax.plot(z, regVH, color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(z, regVAH, color='red', linewidth=3.5, linestyle='--', label='VAH')
    ax.axvspan(zF02, zmax, alpha=0.25, color=colorFO)
    ax.axvspan(zF01, -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\mathrm{tanh}\rho/\rho$')
    #text(0.1, 0.1, '(d)', ha='center', va='center', transform=ax.transAxes,fontsize=16)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/reg_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xp])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    ####################################
    ### plots for Kn ###################
    ####################################

    # KnTaupi
    fig, ax = plt.subplots()
    ax.plot(z, np.ones((nz,1)), color='black', linewidth=3.5, linestyle='-')
    ax.plot(z, np.multiply(taupiVH,np.abs(thetaVH)), color='green', linewidth=3.5, linestyle='--', label='VH')
    ax.plot(z, np.multiply(taupiVAH,np.abs(thetaVAH)), color='red', linewidth=3.5, linestyle='--', label='VAH')
    plt.legend(loc='best', frameon=False)
    ax.axvspan(zF02, zmax, alpha=0.25, color=colorFO)
    ax.axvspan(zF01, -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\mathrm{Kn}\equiv\theta\tau_\pi$')
    #text(0.1, 0.9, '(d)', ha='center', va='center', transform=ax.transAxes,fontsize=16)
    plt.tight_layout(pad=0.5, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/Kn_t-' + '{:.1f}'.format(t) + '_x-' + '{:.1f}'.format(x[xp])
            + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    ####################################
    ### plots for pimunu ###############
    ####################################
    colors = np.array(['black', 'green', 'blue', 'red', 'cyan', 'orange', 'purple'])

    fig1, ax1 = plt.subplots()
    ax1.plot(z, pittVH, color=colors[0], linewidth=3.5, linestyle='-', label=r'$\pi^{\tau\tau}$')
    ax1.plot(z, pitxVH, color=colors[1], linewidth=3.5, linestyle='--', label=r'$\pi^{\tau x}$')
    ax1.plot(z, pityVH, color=colors[2], linewidth=3.5, linestyle='--', label=r'$\pi^{\tau y}$')
    ax1.plot(z, pixxVH, color=colors[3], linewidth=3.5, linestyle='--', label=r'$\pi^{xx}$')
    #first_legend = plt.legend(loc='lower left', frameon=False)
    #plt.gca().add_artist(first_legend)
    line1_5, = ax1.plot(z, pixyVH, color=colors[4], linewidth=3.5, linestyle='--', label=r'$\pi^{xy}$')
    line1_6, = ax1.plot(z, piyyVH, color=colors[5], linewidth=3.5, linestyle='--', label=r'$\pi^{yy}$')
    line1_7, = ax1.plot(z, pinnVH, color=colors[6], linewidth=3.5, linestyle='--', label=r'$\tau^2\pi^{nn}$')
    #plt.legend(handles=[line1_5,line1_6,line1_7], loc='lower right', frameon=False)
    ax1.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\pi^{\mu\nu}/(\mathcal{E}+\mathcal{P}_\mathrm{0})$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig2, ax2 = plt.subplots()
    ax2.plot(z, pitt2, color='black', linewidth=3.5, linestyle='', marker='o' ,label=r'$\tilde{\pi}_{\perp}^{\tau\tau}$')
    ax2.plot(z, pitx2, color='green', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{\tau x}$')
    ax2.plot(z, pity2, color='blue', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{\tau y}$')
    ax2.plot(z, pixx2, color='red', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{xx}$')
    #first_legend = plt.legend(loc='lower left', frameon=False)
    #plt.gca().add_artist(first_legend)
    first_legend = plt.legend(loc='upper center', frameon=False)
    plt.gca().add_artist(first_legend)
    line2_5, = ax2.plot(z, pixy2, color='cyan', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{xy}$')
    line2_6, = ax2.plot(z, piyy2, color='orange', linewidth=3.5, linestyle='--', label=r'$\tilde{\pi}_{\perp}^{yy}$')
    line2_7, = ax2.plot(z, pinn2, color='purple', linewidth=3.5, linestyle='--', label=r'$\tau^2\tilde{\pi}_{\perp}^{nn}$')
    #plt.legend(handles=[line2_5,line2_6,line2_7], loc='lower right', frameon=False)
    plt.legend(handles=[line2_5,line2_6,line2_7], loc='lower center', frameon=False)
    ax2.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\tilde{\pi}_{\perp}^{\mu\nu}/(\mathcal{E}+\mathcal{P}_\mathrm{0})$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    y1 = ax1.get_ylim()
    y2 = ax2.get_ylim()
    minVal = np.min([y1[0],y2[0]])
    maxVal = np.max([y1[1],y2[1]])
    ax1.set_ylim([minVal,maxVal])
    ax2.set_ylim([minVal,maxVal])

    ax1.text(0.1, 0.9, '(b)', ha='center', va='center', transform=ax.transAxes,fontsize=16)
    ax2.text(0.1, 0.9, '(d)', ha='center', va='center', transform=ax.transAxes,fontsize=16)

    fig1, ax1 = plt.subplots()
    ms = 2
    ax1.plot(z, pittVH, color=colors[0], markersize=ms, linestyle='', marker='o', label=r'$\pi^{\tau\tau}$')
    ax1.plot(z, pitnVH, color=colors[1], markersize=ms, linestyle='', marker='o', label=r'$\pi^{\tau\eta}$')
    ax1.plot(z, pixxVH, color=colors[2], markersize=ms, linestyle='', marker='o', label=r'$\pi^{xx}$')
    ax1.plot(z, pixnVH, color=colors[3], markersize=ms, linestyle='', marker='o', label=r'$\pi^{x\eta}$')
    #first_legend = plt.legend(loc='lower left', frameon=False)
    #plt.gca().add_artist(first_legend)
    line1_5, = ax1.plot(z, piyyVH, color=colors[4], markersize=ms, linestyle='', marker='o', label=r'$\pi^{yy}$')
    line1_6, = ax1.plot(z, piynVH, color=colors[5], markersize=ms, linestyle='', marker='o', label=r'$\pi^{y\eta}$')
    line1_7, = ax1.plot(z, pinnVH, color=colors[6], markersize=ms, linestyle='', marker='o', label=r'$\tau^2\pi^{nn}$')
    plt.legend(loc='best', frameon=False)
    #plt.legend(handles=[line1_5,line1_6,line1_7], loc='lower right', frameon=False)
    #### VAH ####
    lw = 2.5
    ax1.plot(z, pitt2, color='black', linewidth=lw, linestyle='-',label=r'$\tilde{\pi}_{\perp}^{\tau\tau}$')
    ax1.plot(z, pitn2, color='green', linewidth=lw, linestyle='-', label=r'$\tilde{\pi}_{\perp}^{\tau x}$')
    ax1.plot(z, pixx2, color='blue', linewidth=lw, linestyle='-', label=r'$\tilde{\pi}_{\perp}^{\tau y}$')
    ax1.plot(z, pixn2, color='red', linewidth=lw, linestyle='-', label=r'$\tilde{\pi}_{\perp}^{xx}$')
    ax1.plot(z, piyy2, color='cyan', linewidth=lw, linestyle='-', label=r'$\tilde{\pi}_{\perp}^{xy}$')
    ax1.plot(z, piyn2, color='orange', linewidth=lw, linestyle='-', label=r'$\tilde{\pi}_{\perp}^{yy}$')
    ax1.plot(z, pinn2, color='purple', linewidth=lw, linestyle='-', label=r'$\tau^2\tilde{\pi}_{\perp}^{nn}$')
    ###################
    ax1.axvspan(zF02, zmax, alpha=0.25, color=colorFO)
    ax1.axvspan(zF01, -zmax, alpha=0.25, color=colorFO)
    pylab.xlim([-zmax,zmax])
    ax1.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\eta_s$')
    plt.ylabel(r'$\pi^{\mu\nu}$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    plt.show()
