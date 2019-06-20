#!/usr/bin/env python3
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.ndimage
import matplotlib.pyplot as plt
import pylab
from scipy import interpolate
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import brewer2mpl

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
    return (2.473173363908116e-10 - 4.899839370307281e-6*a + 155055.91462124084*a10 - 275435.45350226434*a11 + 350689.68825705117*a12 -
     299725.38986957155*a13 + 151477.08809203724*a14 - 33196.47417939176*a15 - 0.004301975027942015*a2 - 0.14858206981041563*a3 +
     6.249255189587875*a4 - 92.79641927240235*a5 + 807.175057749925*a6 - 4760.015905266286*a7 + 20324.533122685436*a8 -
     64758.869552496515*a9)/(0.00008222793468208523 + 0.03411917870833943*a + 4.895969276094396e6*a10 - 8.84162305829353e6*a11 +
     1.1445063656613324e7*a12 - 9.918442713390596e6*a13 + 5.065882388219598e6*a14 - 1.1181016364928822e6*a15 + 0.23871740573818725*a2 -
     23.50912574236691*a3 + 417.4953123877312*a4 - 4234.215775452717*a5 + 29824.022790048104*a6 - 157419.8447785501*a7 +
     641300.6529027821*a8 - 2.0248032895288002e6*a9)

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
    dataDirIdeal = root + 'gpu-vh/Glb_1d_ideal_isotropicIC_conformalEOS'
    # VH
    dataDirVH = root + 'gpu-vh/Glb_1d_shear_isotropicIC_conformalEOS_etaOverS0p2_noreg'
    dataDirVHReg = root + '/gpu-vh/Glb_1d_shear_isotropicIC_conformalEOS_etaOverS0p2_reg0p1-1'
    dataDirVHReg10 = root + 'gpu-vh/Glb_1d_shear_isotropicIC_conformalEOS_etaOverS0p2_reg1-10'
    # AH
    dataDirAH = root + 'cpu-vah/Glb_1d_shear_LOAhydro_isotropicIC_conformalEOS_etaOverS0p2'
    # VAH
    dataDirVAH = root + 'cpu-vah/Glb_1d_shear_withPimunu_isotropicIC_conformalEOS_etaOverS0p2_noreg'
    dataDirVAHReg = root + 'cpu-vah/Glb_1d_shear_withPimunu_isotropicIC_conformalEOS_etaOverS0p2_reg0p1-1'
    dataDirVAHReg10 = root + 'cpu-vah/Glb_1d_shear_withPimunu_WTz_isotropicIC_conformalEOS_etaOverS0p2_reg1-10'

    plotDir = 'tests/figs/qgp/Glb-reg-comp'

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

    # For t=1 -- mcGlb
#    xF01 = -5.57862
#    xF02 = 5.4526
    # For t=3 -- mcGlb
#    xF01 = -6.34775
#    xF02 = 6.04202

    e = np.zeros((1, nx))
    pl = np.zeros((1, nx))

    # ed
    eIdeal = load_var_int(dataDirIdeal, t, 'e', nx, ny, nz, xx)
    eVH = load_var_int(dataDirVH, t, 'e', nx, ny, nz, xx)
    eVHReg = load_var_int(dataDirVHReg, t, 'e', nx, ny, nz, xx)
    eVHReg10 = load_var_int(dataDirVHReg10, t, 'e', nx, ny, nz, xx)
    eAH = load_var_int(dataDirAH, t, 'e', nx, ny, nz, xx)
    eVAH = load_var_int(dataDirVAH, t, 'e', nx, ny, nz, xx)
    eVAHReg = load_var_int(dataDirVAHReg, t, 'e', nx, ny, nz, xx)
    eVAHReg10 = load_var_int(dataDirVAHReg10, t, 'e', nx, ny, nz, xx)
    # ux
    uxIdeal = load_var_int(dataDirIdeal, t, 'ux', nx, ny, nz, xx)
    uxVH = load_var_int(dataDirVH, t, 'ux', nx, ny, nz, xx)
    uxVHReg = load_var_int(dataDirVHReg, t, 'ux', nx, ny, nz, xx)
    uxVHReg10 = load_var_int(dataDirVHReg10, t, 'ux', nx, ny, nz, xx)
    uxAH = load_var_int(dataDirAH, t, 'ux', nx, ny, nz, xx)
    uxVAH = load_var_int(dataDirVAH, t, 'ux', nx, ny, nz, xx)
    uxVAHReg = load_var_int(dataDirVAHReg, t, 'ux', nx, ny, nz, xx)
    uxVAHReg10 = load_var_int(dataDirVAHReg10, t, 'ux', nx, ny, nz, xx)
    # p
    pIdeal = load_var_int(dataDirIdeal, t, 'p', nx, ny, nz, xx)
    pVH = load_var_int(dataDirVH, t, 'p', nx, ny, nz, xx)
    pVHReg = load_var_int(dataDirVHReg, t, 'p', nx, ny, nz, xx)
    pVHReg10 = load_var_int(dataDirVHReg10, t, 'p', nx, ny, nz, xx)
    pAH = load_var_int(dataDirAH, t, 'p', nx, ny, nz, xx)
    pVAH = load_var_int(dataDirVAH, t, 'p', nx, ny, nz, xx)
    pVAHReg = load_var_int(dataDirVAHReg, t, 'p', nx, ny, nz, xx)
    pVAHReg10 = load_var_int(dataDirVAHReg10, t, 'p', nx, ny, nz, xx)
    # pl
    plAH = load_var_int(dataDirAH, t, 'pl', nx, ny, nz, xx)
    plVAH = load_var_int(dataDirVAH, t, 'pl', nx, ny, nz, xx)
    plVAHReg = load_var_int(dataDirVAHReg, t, 'pl', nx, ny, nz, xx)
    plVAHReg10 = load_var_int(dataDirVAHReg10, t, 'pl', nx, ny, nz, xx)
    # pinn
    pinnVH = load_var_int(dataDirVH, t, 'pinn', nx, ny, nz, xx)*t2
    pinnVHReg = load_var_int(dataDirVHReg, t, 'pinn', nx, ny, nz, xx)*t2
    pinnVHReg10 = load_var_int(dataDirVHReg10, t, 'pinn', nx, ny, nz, xx)*t2
#    pinnVAH = load_var_int(dataDirVAH, t, 'pinn', nx, ny, nz, xx)*t2
#    pinnVAHReg = load_var_int(dataDirVAHReg, t, 'pinn', nx, ny, nz, xx)*t2
#    pinnVAHReg10 = load_var_int(dataDirVAHReg10, t, 'pinn', nx, ny, nz, xx)*t2
    # reg
    regVH = load_var(dataDirVH, t, 'regulations', nx, ny, nz)
    regVHReg = load_var(dataDirVHReg, t, 'regulations', nx, ny, nz)
    regVHReg10 = load_var(dataDirVHReg10, t, 'regulations', nx, ny, nz)
    regVAH = load_var(dataDirVAH, t, 'regulations', nx, ny, nz)
    regVAHReg = load_var(dataDirVAHReg, t, 'regulations', nx, ny, nz)
    regVAHReg10 = load_var(dataDirVAHReg10, t, 'regulations', nx, ny, nz)
    '''
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
    '''

    # ptHat
    ptHatAH = (eAH-plAH-Rbar(eAH,pAH,plAH)*(eAH-3*plAH))/2
    ptHatVAH = (eVAH-plVAH-Rbar(eVAH,pVAH,plVAH)*(eVAH-3*plVAH))/2
    ptHatVAHReg = (eVAHReg-plVAHReg-Rbar(eVAHReg,pVAHReg,plVAHReg)*(eVAHReg-3*plVAHReg))/2
    ptHatVAHReg10 = (eVAHReg10-plVAHReg10-Rbar(eVAHReg10,pVAHReg10,plVAHReg10)*(eVAHReg10-3*plVAHReg10))/2
    # pratio
    pratioVH = np.divide(pVH+pinnVH,pVH-pinnVH/2)
    pratioVHReg = np.divide(pVHReg+pinnVHReg,pVHReg-pinnVHReg/2)
    pratioVHReg10 = np.divide(pVHReg10+pinnVHReg10,pVHReg10-pinnVHReg10/2)
    pratioAH = np.divide(plAH,ptHatAH)
    pratioVAH = np.divide(plVAH,ptHatVAH)
    pratioVAHReg = np.divide(plVAHReg,ptHatVAHReg)
    pratioVAHReg10 = np.divide(plVAHReg10,ptHatVAHReg10)

    #####################################################################################################
    # Plots
    #####################################################################################################
    plt.style.use('fivethirtyeight')
    plt.style.use('seaborn-colorblind')
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
    ax.plot(xx, uxVH, color='red', linewidth=3.5, linestyle='--', label='VH-no reg.')
    ax.plot(xx, uxVHReg, color='blue', linewidth=3.5, linestyle='--', label='VH-(0.1,1)')
    ax.plot(xx, uxVHReg10, color='green', linewidth=3.5, linestyle='--', label='VH-(1,10)')
    ax.plot(xx, uxAH, color='orange', linewidth=3.5, linestyle='--', label='AH')
    ax.plot(xx, uxVAH, color='salmon', linewidth=3.5, linestyle='--', label='VAH-no reg.')
    ax.plot(xx, uxVAHReg, color='cyan', linewidth=3.5, linestyle='--', label='VAH-(0.1,1)')
    ax.plot(xx, uxVAHReg10, color='purple', linewidth=3.5, linestyle='--', label='VAH-(1,10)')
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
    ax.plot(xx, pratioVH, color='red', linewidth=3.5, linestyle='--', label='VH-no reg.')
    ax.plot(xx, pratioVHReg, color='blue', linewidth=3.5, linestyle='--', label='VH-(0.1,1)')
    ax.plot(xx, pratioVHReg10, color='green', linewidth=3.5, linestyle='--', label='VH-(1,10)')
    ax.plot(xx, pratioAH, color='orange', linewidth=3.5, linestyle='--', label='AH')
    ax.plot(xx, pratioVAH, color='salmon', linewidth=3.5, linestyle='--', label='VAH-no reg.')
    ax.plot(xx, pratioVAHReg, color='cyan', linewidth=3.5, linestyle='--', label='VAH-(0.1,1)')
    ax.plot(xx, pratioVAHReg10, color='purple', linewidth=3.5, linestyle='--', label='VAH-(1,10)')
    ax.axvspan(xF02, xmax, alpha=0.25, color='dodgerblue')
    ax.axvspan(xF01, -xmax, alpha=0.25, color='dodgerblue')
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathcal{P}_L/\mathcal{P}_\perp$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/pl_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    #### TEMP
    # Valid names are: ['Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3']
    bmap = brewer2mpl.get_map('Set2', 'qualitative', 5)
    colors = bmap.mpl_colors

    fig, ax = plt.subplots()
    ax.plot(xx, np.ones((nxx, 1)), color=colors[3], linewidth=2, linestyle='-', label='Ideal')
    ax.fill_between(xx, pratioVH, pratioVHReg, alpha=0.65, linewidth=0, color=colors[0])
    ax.plot(xx, pratioVHReg10, color=colors[0], linewidth=2.0, linestyle='-', label='VH-(1,10)')
    ax.fill_between(xx, pratioVAH, pratioVAHReg, alpha=0.65, linewidth=0, color=colors[1])
    ax.plot(xx, pratioAH, color=colors[3], linewidth=2.0, linestyle='--', label='AH')
    ax.axvspan(xF02, xmax, alpha=0.25, color=colors[2])
    ax.axvspan(xF01, -xmax, alpha=0.25, color=colors[2])
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathcal{P}_L/\mathcal{P}_\perp$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/pl_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    ax.plot(x, np.ones((nx, 1)), color='black', linewidth=2.0, linestyle='-', label='Ideal')
    ax.fill_between(x, regVHReg, regVHReg10, alpha=0.25, linewidth=0, color=colors[0])
    ax.fill_between(x, regVAHReg, regVAHReg10, alpha=0.25, linewidth=0, color=colors[1])
    ax.axvspan(xF02, xmax, alpha=0.25, color=colors[2])
    ax.axvspan(xF01, -xmax, alpha=0.25, color=colors[2])
    pylab.xlim([-xmax,xmax])
    xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$x\,[\mathrm{fm}]$')
    plt.ylabel(r'$\mathrm{tanh}\rho/\rho$')
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/reg_regComp_t-' + '{:.1f}'.format(t) + '.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    '''
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
    '''
    ####################################
    ### plots for Kn ###################
    ####################################
    '''
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
    '''
    plt.show()
