#!/usr/bin/env python3
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.ndimage

def load_var(dir, t, var, nx, ny, nz):
    x, y, n, dataRaw = np.loadtxt(dir + '/' + var + '_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    data = np.reshape(dataRaw, (nx, ny, nz))
    return squeeze(data)

def make_plots(dir, ifig, var, varstr, titleTxt, cmap, t, x, e):
    labelSize = 16
    fig = plt.figure(ifig)
    im = plt.imshow(var, vmin=0, vmax=1, extent=[-xmax, xmax, t0, tf], cmap=cmap, origin='lower')
    e = scipy.ndimage.zoom(e, 3)
    x = scipy.ndimage.zoom(x, 3)
    t = scipy.ndimage.zoom(t, 3)
    plt.contour(x,t,e,efList,colors='k',linewidths=1.5)
    plt.ylim(ymax = tf, ymin = t0)
    plt.tick_params(top='off', right='off')
    plt.title(titleTxt, fontsize=22, y=1.01)
    plt.xlabel(r'$x\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.ylabel(r'$\tau\,[\mathrm{fm/c}]$', fontsize=labelSize)
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(im, cax=cax)
    plt.tight_layout()
    savefig(dir + '/' + varstr + '.pdf', bbox_inches='tight')

def make_plots_e(dir, ifig, var, varstr, titleTxt, cmap, t, x, e):
    labelSize = 16
    fig = plt.figure(ifig)
    im = plt.imshow(var, extent=[-xmax, xmax, t0, tf], cmap=cmap, origin='lower')
    e = scipy.ndimage.zoom(e, 3)
    x = scipy.ndimage.zoom(x, 3)
    t = scipy.ndimage.zoom(t, 3)
    plt.contour(x,t,e,efList,colors='k',linewidths=1.5)
    plt.ylim(ymax = tf, ymin = t0)
    plt.tick_params(top='off', right='off')
    plt.title(titleTxt, fontsize=22, y=1.01)
    plt.xlabel(r'$x\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.ylabel(r'$\tau\,[\mathrm{fm/c}]$', fontsize=labelSize)
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(im, cax=cax)
    plt.tight_layout()
    savefig(dir + '/' + varstr + '.pdf', bbox_inches='tight')

def make_plots_reg(dir, ifig, var, varstr, titleTxt, cmap, t, x, e):
    labelSize = 16
    fig = plt.figure(ifig)
    im = plt.imshow(var, vmin=np.min(np. min(var)), vmax=np.max(np. max(var)), extent=[-xmax, xmax, t0, tf], cmap=cmap, origin='lower')
    e = scipy.ndimage.zoom(e, 3)
    x = scipy.ndimage.zoom(x, 3)
    t = scipy.ndimage.zoom(t, 3)
    plt.contour(x,t,e,efList,colors='w',linewidths=1.5)
    plt.ylim(ymax = tf, ymin = t0)
    plt.tick_params(top='off', right='off')
    plt.title(titleTxt, fontsize=22, y=1.01)
    plt.xlabel(r'$x\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.ylabel(r'$\tau\,[\mathrm{fm/c}]$', fontsize=labelSize)
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(im, cax=cax)
    plt.tight_layout()
    savefig(dir + '/' + varstr + '.pdf', bbox_inches='tight')

if __name__ == '__main__':
    root = '/media/bazow/Data/fluid_dynamic_output_for_thesis/cpu-vah'
#    dataDir = root + '/mcGlb_1d_pimunu_Pi_isotropicIC_QCDEOS_etaOverS0p2_reg0p1-1'
#    plotDir = 'tests/figs/qgp/validity/reg0p1-1'
    dataDir = root + '/mcGlb_1d_pimunu_Pi_isotropicIC_QCDEOS_etaOverS0p2_reg1-10'
    #dataDir = root + '/mcGlb_1d_pimunu_Pi_isotropicIC_QCDEOS_etaOverS0p2_reg1-10_theta1p8'
    plotDir = 'tests/figs/qgp/validity/reg1-10'
    #plotDir = 'tests/figs/qgp/validity/reg1-10_theta1p8'

    nx = 201
    ny = 1
    nz = 1
    dx = 0.1
    dy = 0.1
    dz = 0.1
    xmax = (nx - 1) / 2 * dx
    x = np.linspace(-xmax, xmax, num=nx)

    t0 = 0.5
    tf = 20.5
    tf=20.5
    dt = 0.05
    nt = int((tf-t0)/dt)+1
    t = np.linspace(t0, tf, num=nt)

    efList = [1.8367,9.478]

    reg= np.zeros((nt, nx))
    RPi= np.zeros((nt, nx))
    Rpi= np.zeros((nt, nx))
    RPi2= np.zeros((nt, nx))
    Rpi2= np.zeros((nt, nx))
    Rw= np.zeros((nt, nx))
    fTSol= np.zeros((nt, nx))
    KnTauPi= np.zeros((nt, nx))
    KnTaupi= np.zeros((nt, nx))
    e = np.zeros((nt, nx))
    pl = np.zeros((nt, nx))

    for i in range(0, nt):
        ti = t[i]

        reg[i, :] = load_var(dataDir, ti, 'regulations', nx, ny, nz)
        #fTSol[i, :] = load_var(dataDir, ti, 'fTSol_2', nx, ny, nz)
        e[i, :] = load_var(dataDir, ti, 'e', nx, ny, nz)

    cmap = 'jet'
    cmap2 = 'jet_r'

    plt.rcParams['image.interpolation'] = 'none'

    #####################################################################################################
    # Plots
    #####################################################################################################
#    make_plots_reg(plotDir, 2, reg, 'regulations', r'$\mathrm{tanh}\rho/\rho$', 'Greys', t, x, e)
    make_plots_reg(plotDir, 3, reg, 'regulations', r'$\mathrm{tanh}\rho/\rho$', 'Spectral', t, x, e)
#    make_plots_reg(plotDir, 4, reg, 'regulations', r'$\mathrm{tanh}\rho/\rho$', 'magma', t, x, e)
#    make_plots_reg(plotDir, 5, reg, 'regulations', r'$\mathrm{tanh}\rho/\rho$', 'BrBG_r', t, x, e)

    plt.show()
