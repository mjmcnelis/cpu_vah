import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

nx = 201
dx = 0.1
xmax = (nx - 1) / 2 * dx

t0 = 0.5
tf = 10.5


def make_plots(dir, ifig, var, varstr, titleTxt, cmap, xtc, ttc, xtf1, ttf1, xtf2, ttf2, xtf3, ttf3):
    labelSize = 16
    fig = plt.figure(ifig)
    im = plt.imshow(var, vmin=0, vmax=1, extent=[-xmax, xmax, t0, tf], cmap=cmap, origin='lower')
    plt.plot(xtc,ttc,'w-')
    plt.plot(xtf1, ttf1, 'w-')
    plt.plot(xtf2, ttf2, 'w-')
    plt.plot(xtf3, ttf3, 'w-')
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
    outputDir = 'tests/output/qgp/validity'
    plotDir = 'tests/figs/qgp/validity'

    xtc, ttc = np.loadtxt('/home/bazow/Documents/tc_xy_surface.dat', unpack=True)
    xtf1, ttf1 = np.loadtxt('/home/bazow/Documents/tf_surface_1.dat', unpack=True)
    xtf2, ttf2 = np.loadtxt('/home/bazow/Documents/tf_surface_2.dat', unpack=True)
    xtf3, ttf3 = np.loadtxt('/home/bazow/Documents/tf_surface_3.dat', unpack=True)

    reg = np.loadtxt(outputDir + '/reg.dat', unpack=True)
#    R2pi = np.loadtxt(outputDir + '/R2pi.dat', unpack=True)
    Rpi = np.loadtxt(outputDir + '/Rpi.dat', unpack=True)
#    R2Pi = np.loadtxt(outputDir + '/R2Pi.dat', unpack=True)
#    RPi = np.loadtxt(outputDir + '/RPi.dat', unpack=True)
    KnTaupi = np.loadtxt(outputDir + '/KnTaupi.dat', unpack=True)
 #   KnTauPi = np.loadtxt(outputDir + '/KnTauPi.dat', unpack=True)
    fTSol = np.loadtxt(outputDir + '/fTSol.dat', unpack=True)

    cmap = 'jet'
    cmap2 = 'jet_r'

    plt.rcParams['image.interpolation'] = 'none'

    #####################################################################################################
    # Plots
    #####################################################################################################
    make_plots(plotDir, 0, KnTaupi, 'KnTaupi', r'$\mathrm{Kn}_{\theta\pi}\equiv\theta\tau_{\pi}$', cmap, xtc, ttc, xtf1, ttf1, xtf2, ttf2, xtf3, ttf3)
#    make_plots(plotDir, 1, KnTauPi, 'KnTauPi', r'$\mathrm{Kn}_{\theta\Pi}\equiv\theta\tau_{\Pi}$', cmap, xtc, ttc, xtf1, ttf1, xtf2, ttf2, xtf3, ttf3)
#    make_plots(plotDir, 2, RPi, 'RPi', r'$\tilde{\mathrm{R}}^{-1}_{\Pi}$', cmap, xtc, ttc, xtf1, ttf1, xtf2, ttf2, xtf3, ttf3)
    make_plots(plotDir, 3, Rpi, 'Rpi', r'$\tilde{\mathrm{R}}^{-1}_{\pi}$', cmap, xtc, ttc, xtf1, ttf1, xtf2, ttf2, xtf3, ttf3)
    make_plots(plotDir, 4, fTSol, 'fTSol', r'$f_\perp(u_\perp)\neq 0\,\forall u_\perp\geq 0$', cmap, xtc, ttc, xtf1, ttf1, xtf2, ttf2, xtf3, ttf3)
    make_plots(plotDir, 5, reg, 'regulations', r'$\mathrm{tanh}\rho/\rho$', cmap2, xtc, ttc, xtf1, ttf1, xtf2, ttf2, xtf3, ttf3)

    plt.show()
