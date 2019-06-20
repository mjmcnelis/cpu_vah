#!/usr/bin/env python3
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

def load_var(dir, t, var, nx, ny, nz):
    x, y, n, dataRaw = np.loadtxt(dir + '/' + var + '_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    data = np.reshape(dataRaw, (nx, ny, nz))
    return squeeze(data)

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def make_plots(dir, ifig, var, varstr, titleTxt, cmap, vmin, vmax):
    labelSize = 16
    fig = plt.figure(ifig)
    #ax = fig.add_subplot(111)
    im = plt.imshow(var, vmin=vmin, vmax=vmax, extent=[-xmax, xmax, t0, tf], cmap=cmap, origin='lower')
    plt.ylim(ymax = tf, ymin = t0)
    plt.tick_params(top='off', right='off')
    plt.title(titleTxt, fontsize=22, y=1.01)
    plt.xlabel(r'$x\,[\mathrm{fm}]$', fontsize=labelSize)
    plt.ylabel(r'$\tau\,[\mathrm{fm/c}]$', fontsize=labelSize)
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(im, cax=cax)
    #############################################
    #forceAspect(ax,aspect=1)
    ############################################
    plt.tight_layout()
    savefig(dir + '/' + varstr + '.pdf', bbox_inches='tight')

if __name__ == '__main__':
    root = '/media/bazow/Data/fluid_dynamic_output_for_thesis/'
    type_vh = 'gpu-vh/'
    type_vah = 'cpu-vah/'
    tag = 'Glb_1d_shear_isotropicIC_conformalEOS_etaOverS0p2_reg1-10'
    dataDir_vh = root + type_vh + tag
    tag = 'Glb_1d_shear_withPimunu_WTz_isotropicIC_conformalEOS_etaOverS0p2_reg1-10'
    dataDir_vah = root + type_vah + tag
    plotDir_vh = 'tests/figs/qgp/Glb/vh'
    plotDir_vah = 'tests/figs/qgp/Glb/vah'

    nx = 201
    ny = 1
    nz = 1
    dx = 0.1
    dy = 0.1
    dz = 0.1
    xmax = (nx - 1) / 2 * dx
    x = np.linspace(-xmax, xmax, num=nx)

    t0 = 0.5
    tf = 10.0
    dt = 0.05
    nt = int((tf-t0)/dt)+1
    t = np.linspace(t0, tf, num=nt)

    hbarc = 0.197327

    # Viscous hydro
    KnTaupi_VH= np.zeros((nt, nx))
    e_VH= np.zeros((nt, nx))
    p_VH= np.zeros((nt, nx))
    t2pinn_VH= np.zeros((nt, nx))
    reg_VH= np.zeros((nt, nx))
    # Anisotropic hydro
    KnTaupi_VAH= np.zeros((nt, nx))
    e_VAH= np.zeros((nt, nx))
    pl_VAH= np.zeros((nt, nx))
    fTSol_VAH= np.zeros((nt, nx))
    reg_VAH= np.zeros((nt, nx))

    for i in range(0, nt):
        ti = t[i]

        # Viscous hydro
        KnTaupi_VH[i, :] = load_var(dataDir_vh, ti, 'KnTaupi', nx, ny, nz)
        e_VH[i, :] = load_var(dataDir_vh, ti, 'e', nx, ny, nz)
        p_VH[i, :] = load_var(dataDir_vh, ti, 'p', nx, ny, nz)
        t2pinn_VH[i, :] = ti*ti*load_var(dataDir_vh, ti, 'pinn', nx, ny, nz)
        reg_VH[i, :] = load_var(dataDir_vh, ti, 'regulations', nx, ny, nz)
        # Anisotropic hydro
        KnTaupi_VAH[i, :] = load_var(dataDir_vah, ti, 'knTaupi', nx, ny, nz)
        e_VAH[i, :] = load_var(dataDir_vah, ti, 'e', nx, ny, nz)
        pl_VAH[i, :] = load_var(dataDir_vah, ti, 'pl', nx, ny, nz)
        fTSol_VAH[i, :] = load_var(dataDir_vah, ti, 'fTSol_2', nx, ny, nz)
        reg_VAH[i, :] = load_var(dataDir_vah, ti, 'regulations', nx, ny, nz)

    vmin_e = np.min(np.min(np.min(e_VH)))
    vmax_e = np.max(np.max(np.max(e_VH)))

    cmap = 'jet'
    cmap2 = 'jet_r'

    plt.rcParams['image.interpolation'] = 'none'

    #####################################################################################################
    # Plots
    #####################################################################################################
    # Viscous hydro
    make_plots(plotDir_vh, 0, KnTaupi_VH, 'KnTaupi', r'$\mathrm{Kn}_{\theta\pi}\equiv\theta\tau_{\pi}$', cmap, 0, 1)
    make_plots(plotDir_vh, 1, e_VH, 'e', r'$\mathcal{E}$', cmap, vmin_e, vmax_e)
    make_plots(plotDir_vh, 2, np.divide(p_VH+t2pinn_VH,p_VH-t2pinn_VH/2), 'pLpT', r'$\mathcal{P}_{L}/\mathcal{P}_{\perp}$', cmap, 0, 1)
    make_plots(plotDir_vh, 3, reg_VH, 'reg', r'\mathrm{tanh}\rho/\rho$', cmap2, 0, 1)
    # Anisotropic hydro
    make_plots(plotDir_vah, 4, KnTaupi_VAH, 'KnTaupi', r'$\mathrm{Kn}_{\theta\pi}\equiv\theta\tau_{\pi}$', cmap, 0, 1)
    make_plots(plotDir_vah, 5, e_VAH, 'e', r'$\mathcal{E}$', cmap, vmin_e, vmax_e)
    make_plots(plotDir_vah, 6, np.divide(pl_VAH,(e_VAH-pl_VAH)/2), 'pLpT', r'$\mathcal{P}_{L}/\mathcal{P}_{\perp}$', cmap, 0, 1)
    make_plots(plotDir_vah, 7, fTSol_VAH, 'fTSol', r'$f_\perp(u_\perp)\neq 0\,\forall u_\perp\geq 0$', cmap, 0, 1)

    plt.show()
