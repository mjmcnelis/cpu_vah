#!/usr/bin/env python3
from pylab import *

def load_var(dir, t, var, nx, ny, nz):
    x, y, n, dataRaw = np.loadtxt(dir + '/' + var + '_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    data = np.reshape(dataRaw, (nx, ny, nz))
    return squeeze(data)


if __name__ == '__main__':
    root = '/media/bazow/Data/fluid_dynamic_output_for_thesis/cpu-vah'
    dataDir = root + '/mcGlb_1d_pimunu_WTz_Pi_isotropicIC_QCDEOS_etaOverS0p2_reg1-10'
    dataDir = root + '/Glb_1d_pimunu_Pi_isotropicIC_QCDEOS_etaOverS0p2_reg1-10'
    outputDir = 'tests/output/qgp/validity'

    nx = 201
    ny = 1
    nz = 1
    dx = 0.1
    dy = 0.1
    dz = 0.1
    xmax = (nx - 1) / 2 * dx
    x = np.linspace(-xmax, xmax, num=nx)

    t0 = 0.5
    tf = 10.5
    dt = 0.05
    nt = int((tf-t0)/dt)+1
    t = np.linspace(t0, tf, num=nt)

    reg= np.zeros((nt, nx))
#    RPi= np.zeros((nt, nx))
    Rpi= np.zeros((nt, nx))
    fTSol= np.zeros((nt, nx))
#    KnTauPi= np.zeros((nt, nx))
    KnTaupi= np.zeros((nt, nx))
    e = np.zeros((nt, nx))

    for i in range(0, nt):
        ti = t[i]

        reg[i, :] = load_var(dataDir, ti, 'regulations', nx, ny, nz)
 #       RPi[i, :] = load_var(dataDir, ti, 'RPi', nx, ny, nz)
        Rpi[i, :] = load_var(dataDir, ti, 'Rpi', nx, ny, nz)
        fTSol[i, :] = load_var(dataDir, ti, 'fTSol_2', nx, ny, nz)
 #       KnTauPi[i, :] = load_var(dataDir, ti, 'knTauPi', nx, ny, nz)
        KnTaupi[i, :] = load_var(dataDir, ti, 'knTaupi', nx, ny, nz)
        e[i, :] = load_var(dataDir, ti, 'e', nx, ny, nz)

    np.savetxt(outputDir + '/reg.dat', np.transpose(reg), fmt="%.16f")
#    np.savetxt(outputDir + '/RPi.dat', np.transpose(RPi), fmt="%.16f")
    np.savetxt(outputDir + '/Rpi.dat', np.transpose(Rpi), fmt="%.16f")
    np.savetxt(outputDir + '/fTSol.dat', np.transpose(fTSol), fmt="%.16f")
#    np.savetxt(outputDir + '/KnTauPi.dat', np.transpose(KnTauPi), fmt="%.16f")
    np.savetxt(outputDir + '/KnTaupi.dat', np.transpose(KnTaupi), fmt="%.16f")
    np.savetxt(outputDir + '/e.dat', np.transpose(KnTaupi), fmt="%.16f")
