import scipy.integrate as integrate
from numpy import sqrt, sin, cos, pi
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy import interpolate
from scipy.integrate import simps

def load_var(dir, t, var, nx, ny, nz, xx, yy):
    x, y, n, dataRaw = np.loadtxt(dir + '/' + var + '_' + '{:.3f}'.format(t) + '.dat', unpack=True)
    data = np.reshape(dataRaw, (nx, ny, nz))
    dataInt = interpolate.interp2d(x, y, squeeze(data), kind='cubic')
    return dataInt(xx,yy)

root = '/media/bazow/Data/fluid_dynamic_output_for_thesis/'
type_vh = 'gpu-vh/'
tag = 'Glb_2d_shear_isotropicIC_conformalEOS_etaOverS0p2'
dataDir_vh = root + type_vh + tag

outputDir = 'tests/output/qgp/glauber_2d_vh'

nx = 161
ny = 161
nz = 1
dx = 0.1
dy = 0.1
dz = 0.1
xmax = (nx - 1) / 2 * dx
x = np.linspace(-xmax, xmax, num=nx)
y = np.linspace(-xmax, xmax, num=nx)

nxx = 1*nx
xx = np.linspace(-xmax, xmax, num=nxx)
yy = np.linspace(-xmax, xmax, num=nxx)

t0 = 0.5
tf = 1.0
dt = 0.05
nt = int((tf - t0) / dt) + 1
t = np.linspace(t0, tf, num=nt)

pitt = np.zeros((nt, 1))
pitx = np.zeros((nt, 1))
pity = np.zeros((nt, 1))
pixx = np.zeros((nt, 1))
pixy = np.zeros((nt, 1))
piyy = np.zeros((nt, 1))
pinn = np.zeros((nt, 1))

for i in range(0, nt):
    ti = t[i]
    e2d = load_var(dataDir_vh, ti, 'e', nx, ny, nz, xx, yy)
    p2d = load_var(dataDir_vh, ti, 'p', nx, ny, nz, xx, yy)
    pitt2d = load_var(dataDir_vh, ti, 'pitt', nx, ny, nz, xx, yy)
    pitx2d = load_var(dataDir_vh, ti, 'pitx', nx, ny, nz, xx, yy)
    pity2d = load_var(dataDir_vh, ti, 'pity', nx, ny, nz, xx, yy)
    pixx2d = load_var(dataDir_vh, ti, 'pixx', nx, ny, nz, xx, yy)
    pixy2d = load_var(dataDir_vh, ti, 'pixy', nx, ny, nz, xx, yy)
    piyy2d = load_var(dataDir_vh, ti, 'piyy', nx, ny, nz, xx, yy)
    pinn2d = load_var(dataDir_vh, ti, 'pinn', nx, ny, nz, xx, yy)

    pitt[i] = simps(simps(np.divide(pitt2d,e2d+p2d), yy), xx)/(2*xmax)/(2*xmax)
    pitx[i] = simps(simps(np.divide(pitx2d,e2d+p2d), yy), xx)/(2*xmax)/(2*xmax)
    pity[i] = simps(simps(np.divide(pity2d,e2d+p2d), yy), xx)/(2*xmax)/(2*xmax)
    pixx[i] = simps(simps(np.divide(pixx2d,e2d+p2d), yy), xx)/(2*xmax)/(2*xmax)
    pixy[i] = simps(simps(np.divide(pixy2d,e2d+p2d), yy), xx)/(2*xmax)/(2*xmax)
    piyy[i] = simps(simps(np.divide(piyy2d,e2d+p2d), yy), xx)/(2*xmax)/(2*xmax)
    pinn[i] = simps(simps(np.divide(pinn2d,e2d+p2d), yy), xx)/(2*xmax)/(2*xmax)

np.savetxt(outputDir+'/pitt.dat', np.c_[t,pitt], fmt="%.16f")
np.savetxt(outputDir+'/pitx.dat', np.c_[t,pitx], fmt="%.16f")
np.savetxt(outputDir+'/pity.dat', np.c_[t,pity], fmt="%.16f")
np.savetxt(outputDir+'/pixx.dat', np.c_[t,pixx], fmt="%.16f")
np.savetxt(outputDir+'/pixy.dat', np.c_[t,pixy], fmt="%.16f")
np.savetxt(outputDir+'/piyy.dat', np.c_[t,piyy], fmt="%.16f")
np.savetxt(outputDir+'/pinn.dat', np.c_[t,pinn], fmt="%.16f")

#plt.figure(1)
#plt.plot(t,pitt,t,pitx,t,pity,t,pixx,t,pixy,t,piyy,t,pinn)

#plt.show()
