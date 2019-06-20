import scipy.integrate as integrate
from numpy import sqrt, sin, cos, pi
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy import interpolate
from scipy.integrate import simps

dataDir = 'tests/output/qgp/glauber_2d_vh'

t, pitt = np.loadtxt(dataDir+'/pitt.dat', unpack=True)
t, pitx = np.loadtxt(dataDir+'/pitx.dat', unpack=True)
t, pity = np.loadtxt(dataDir+'/pity.dat', unpack=True)
t, pixx = np.loadtxt(dataDir+'/pixx.dat', unpack=True)
t, pixy = np.loadtxt(dataDir+'/pixy.dat', unpack=True)
t, piyy = np.loadtxt(dataDir+'/piyy.dat', unpack=True)
t, pinn = np.loadtxt(dataDir+'/pinn.dat', unpack=True)

labelSize = 16

plt.figure(1)
#plt.plot(t,pitt)
#plt.plot(t,pitx)
#plt.plot(t,pity)
plt.plot(t,pitt-pixx-piyy)
#plt.plot(t,pixx)
#plt.plot(t,pixy)
#plt.plot(t,piyy)
plt.plot(t,t**2*pinn)
plt.tight_layout(pad=0.1, h_pad=None, w_pad=None, rect=[0, 0, 1, 1])
plt.legend(loc='best', frameon=False, fontsize=labelSize)

plt.show()
