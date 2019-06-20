import matplotlib.pyplot as plt
from matplotlib.pylab import *

import equation_of_state as eos


def bulkViscosityToEntropyDensity(T):
    A1 = -13.77
    A2 = 27.55
    A3 = 13.45
    LAMBDA1 = 0.9
    LAMBDA2 = 0.25
    LAMBDA3 = 0.9
    LAMBDA4 = 0.22
    SIGMA1 = 0.025
    SIGMA2 = 0.13
    SIGMA3 = 0.0025
    SIGMA4 = 0.022
    x = T / 1.01355;
    zetabar = A1 * x * x + A2 * x - A3
    if x > 1.05:
        zetabar = LAMBDA1 * exp(-(x - 1) / SIGMA1) + LAMBDA2 * exp(-(x - 1) / SIGMA2) + 0.001
    elif x < 0.995:
        zetabar = LAMBDA3 * exp((x - 1) / SIGMA3) + LAMBDA4 * exp((x - 1) / SIGMA4) + 0.03
    return zetabar

def main():
    hbarc=0.197327

    T = np.linspace(0.75, 1.6, num=200)
    len=T.size
    e = np.zeros((len, 1))
    zetaBar=np.linspace(0.01, 2, num=200)
    zetaBarLower=np.linspace(0.01, 2, num=200)
    zetaBarUpper=np.linspace(0.01, 2, num=200)
    for i in range(0,len):
        Ti=T[i]
        e[i]=eos.equilibriumEnergyDensity(Ti)
        zetaBar[i]=bulkViscosityToEntropyDensity(Ti)
        zetaBarLower[i]=0.001*bulkViscosityToEntropyDensity(Ti)
        zetaBarUpper[i]=3*bulkViscosityToEntropyDensity(Ti)

    plotDir = 'tests/figs'

#    plt.style.use('bmh')
    #    mpl.rcParams['font.family'] = 'Ubuntu'
    #    plt.rcParams['text.color'] = 'k'
    #    plt.rcParams['xtick.color'] = 'k'
    #    plt.rcParams['ytick.color'] = 'k'
    #    plt.rcParams['axes.labelcolor'] = 'k'
#    plt.rcParams['axes.facecolor'] = 'white'
#    plt.rcParams['axes.grid'] = 'False'
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    labelSize = 16

    plt.figure(1)
    plt.plot(T*hbarc,zetaBar)
    plt.plot(T*hbarc,zetaBarLower,'k-', linewidth=0.5)
    plt.plot(T*hbarc,zetaBarUpper,'k-', linewidth=0.5)
    plt.fill_between(T*hbarc, zetaBarLower, zetaBarUpper, color='dodgerblue', alpha='0.15')
    plt.ylabel(r'$\zeta/s$', fontsize=labelSize)
    plt.xlabel('T [GeV]', fontsize=labelSize)
    plt.tick_params(top='off', right='off')
    plt.tight_layout()
    savefig(plotDir+'/specificBulkViscosityPlot.pdf', bbox_inches='tight')

    plt.show()

if __name__ == '__main__':
    main()



