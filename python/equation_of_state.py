from scipy import integrate
from matplotlib.pylab import *

def Power(a, b):
    return pow(a, b)


def equilibriumPressure(e):
    return (-0.25181736420168666 + 9737.845799644809 * e + 1.077580993288114e6 * Power(e, 2) + 3.1729694865420084e6 * Power(e,3) + 1.6357487344679043e6 * Power(e, 4) +
     334334.4309240126 * Power(e, 5) + 41913.439282708554 * Power(e, 6) + 6340.448389300905 * Power(e,7) + 141.5073484468774 * Power(e, 8) +
     0.7158279081255019 * Power(e, 9) + 0.0009417586777847889 * Power(e, 10) + 3.1188455176941583e-7 * Power(e,11) + 1.9531729608963267e-11 * Power(e, 12)) / \
    (45829.44617893836 + 4.0574329080826794e6 * e + 2.0931169138134286e7 * Power(e, 2) + 1.3512402226067686e7 * Power(e,3) + 1.7851642641834426e6 * Power(e, 4) +
     278581.2989342773 * Power(e, 5) + 26452.34905933697 * Power(e, 6) + 499.04919730607065 * Power(e,7) + 2.3405487982094204 * Power(e, 8) +
     0.002962497695527404 * Power(e, 9) + 9.601103399348206e-7 * Power(e, 10) + 5.928138360995685e-11 * Power(e,11) - 3.2581066229887368e-18 * Power(e, 12))

def speedOfSoundSquared(e):
    return (1.2625814949136498e-42 + 2.4439030393295734e-33 * e + 9.322641344497548e-26 * Power(e, 2) +
            1.9540661873302495e-19 * Power(e, 3) + 5.0244280894342345e-14 * Power(e, 4) +
            2.369176131226614e-9 * Power(e, 5) + 0.00002467595835925649 * Power(e, 6) +
            0.0608709022174928 * Power(e, 7) + 35.13582886022962 * Power(e, 8) +
            4222.887394189519 * Power(e, 9) + 72483.09074617726 * Power(e, 10) +
            417317.67297479673 * Power(e, 11) + 585334.4710663276 * Power(e, 12) +
            57539.2198687339 * Power(e, 13) + 14441.972045822631 * Power(e, 14) +
            1216.0067118006687 * Power(e, 15) + 13.489236656837132 * Power(e, 16) +
            14.786577881305382 * Power(e, 17) - 2.5565796031645758 * Power(e, 18) +
            0.05603218920634199 * Power(e, 19) + 0.000984102224155054 * Power(e, 20) +
            2.1161830800816893e-6 * Power(e, 21) + 8.101199742688298e-10 * Power(e, 22) +
            3.89552175647867e-14 * Power(e, 23)) / (
           3.56558833745313e-41 + 4.6403916630291216e-32 * e + 1.246926715508643e-24 * Power(e, 2) +
           1.9406804331642527e-18 * Power(e, 3) + 3.8792079821435393e-13 * Power(e, 4) +
           1.473941085930358e-8 * Power(e, 5) + 0.0001283715493908053 * Power(e, 6) +
           0.27578987407916083 * Power(e, 7) + 144.29417387657588 * Power(e, 8) +
           16436.520976883963 * Power(e, 9) + 315466.19533599314 * Power(e, 10) +
           2.8901807318940037e6 * Power(e, 11) + 5.506960426273549e6 * Power(e, 12) -
           355776.5733809447 * Power(e, 13) + 156209.84678369848 * Power(e, 14) -
           5701.991686556281 * Power(e, 15) + 832.8660760210314 * Power(e, 16) +
           1.9784144986425583 * Power(e, 17) - 7.706332342730179 * Power(e, 18) +
           0.20216864185984038 * Power(e, 19) + 0.00313721024057907 * Power(e, 20) +
           6.541440626577012e-6 * Power(e, 21) + 2.4689937881750025e-9 * Power(e, 22) +
           1.1758507795621395e-13 * Power(e, 23))


def effectiveTemperature(e):
    return (-4.873400811469563e-10 - 0.33967898165285937 * e - 243.87315076666465 * Power(e, 2) -
            7950.4417656525575 * Power(e, 3) - 15830.731725359203 * Power(e, 4) -
            1638.2165592489832 * Power(e, 5) - 15.708990537469983 * Power(e, 6) -
            0.017506852476184522 * Power(e, 7) - 1.836593645164133e-6 * Power(e, 8)) / (
           -4.873378567876981e-6 - 1.9735695129421675 * e - 691.5144322242182 * Power(e, 2) -
           13185.130839434429 * Power(e, 3) - 20098.669864541644 * Power(e, 4) -
           1169.8522357611525 * Power(e, 5) - 6.630649855145946 * Power(e, 6) -
           0.004312673981317225 * Power(e, 7) - 2.1655522622011757e-7 * Power(e, 8))


def equilibriumEnergyDensity(T):
    return (-0.011958188410851651 + 119.89423098138208 * T - 3156.9475699248055 * Power(T, 2) +
            32732.86844374939 * Power(T, 3) - 187899.8994764422 * Power(T, 4) + 712537.3610845465 * Power(T, 5) -
            1.557049803609345e6 * Power(T, 6) + 1.4852519861308339e6 * Power(T, 7) +
            532132.6079941876 * Power(T, 8) - 1.963099445042592e6 * Power(T, 9) -
            4484.44579242679 * Power(T, 10) + 1.7984228830058286e6 * Power(T, 11) +
            119345.25619517374 * Power(T, 12) - 1.3499773937058165e6 * Power(T, 13) -
            207838.4995663606 * Power(T, 14) + 654970.2138652403 * Power(T, 15) -
            78643.00334616247 * Power(T, 16) + 40274.00078068926 * Power(T, 17) +
            422619.58977657766 * Power(T, 18) - 409688.07836393174 * Power(T, 19) -
            62005.75915066359 * Power(T, 20) + 46788.14270090656 * Power(T, 21) +
            40784.330477857235 * Power(T, 22) - 12589.47744840392 * Power(T, 23)) / (
           31630.074365558292 - 127100.88940643385 * T + 173528.1225422275 * Power(T, 2) -
           39403.297956865215 * Power(T, 3) - 85582.57873541754 * Power(T, 4) +
           9320.560804233442 * Power(T, 5) + 50882.74198960172 * Power(T, 6) +
           20335.926473421183 * Power(T, 7) - 14897.725710713818 * Power(T, 8) - 23836.484117457 * Power(T, 9) -
           13726.013896090335 * Power(T, 10) + 4517.908673107615 * Power(T, 11) +
           18056.19917986404 * Power(T, 12) + 14954.82860467155 * Power(T, 13) +
           2569.623976952738 * Power(T, 14) - 9304.046211514986 * Power(T, 15) -
           15606.429173842751 * Power(T, 16) + 8383.710735812094 * Power(T, 17) +
           1591.3177623932843 * Power(T, 18) - 678.748230997762 * Power(T, 19) -
           33.58687934953277 * Power(T, 20) + 3.2520554133126285 * Power(T, 21) -
           0.19647288043440464 * Power(T, 22) + 0.005443394551264717 * Power(T, 23))

def main():
    hbarc = 0.197327

    Nf = 3
    eg = 16 * math.pi ** 2 / 30
    eq = 6 * 3 * 7 * math.pi ** 2 / 120
    sFac = eg + eq

    len = 500
    T = np.linspace(0.1, 5, num=len)
    e = np.zeros((len, 1))
    p = np.linspace(0.1, 5, num=len)
    cs2 = np.linspace(0.1, 5, num=len)
    for i in range(0, len):
        Ti = T[i]
        ei = equilibriumEnergyDensity(Ti)
        e[i] = ei / (sFac * pow(Ti, 4))
        p[i] = equilibriumPressure(ei) / (sFac * pow(Ti, 4) / 3)
        cs2[i] = speedOfSoundSquared(ei)

    hbarc = 0.197327

    cs2SB = np.ones((len, 1)) / 3

#    plt.style.use('bmh')
    #    mpl.rcParams['font.family'] = 'Ubuntu'
    #    plt.rcParams['text.color'] = 'k'
    #    plt.rcParams['xtick.color'] = 'k'
    #    plt.rcParams['ytick.color'] = 'k'
    #    plt.rcParams['axes.labelcolor'] = 'k'
#    plt.rcParams['axes.facecolor'] = 'white'
#    plt.rcParams['axes.grid'] = 'False'
    #matplotlib.rcParams.update({'font.size': 14})
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    plotDir = 'tests/figs'

    labelSize = 16

    plt.figure(1)
    plt.plot(T * hbarc, e, label=r'$\mathcal{E}/\mathcal{E}_\mathrm{ideal}$')
    plt.plot(T * hbarc, p, 'r-',label=r'$\mathcal{P}_{0}/\mathcal{P}_\mathrm{ideal}$')
    plt.xlabel('T [GeV]', fontsize=labelSize)
    plt.tick_params(top='off', right='off')
    plt.tight_layout(pad=0.1, h_pad=None, w_pad=None, rect=[0, 0, 1, 1])
    plt.legend(loc='best', frameon=False, fontsize=22)
    savefig(plotDir + '/pressureEoS.pdf', pad=0.1, h_pad=None, w_pad=None, rect=[0, 0, 1, 1])

    plt.figure(2)
    plt.plot(T * hbarc, cs2)
    plt.plot(T * hbarc, cs2SB, 'k-', label='Stefan-Boltzmann limit')
    plt.ylabel(r'$c^{2}_{s}$', fontsize=22)
    plt.xlabel('T [GeV]', fontsize=labelSize)
    plt.tick_params(top='off', right='off')
    plt.tight_layout(pad=0.1, h_pad=None, w_pad=None, rect=[0, 0, 1, 1])
    plt.legend(loc='best', frameon=False, fontsize=22)
    savefig(plotDir + '/speedOfSoundSquared.pdf', pad=0.1, h_pad=None, w_pad=None, rect=[0, 0, 1, 1])

    plt.show()

if __name__ == '__main__':
    main()
