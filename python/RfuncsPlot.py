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

def Rbar(a):
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

def Rtilde(a):
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
    return (-6.674731906076046e-6 + 0.004617789933500251*a + 0.7207562721999754*a2 + 9.097427250602184*a3 - 4.475814747302824*a4 - 36.37501529319408*a5 +
     46.868405146729316*a6 - 15.833867583743228*a7)/(
        0.06856675185266 + 2.9181587012768597*a + 11.951184087839218*a2 - 29.708257843442173*a3 - 2.618233802059826*a4 + 34.646239784689065*a5 -
     19.62596366454439*a6 + 2.374808442453899*a7)

def Rhat(a):
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
    return (0.0024792827625583747 + 1.943027171680747*a + 53.46970495217282*a2 + 19.989171951866325*a3 - 347.1285593126723*a4 + 412.2647882672885*a5 -
     140.53693383827797*a6)/(0.5061402347582388 + 29.466067530916984*a + 126.07947638942892*a2 - 334.420268508072*a3 + 86.57706367583984*a4 +
     183.53625188578846*a5 - 91.68259808111912*a6)

def RfuncsFromSecondOrderTransportCoefficients(a):
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

    # m^2 * \gamma^{(0)}_{-2,0,0,0}
    m2_R_Z2lnZ_g_0_m2000 = (-3.2479315686296945e-7 - 0.0024818436417066627*a - 1.6111984012761296*a2 - 127.73612771161369*a3 - 1008.9938164526582*a4 + 6967.254906666843*a5 -
     18324.948017843824*a6 + 25099.55628167157*a7 - 17270.922580416685*a8 + 4667.403001315362*a9)/ (
        9.890948008782298e-9 + 0.00016223744942607356*a + 0.1757897444655473*a2 + 20.44234240598682*a3 + 201.40183830725502*a4 - 1449.0506870352772*a5 +
     4018.683126779512*a6 - 5675.330193067643*a7 + 3959.9290468807912*a8 - 1076.2514264594497*a9)
    m2_R_Z2_g_0_m2000 = (-1.5559930868914822e-6 - 0.014981709826239044*a - 2180.1717646643933*a10 - 13.100395887657402*a2 - 1270.1133038412427*a3 + 6483.9906047206305*a4 -
     12792.244127833776*a5 + 6594.067183933232*a6 + 13564.962858734329*a7 - 22630.921246084006*a8 + 12243.54516222993*a9)/(
        4.9561108019411677e-8 + 0.0009652096929343721*a - 459.95391162217646*a10 + 1.3059657887008547*a2 + 171.9024694197664*a3 - 1040.5683376522913*a4 +
     2838.1416254430032*a5 - 4271.725089036778*a6 + 4052.707851388653*a7 - 2940.478561481007*a8 + 1648.6670223776266*a9)

    # m^(-2) * \gamma^{(0)}_{-2,0,4,0}
    mm2_R_0_g_0_m2040 =(0.29828706785393566 - 171.02025599125028*a - 14439.825541817447*a10 + 5056.403984218289*a2 - 32913.7554230461*a3 + 84210.31555291355*a4 - 79148.33364101285*a5 -
     49171.79759326144*a6 + 187041.18039643095*a7 - 182011.85508347143*a8 + 81548.38931795595*a9)/(
        127.28850795911887 - 733.4476008983003*a - 311.1521318940736*a10 + 1790.201010676183*a2 - 1994.5901505042752*a3 + 716.7012451042981*a4 -
     167.57357287690303*a5 + 821.3682997206282*a6 - 35.707655340707326*a7 - 1458.9355808622663*a8 + 1245.8476289159084*a9)
    mm2_R_Zm2_g_0_m2040 = (0.5503730699440683 + 116.12436316929069*a - 26724.912117163938*a10 - 12296.200600281796*a2 + 91390.00365346204*a3 - 305389.2113354753*a4 + 590964.5729874016*a5 -
     741756.2618562379*a6 + 644188.3945881774*a7 - 389945.72443247156*a8 + 149452.6643760904*a9)/ (
        212.4342766621021 - 972.205363010775*a - 452.4713348428893*a10 + 1551.6466075450276*a2 + 28.500987927598654*a3 - 2443.564066861641*a4 + 842.0321421896035*a5 +
     2382.6176950147824*a6 - 749.7704231352149*a7 - 2417.889536323947*a8 + 2018.669014833873*a9)
    mm2_R_LnZ_g_0_m2040 = (0.0020056853283532314 - 1.3550103585185087*a + 226.56390846375785*a10 - 44.449676579553895*a11 + 30.27602935864644*a2 - 165.3307114623337*a3 +
     357.8168925552034*a4 - 262.6276709033351*a5 - 155.2484363483795*a6 + 252.16638467851257*a7 + 151.3774157271213*a8 - 389.19113081630536*a9)/ (
        0.27036096424857836 - 1.3334833291297246*a - 0.25512595845135466*a10 + 0.07919845582019089*a11 + 2.606528106886374*a2 - 1.3981000048625831*a3 -
     1.882929319557247*a4 + 1.557102712196459*a5 + 1.437424923433941*a6 - 0.9689179954020858*a7 - 0.9074107943793904*a8 + 0.7953522391949045*a9)

    # m^2  \gamma^{(1)}_{-2,0,1,1}
    m2_R_2_g_1_m2011 = (0.00003063298763130019 + 0.010032700587292002*a - 0.11684107276326144*a2 - 13.34422667030812*a3 + 31.063751017661648*a4 - 11.296976111125272*a5 -
     20.364862088598894*a6 + 17.652385794980933*a7 - 3.6034202534002184*a8)/ (
        -0.0010721674876933715 - 0.1840204637199662*a + 3.7173267186679557*a2 + 207.90617859636498*a3 - 333.74586079841345*a4 - 226.95976766424184*a5 +
     660.0880880021133*a6 - 360.48927623536525*a7 + 49.694216832228015*a8)

    # m^2  \gamma^{(2)}_{-2,0,0,0}
    m2_R_2_g_2_m2000 = (0.0006580997236360711 + 0.30391474362309856*a + 15.042463311295156*a2 + 107.84132342439199*a3 - 231.79170934377856*a4 - 346.2926959207636*a5 +
     1181.56819743041*a6 - 1076.194286422169*a7 + 397.37538455708795*a8 - 47.85179965880513*a9)/ (
        0.13664546375128556 + 31.767030408850044*a + 968.3328377478067*a2 + 3478.9524703737707*a3 - 14681.332850820749*a4 + 3039.409829042078*a5 +
     32680.026784745856*a6 - 42704.12594765116*a7 + 20711.97612333118*a8 - 3525.1399850211346*a9)

    R_0_g_0_m2020 = (9.38649978486803e-6 + 0.06187055601680752*a + 12.75009989755989*a2 + 97.45033850687406*a3 - 773.8028659980835*a4 + 2312.148214043987*a5 - 3734.929336102941*a6 +
     3419.595878280452*a7 - 1687.5301976042288*a8 + 354.2562840062321*a9)/ (
        0.0006942401750132732 + 0.7753677259572115*a + 52.5797305539857*a2 + 80.6761035743086*a3 - 566.5441194717436*a4 - 224.2580378352535*a5 +
     4523.981522596095*a6 - 8733.812703192234*a7 + 6768.412948782072*a8 - 1901.8111492250273*a9)

    # \gamma^{(1)}_{-2,0,3,1}
    R_0_g_1_m2031 = (-1.0233483506421896e-7 - 0.005510394233958543*a - 0.5161308003737349*a2 - 4.115511461930346*a3 + 6.378431203946746*a4 + 3.926438723664259*a5 -
     8.465485699618803*a6 + 2.7925630611642154*a7)/ (0.001958161993306958 + 0.22517859370360388*a + 2.883216830325076*a2 + 0.1905363935371778*a3 - 12.192584184275201*a4 + 10.729468548753893*a5 -
     0.8635431725599291*a6 - 0.9690254375998808*a7)

    # \gamma^{(2)}_{-2,0,2,0}
    R_0_g_2_m2020 = (-6.632573230669451e-7 + 0.03213700809320489*a + 2.956808785961552*a2 + 25.852318177177626*a3 - 88.72830609932925*a4 + 91.14789347649632*a5 -
     34.03754591398221*a6 + 2.779419987826611*a7)/(
        0.18449848515718037 + 15.3508744250664*a + 115.64291110511445*a2 - 411.5208177159956*a3 + 400.73086985704936*a4 - 87.67260855227059*a5 -
     51.03507886676561*a6 + 18.32629361855308*a7)

    return {
            'm2_R_Z2lnZ_g_0_m2000':m2_R_Z2lnZ_g_0_m2000,
            'm2_R_Z2_g_0_m2000':m2_R_Z2_g_0_m2000,
            'mm2_R_0_g_0_m2040':mm2_R_0_g_0_m2040,
            'mm2_R_Zm2_g_0_m2040':mm2_R_Zm2_g_0_m2040,
            'mm2_R_LnZ_g_0_m2040':mm2_R_LnZ_g_0_m2040,
            'm2_R_2_g_1_m2011':m2_R_2_g_1_m2011,
            'm2_R_2_g_2_m2000':m2_R_2_g_2_m2000,
            'R_0_g_0_m2020':R_0_g_0_m2020,
            'R_0_g_1_m2031':R_0_g_1_m2031,
            'R_0_g_2_m2020':R_0_g_2_m2020
            }

if __name__ == '__main__':
    na = 1000
    a = np.linspace(0, 1, num=na)
    a2 = np.linspace(0.001, .999, num=na)

    plotDir = 'tests/figs/qgp/'

    rxi = zeros((na,1))
    rbar = zeros((na,1))
    rtilde = zeros((na,1))
    rhat = zeros((na,1))
    # R functions from second-order transport coefficients
    m2_R_Z2lnZ_g_0_m2000 = zeros((na,1))
    m2_R_Z2_g_0_m2000 = zeros((na,1))
    mm2_R_0_g_0_m2040 = zeros((na,1))
    mm2_R_Zm2_g_0_m2040 = zeros((na,1))
    mm2_R_LnZ_g_0_m2040 = zeros((na,1))
    m2_R_2_g_1_m2011 = zeros((na,1))
    m2_R_2_g_2_m2000 = zeros((na,1))

    R_0_g_0_m2020 = zeros((na,1))
    R_0_g_1_m2031 = zeros((na,1))
    R_0_g_2_m2020 = zeros((na,1))

    for i in range(0, na):
        ai = a[i]
        rxi[i] = ai
        rbar[i] = Rbar(ai)
        rtilde[i] = Rtilde(ai)
        rhat[i] = Rhat(ai)

        # second-order
        ai2 = a2[i]
        RsecondOrder = RfuncsFromSecondOrderTransportCoefficients(ai2)
        m2_R_Z2lnZ_g_0_m2000[i] = RsecondOrder['m2_R_Z2lnZ_g_0_m2000']
        m2_R_Z2_g_0_m2000[i] = RsecondOrder['m2_R_Z2_g_0_m2000']
        mm2_R_0_g_0_m2040[i] = RsecondOrder['mm2_R_0_g_0_m2040']
        mm2_R_Zm2_g_0_m2040[i] = RsecondOrder['mm2_R_Zm2_g_0_m2040']
        mm2_R_LnZ_g_0_m2040[i] = RsecondOrder['mm2_R_LnZ_g_0_m2040']
        m2_R_2_g_1_m2011[i] = RsecondOrder['m2_R_2_g_1_m2011']
        m2_R_2_g_2_m2000[i] = RsecondOrder['m2_R_2_g_2_m2000']

        R_0_g_0_m2020[i] = RsecondOrder['R_0_g_0_m2020']
        R_0_g_1_m2031[i] = RsecondOrder['R_0_g_1_m2031']
        R_0_g_2_m2020[i] = RsecondOrder['R_0_g_2_m2020']

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

    fig, ax = plt.subplots()
    ax.plot(a, rxi, color='black', linewidth=3.5, linestyle='-', label=r'$\mathcal{R}_\xi$')
    ax.plot(a, rbar, color='red', linewidth=3.5, linestyle='--', label=r'$\bar{\mathcal{R}}$')
    ax.plot(a, rtilde, color='blue', linewidth=3.5, linestyle='-.', label=r'$\tilde{\mathcal{R}}$')
    ax.plot(a, rhat, color='green', linewidth=3.5, linestyle=':', label=r'$\hat{\mathcal{R}}$')
    ax.axvspan(0, 0.33333, alpha=0.25, color='dodgerblue')
    pylab.xlim([a[0],a[na-1]])
    #xticks(np.arange(-10, 11, 5))
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\hat{\mathcal{P}}_L/\mathcal{E}$')
    plt.ylabel('')
    plt.legend(loc='best', frameon=False)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    savefig(plotDir+'/RfuncsPlot.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    ## split into two
    fig, ax = plt.subplots()
    ax.plot(a2, mm2_R_Zm2_g_0_m2040, color='green', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,4,0}^{(0)}$')
    ax.plot(a2, mm2_R_LnZ_g_0_m2040, color='magenta', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(2,1)})_{-2,0,4,0}^{(0)}$')
    ax.plot(a2, mm2_R_0_g_0_m2040, color='purple', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(2,0)})_{-2,0,4,0}^{(0)}$')
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\hat{\mathcal{P}}_L/\mathcal{E}$')
    plt.ylabel('')
    plt.legend(loc='best', frameon=False)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    #savefig(plotDir+'/RfuncsPlot.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    line1, = ax.plot(a2, R_0_g_0_m2020, color='orange', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,2,0}^{(0)}$')
    line2, = ax.plot(a2, m2_R_2_g_1_m2011, color='salmon', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,1,1}^{(1)}$')
    line3, = ax.plot(a2, R_0_g_1_m2031, color='blue', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,3,1}^{(1)}$')
    plt.legend(handles=[line1,line2,line3], loc='lower left', frameon=False)
    line4, = ax.plot(a2, m2_R_2_g_2_m2000, color='brown', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,0,0}^{(2)}$')
    line5, = ax.plot(a2, R_0_g_2_m2020, color='red', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,2,0}^{(2)}$')
    line6, = ax.plot(a2, m2_R_Z2_g_0_m2000, color='brown', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,0,0}^{(2)}$')
    #plt.legend(handles=[line4,line5,line6], loc='lower right', frameon=False)
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\hat{\mathcal{P}}_L/\mathcal{E}$')
    plt.ylabel('')
    #plt.legend(loc='best', frameon=False)
    #plt.legend(handles=[line1,line2,line3], loc='lower right', frameon=False)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    #savefig(plotDir+'/RfuncsPlot.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    fig, ax = plt.subplots()
    # -2,0,2,0
    ax.plot(a2, R_0_g_0_m2020, color='orange', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,2,0}^{(0)}$')
    ax.plot(a2, m2_R_2_g_1_m2011, color='salmon', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,1,1}^{(1)}$')
    ax.plot(a2, m2_R_2_g_2_m2000, color='brown', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,0,0}^{(2)}$')
    ax.plot(a2, R_0_g_1_m2031, color='blue', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,3,1}^{(1)}$')
    ax.plot(a2, R_0_g_2_m2020, color='red', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,2,0}^{(2)}$')
    # R_{-2,0,4,0}
    ax.plot(a2, mm2_R_Zm2_g_0_m2040, color='green', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,4,0}^{(0)}$')
    ax.plot(a2, mm2_R_LnZ_g_0_m2040, color='magenta', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(2,1)})_{-2,0,4,0}^{(0)}$')
    ax.plot(a2, mm2_R_0_g_0_m2040, color='purple', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(2,0)})_{-2,0,4,0}^{(0)}$')
    # R_{-2,0,0,0}^{(2)}
    ax.plot(a2, m2_R_2_g_2_m2000, color='brown', linewidth=3.5, linestyle='-', label=r'$(\mathcal{R}_\gamma^{(0,0)})_{-2,0,0,0}^{(2)}$')
    ax.xaxis.set_minor_locator(minorLocator)
    plt.xlabel(r'$\hat{\mathcal{P}}_L/\mathcal{E}$')
    plt.ylabel('')
    plt.legend(loc='best', frameon=False)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)
    #savefig(plotDir+'/RfuncsPlot.pdf', pad=pad, h_pad=h_pad, w_pad=w_pad, rect=rect)

    plt.show()
