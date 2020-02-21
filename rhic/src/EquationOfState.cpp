#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include "../include/Macros.h"
#include "../include/Precision.h"
#include "../include/EquationOfState.h"

precision energy_density_cutoff(precision e_min, precision e)
{
	precision e_cut = fmax(0., e);

	return e_cut  +  e_min * exp(- e_cut / e_min);	// regulated energy density asymptotes to e_min
													// as e -> 0 (avoids discontinuites in energy profile)
	//return fmax(e_min, e);						// hard cutoff
}




equation_of_state_new::equation_of_state_new(precision e1_in, precision conformal_prefactor_in)
{
	e1 = e1_in;
	conformal_prefactor = conformal_prefactor_in;

	precision hotqcd_e_min = 0.00175;		// min energy density in fm^-4 rounded up to 3 SFs (see notebook in eos/hotqcd_smash)

	if(e1 < hotqcd_e_min)
	{
		printf("equation_of_state_new error: e1 = %.3f is smaller than minimum energy density = %lf in hotqcd eos table\n", e1, hotqcd_e_min);
		exit(-1);
	}

#ifdef LATTICE_QCD
#ifndef CONFORMAL_EOS
	precision e2  = e1  * e1;
	precision e3  = e2  * e1;
	precision e4  = e3  * e1;
	precision e5  = e4  * e1;
	precision e6  = e5  * e1;
	precision e7  = e6  * e1;
	precision e8  = e7  * e1;
	precision e9  = e8  * e1;
	precision e10 = e9  * e1;
	precision e11 = e10 * e1;
	precision e12 = e11 * e1;
	precision e13 = e12 * e1;
	precision e14 = e13 * e1;
	precision e15 = e14 * e1;
	precision e16 = e15 * e1;
	precision e17 = e16 * e1;
	precision e18 = e17 * e1;

	T = (0.000017501666790163825 + 0.03699450264229552*e1 + 9.53079264346692*e2 + 478.0574558421746*e3 + 4815.265814410284*e4 + 6071.818762824047*e5 - 3242.07917034678*e6 + 3586.8908905829185*e7 - 1492.4308460848283*e8 + 260.90789608163226*e9 + 2.9614572565180657*e10 - 5.4725878573028615*e11 + 0.3606239543966032*e12 + 0.01307881666327972*e13 + 0.00008243332949376859*e14 + 1.3829509318667127e-7*e15 + 6.517489614388798e-11*e16 + 7.742208716620397e-15*e17 + 1.5150865275249154e-19*e18)/(0.00011098075835081602 + 0.14893326363389814*e1 + 27.151183301353235*e2 + 1004.4336510123287*e3 + 7661.3372953228045*e4 + 6316.605157020178*e5 - 3679.435967309141*e6 + 4512.891498218059*e7 - 2114.5174086096513*e8 + 464.103808135828*e9 - 30.418708324615274*e10 - 3.133493850682989*e11 + 0.3625655777300751*e12 + 0.007981157499556196*e13 + 0.000036132179693799986*e14 + 4.4597353454908245e-8*e15 + 1.525105645624315e-11*e16 + 1.234200041943215e-15*e17 + 1.2919403060748804e-20*e18);
#endif
#endif

#ifdef CONFORMAL_EOS
#ifndef LATTICE_QCD
	T = pow(e1 / conformal_prefactor, 0.25);
#endif
#endif

#ifdef CONFORMAL_EOS
#ifdef LATTICE_QCD
	printf("equation_of_state_new::equation_of_state_new error: not eos switch yet\n");
	exit(-1);
#endif
#endif

	T1 = T;
	T2 = T1 * T1;
	T3 = T2 * T1;
	T4 = T3 * T1;
	T5 = T4 * T1;
	T6 = T5 * T1;
	T7 = T6 * T1;
	T8 = T7 * T1;
	T9 = T8 * T1;
	T10 = T9 * T1;
	T11 = T10 * T1;
	T12 = T11 * T1;
	T13 = T12 * T1;
	T14 = T13 * T1;
	T15 = T14 * T1;
	T16 = T15 * T1;
	T17 = T16 * T1;
	T18 = T17 * T1;
	T19 = T18 * T1;
	T20 = T19 * T1;
	T21 = T20 * T1;
	T22 = T21 * T1;
}


equation_of_state_new::~equation_of_state_new()
{

}


precision equation_of_state_new::equilibrium_pressure()
{
	precision p;

#ifdef LATTICE_QCD
#ifndef CONFORMAL_EOS
	p = (380.1314073917778 - 7680.033526326507*T1 + 68602.08760067614*T2 - 359106.2417553722*T3 + 1.2389284083604014e6*T4 - 3.04440599633009e6*T5 + 5.709246989662066e6*T6 - 8.66341866186127e6*T7 + 1.0835954146449104e7*T8 - 1.0664248361631993e7*T9 + 7.473285898937746e6*T10 - 3.213922187759718e6*T11 + 626648.5339141666*T12)/(-13940.851869899469 + 87266.66293844217*T1 - 186446.89426907562*T2 + 77592.80031839639*T3 + 373353.71023107687*T4 - 823775.2162777347*T5 + 813419.0979858796*T6 - 413342.5682364177*T7 + 79476.33235382187*T8 + 7351.581336605625*T9 - 829.7870887229601*T10 + 55.35878026074939*T11 - 1.6527010061152714*T12);
#endif
#endif

#ifdef CONFORMAL_EOS
#ifndef LATTICE_QCD
	p = e1 / 3.;
#endif
#endif

#ifdef CONFORMAL_EOS
#ifdef LATTICE_QCD
	printf("equation_of_state_new::equilibrium_pressure error: not eos switch yet\n");
	exit(-1);
#endif
#endif

	if(p < 0)
	{
		printf("equation_of_state_new::equilibrium_pressure error: p = %lf is negative. Enforcing positive equilibrium_pressure\n", p);
		p = fmax(p, 0.);
	}
	return p;
}


precision equation_of_state_new::speed_of_sound_squared()
{
	precision cs2;

#ifdef LATTICE_QCD
#ifndef CONFORMAL_EOS
	cs2 = (-0.28269533732608904 + 9.795100529341218*T1 - 88.02039352141522*T2 + 396.3568943492647*T3 - 1054.1593908781763*T4 + 1675.8066832885186*T5 - 1266.9814782940746*T6 - 550.0957286372637*T7 + 2108.7175083550146*T8 - 1060.4095987288595*T9 - 2032.4961002057378*T10 + 3457.9927779991444*T11 - 1135.821740237249*T12 - 2487.542017457926*T13 + 3922.812440884798*T14 - 2826.984002733641*T15 + 1179.041785013648*T16 - 276.090105534535*T17 + 28.36006120492868*T18)/(1.8215578803318528 - 1.919551715217968*T1 - 108.04292192746387*T2 + 804.1669648293796*T3 - 2897.7934908524808*T4 + 6237.002516638373*T5 - 7989.399116914714*T6 + 4478.6226678142475*T7 + 2977.837996285879*T8 - 7263.759282392347*T9 + 4854.78645806174*T10 - 1936.5615884900997*T11 + 5184.0221841281555*T12 - 11356.100393452223*T13 + 12838.831974996214*T14 - 8611.788755995802*T15 + 3528.5290442427595*T16 - 825.5150041085104*T17 + 85.25874124597529*T18);
#endif
#endif

#ifdef CONFORMAL_EOS
#ifndef LATTICE_QCD
	cs2 = 1./3.;
#endif
#endif

#ifdef CONFORMAL_EOS
#ifdef LATTICE_QCD
	printf("equation_of_state_new::speed_of_sound_squared error: not eos switch yet\n");
	exit(-1);
#endif
#endif

	return cs2;
}


precision equation_of_state_new::z_quasi()
{
	precision z;

#ifdef LATTICE_QCD
#ifndef CONFORMAL_EOS
	z = (1.1378463937729772 - 3.7429882989442396*T1 + 4.760662889790829*T2 - 61.47601935330168*T3 + 362.9948908884739*T4 - 900.3019397857224*T5 + 981.1920460850209*T6 + 6.4475488614815*T7 - 1094.8979490604052*T8 + 625.5675319117105*T9 + 1026.6426927177527*T10 - 1730.4567564488684*T11 + 746.597311849036*T12 + 590.5360262740035*T13 - 1134.5106896432587*T14 + 973.9121475422542*T15 - 584.1714083515602*T16 + 240.85064317284673*T17 - 53.24312520325643*T18 + 0.6733469377693581*T19 + 1.4881828895072788*T20)/(-0.004733051795027734 + 1.7487733955787945*T1 - 13.540592604152945*T2 + 45.07450770996357*T3 - 83.87976680046712*T4 + 99.00413225938306*T5 - 84.29839809728547*T6 + 41.66125043072567*T7 + 62.06898875780188*T8 - 200.04432376326648*T9 + 235.79807242819783*T10 - 137.42412960602465*T11 + 50.68048704490057*T12 - 28.304345373554163*T13 - 13.825711237747374*T14 + 53.490876334068545*T15 - 11.69543394621012*T16 - 54.49721718220276*T17 + 58.87066327140601*T18 - 24.85781073302625*T19 + 3.974711605201351*T20);

	// z = (0.08303600003873625 - 0.31663876548255276*T1 - 0.3437654480323534*T2 + 3.1935438876025226*T3 - 1.8505414089256629*T4 - 20.103250346043314*T5 + 65.63805450196016*T6 - 102.7208024206455*T7 + 94.78859714874518*T8 - 49.90533697850738*T9 + 9.46999942889758*T10 + 4.6616791100921775*T11 - 3.1325717479601183*T12 + 0.521050748486101*T13 + 0.016946659040681687*T14)/(-0.0013528615203112476 + 0.16110349970479615*T1 - 1.5724445522755603*T2 + 7.583977032359246*T3 - 23.86307662651365*T4 + 55.24608845445409*T5 - 99.50129972931308*T6 + 141.4504474648258*T7 - 156.72733620989175*T8 + 131.3886248621511*T9 - 79.65025242023526*T10 + 32.4639587790247*T11 - 7.690588987983213*T12 + 0.654617821652585*T13 + 0.05753361056246534*T14);
#endif
#endif

#ifdef CONFORMAL_EOS
#ifndef LATTICE_QCD
	z = 0;
#endif
#endif

#ifdef CONFORMAL_EOS
#ifdef LATTICE_QCD
	printf("equation_of_state_new::z_quasi error: not eos switch yet\n");
	exit(-1);
#endif
#endif

	return z;
}


precision equation_of_state_new::mdmde_quasi()
{
	precision mdmde;

#ifdef LATTICE_QCD
#ifndef CONFORMAL_EOS
	mdmde = (1.3115525644293986 - 22.285338262587484*T1 + 144.68825773864012*T2 - 438.8278208945311*T3 + 498.94435930741*T4 + 507.0920663664031*T5 - 1547.92410636253*T6 - 1504.4383935281808*T7 + 9297.376947138675*T8 - 12175.92946074074*T9 + 885.018011913933*T10 + 15602.4846429165*T11 - 17542.763995448786*T12 + 1019.2211264838179*T13 + 17349.60570188111*T14 - 22746.28995274643*T15 + 16476.41163700659*T16 - 7910.216620082677*T17 + 2619.820733941211*T18 - 594.3069839126013*T19 + 88.54508224378642*T20 - 7.847075159560508*T21 + 0.3096271428981808*T22)/(-0.0007931211446189383 + 0.011082456284504633*T1 - 0.0817402082737028*T2 + 1.9694854631078555*T3 - 33.483537869885936*T4 + 298.7440672009308*T5 - 1633.2195758462135*T6 + 6022.525486566549*T7 - 15728.60427066111*T8 + 29581.17081929427*T9 - 39387.13461895701*T10 + 33810.14418055694*T11 - 9869.992893482482*T12 - 21563.545449965684*T13 + 44913.94997468314*T14 - 51453.21757425022*T15 + 41356.76207079614*T16 - 21857.986191619893*T17 + 4472.975794963638*T18 + 3165.6425078859966*T19 - 2985.8721816685843*T20 + 1025.9863871663451*T21 - 136.74301958608172*T22);

	// mdmde = (72.6878665548801 - 382.52170027983834*T1 + 644.3762319130602*T2 - 129.19135072153108*T3 - 634.5356306099509*T4 + 356.5209109636042*T5 + 324.89682684648994*T6 - 478.88942323317116*T7 + 332.9743774340843*T8 + 618.0973965604695*T9 - 1171.439436577495*T10 - 628.7476208107144*T11 + 1515.062832254186*T12 + 835.8529297487167*T13 - 1672.7519285091576*T14 - 840.4005410364722*T15 + 1925.9539310948414*T16 - 89.64113010803419*T17 - 1359.2904984693519*T18 + 1075.955079416841*T19 - 372.5115164916778*T20 + 61.625265972066614*T21 - 4.083321698953611*T22)/(0.4136588243656204 - 10.885009497755833*T1 + 117.16425450104548*T2 - 672.2563826169645*T3 + 2364.751507966597*T4 - 5211.1359137398*T5 + 6490.575169088332*T6 - 2208.3477746327376*T7 - 4654.8752588888465*T8 + 3251.1047628382885*T9 + 5476.1400511200545*T10 - 3576.274191301569*T11 - 7423.831257839161*T12 + 2443.1019713623186*T13 + 10262.74306696526*T14 - 385.2515628289653*T15 - 12417.625007308698*T16 - 1873.514054527994*T17 + 14520.545047537851*T18 - 1021.296289119715*T19 - 12124.52644556774*T20 + 8427.140010025205*T21 - 1773.8514147774117*T22);
#endif
#endif

#ifdef CONFORMAL_EOS
#ifndef LATTICE_QCD
	mdmde = 0;
#endif
#endif

#ifdef CONFORMAL_EOS
#ifdef LATTICE_QCD
	printf("equation_of_state_new::mdmde_quasi error: not eos switch yet\n");
	exit(-1);
#endif
#endif

	return mdmde;
}


precision equation_of_state_new::equilibrium_mean_field()
{
	precision Beq;

#ifdef LATTICE_QCD
#ifndef CONFORMAL_EOS
	Beq = (-1491.6371738990576 + 27163.388236990842*T1 - 212236.89063342367*T2 + 927466.5839512877*T3 - 2.4580524285908886e6*T4 + 3.926409286725362e6*T5 - 3.2755581489732354e6*T6 + 240347.64382754304*T7 + 1.909163778636297e6*T8 - 828276.4353817657*T9 - 1.0180016922236403e6*T10 + 724600.41780134*T11 + 523118.1964580224*T12 - 594543.5372984193*T13 - 115848.46275057436*T14 + 418946.66737534944*T15 - 264584.0981957474*T16 + 89772.07640009356*T17 - 23553.261674561967*T18 + 6097.776585093935*T19 - 930.1069258971695*T20)/(13390.31558372141 - 30456.102091940047*T1 + 1454.5335541405343*T2 + 22226.839099744822*T3 + 13350.101408847955*T4 - 8434.386989976701*T5 - 17563.162052045962*T6 - 8606.75840956273*T7 + 5591.979692842073*T8 + 11420.601270099178*T9 + 6034.222871353725*T10 - 3449.1797165455346*T11 - 7527.240010802675*T12 - 2233.637576196*T13 + 5558.831137569498*T14 + 2510.924373254798*T15 - 4890.926376475482*T16 + 1847.1590191373598*T17 - 205.7160783233205*T18 + 15.933973129547796*T19 - 0.5635677842887817*T20);
	
	// Beq = (-527.7198740151986 + 10083.738278232428*T1 - 83215.1096774566*T2 + 390612.31737945497*T3 - 1.156005341824903e6*T4 + 2.275888848621547e6*T5 - 3.061067641368322e6*T6 + 2.821552088783961e6*T7 - 1.734981576383687e6*T8 + 650753.3294530859*T9 - 102948.4969857233*T10 - 20516.9754475802*T11 + 12000.461040175463*T12 - 1654.3346738122705*T13 + 25.18968917586496*T14)/(-138191.2734866683 + 1.011801904770102e6*T1 - 3.3275105856632143e6*T2 + 6.497780831943228e6*T3 - 8.392061222991193e6*T4 + 7.557844214173838e6*T5 - 4.87965141702813e6*T6 + 2.2944968338930537e6*T7 - 797052.503381821*T8 + 208803.7925704919*T9 - 41919.097781420365*T10 + 6233.939818848559*T11 - 615.2405051998173*T12 + 36.26133736071459*T13 - 0.9678722816941575*T14);
#endif
#endif

#ifdef CONFORMAL_EOS
#ifndef LATTICE_QCD
	Beq = 0;
#endif
#endif

#ifdef CONFORMAL_EOS
#ifdef LATTICE_QCD
	printf("equation_of_state_new::equilibrium_mean_field error: not eos switch yet\n");
	exit(-1);
#endif
#endif

	return Beq;
}


precision equation_of_state_new::beta_shear()
{
	precision beta_shear;

#ifdef LATTICE_QCD
#ifndef CONFORMAL_EOS
	beta_shear = (-8.704487593184917 + 219.28918157517424*T1 - 2477.902566778392*T2 + 16530.42025253919*T3 - 71968.8684621151*T4 + 212357.11233567353*T5 - 420519.0026792284*T6 + 510017.5263064133*T7 - 227746.91973532928*T8 - 313315.51121954236*T9 + 470935.63873870415*T10 + 158066.72037613127*T11 - 938961.3539529926*T12 + 912207.4050632303*T13 - 162516.20652738435*T14 - 401315.26100589056*T15 + 385527.65045220935*T16 - 149881.22794496786*T17 + 22849.199812879146*T18)/(836.1605819652495 - 4994.2129627000995*T1 + 10339.445723645553*T2 - 7863.529406530712*T3 + 14517.067742422922*T4 - 85636.74419675485*T5 + 216425.54907636004*T6 - 286522.46068304294*T7 + 230319.1750154616*T8 - 155698.03604967098*T9 + 154838.0327505353*T10 - 162178.64531668008*T11 + 110057.7642019968*T12 - 41004.48870170546*T13 + 6678.652020583494*T14 - 123.27511354433277*T15 + 10.047607464120292*T16 - 0.511781325382738*T17 + 0.012109259631224398*T18);
#endif
#endif


#ifdef CONFORMAL_EOS
#ifndef LATTICE_QCD
   beta_shear = 4./15. * conformal_prefactor * T1 * T1 * T1 * T1;
#endif
#endif

#ifdef CONFORMAL_EOS
#ifdef LATTICE_QCD
	printf("equation_of_state_new::beta_shear error: not eos switch yet\n");
	exit(-1);
#endif
#endif

	return beta_shear;
}


precision equation_of_state_new::beta_bulk()
{
	precision beta_bulk;

#ifdef LATTICE_QCD
#ifndef CONFORMAL_EOS
	beta_bulk = (0.0011079189759896013 - 0.03607032819669769*T1 + 0.4414267814343914*T2 - 2.876909054432641*T3 + 11.416181759176936*T4 - 28.99095963617232*T5 + 46.949029047256914*T6 - 45.04881878936641*T7 + 19.652738082518503*T8 - 2.270789012391187*T9 + 18.02005251525725*T10 - 40.8215770738997*T11 + 31.80655015730408*T12 + 1.2003007553780203*T13 - 23.541763430860975*T14 + 23.51183245427725*T15 - 12.163227346297504*T16 + 0.5780864080028956*T17 + 5.7083888622600805*T18 - 5.768744942470897*T19 + 2.960799729426568*T20 - 0.8295762503448261*T21 + 0.10194139487600093*T22)/(10.9367675856047 - 97.67668689586313*T1 + 358.59297974309567*T2 - 648.0435677635525*T3 + 430.4494897227276*T4 + 415.9498232081677*T5 - 842.5138911823302*T6 + 39.44143106719552*T7 + 835.0988480486614*T8 - 398.0196479689704*T9 - 418.0163763613797*T10 + 74.33541264481028*T11 + 748.3740580290471*T12 - 502.54214599460346*T13 - 602.7442749814361*T14 + 1262.750834047705*T15 - 1058.9191239838513*T16 + 537.2906742929765*T17 - 179.45408838874548*T18 + 40.09993403446954*T19 - 5.904130059726229*T20 + 0.5345210582722899*T21 - 0.020839894199259086*T22);
#endif
#endif

#ifdef CONFORMAL_EOS
#ifndef LATTICE_QCD
   beta_bulk = 0;
#endif
#endif

#ifdef CONFORMAL_EOS
#ifdef LATTICE_QCD
	printf("equation_of_state_new::beta_bulk error: not eos switch yet\n");
	exit(-1);
#endif
#endif

   return beta_bulk;
}



























// it might be better to make an equation class to store the variables
// what would be better to use? a fit function or interpolation?

equation_of_state::equation_of_state(precision e_in)
{
	e1  = e_in;			// compute powers
	e2  = e1  * e1;
	e3  = e2  * e1;
	e4  = e3  * e1;
	e5  = e4  * e1;
	e6  = e5  * e1;
	e7  = e6  * e1;
	e8  = e7  * e1;
	e9  = e8  * e1;
	e10 = e9  * e1;
	e11 = e10 * e1;
	e12 = e11 * e1;
	e13 = e12 * e1;
	e14 = e13 * e1;
	e15 = e14 * e1;
	e16 = e15 * e1;
	e17 = e16 * e1;
	e18 = e17 * e1;
	e19 = e18 * e1;
	e20 = e19 * e1;
	e21 = e20 * e1;
	e22 = e21 * e1;
	e23 = e22 * e1;
	e24 = e23 * e1;
	e25 = e24 * e1;
	e26 = e25 * e1;
}

equation_of_state::~equation_of_state()
{

}


precision equation_of_state::equilibrium_pressure()
{
	precision p;

	// need to think about how to implement eos switching
	// for equation of state switching, should avoid macros (right this needs to be updated eventually)
#ifdef LATTICE_QCD
	p = (3.258203875251529e-11 - 7.747698746076113e-7*e1 - 0.0006981835859774962*e2 - 0.06174685910888912*e3 + 2.275922623600396*e4 - 37.44924095890363*e5 + 233.86414097675518*e6 + 19.464596980830084*e7 - 12874.726137720376*e8 + 64564.18316194397*e9 - 71526.35129303865*e10 + 33676.244329386056*e11 + 14038.245253970734*e12 - 20414.77445592913*e13 + 8726.519534553134*e14 - 81.50944485934848*e15 - 752.9942377524083*e16 + 258.28927588902246*e17 - 17.453939443529997*e18 + 11.355064193208065*e19 + 0.6357364785590592*e20 + 0.0091023656741957*e21 + 0.00003993269353049751*e22 + 4.8771388046509797e-8*e23 + 1.0405466824842022e-11*e24)/(-4.732328155036805e-6 - 0.00325496338546896*e1 - 0.2414677643015335*e2 + 8.446673595606319*e3 - 125.95148874838058*e4 + 478.4320052947472*e5 + 4620.6136514507825*e6 - 77854.87566345902*e7 + 310810.3126290067*e8 - 192579.3199639818*e9 - 159145.94873938535*e10 + 436458.4758128175*e11 - 351005.5467031545*e12 + 154372.28841797868*e13 - 31046.555167702114*e14 + 1842.3266104828965*e15 + 553.6716087575791*e16 + 65.80263807137636*e17 + 59.03967708029564*e18 + 2.5943899381998428*e19 + 0.032381183915126406*e20 + 0.00013101670868370324*e21 + 1.527281015828544e-7*e22 + 3.1747370586121694e-11*e23 + 1.0273281217533355e-17*e24);
#else
	p = e1 / 3.;
#endif
	if(p < 0)
	{
		printf("equation_of_state error: p = %lf is negative. Enforcing positive equilibrium_pressure\n", p);
		p = fmax(p, 0.);
	}
	return p;
}


precision equation_of_state::speed_of_sound_squared()
{
	precision cs2;

#ifdef LATTICE_QCD
	cs2 = (-3.098379231725401e-13 - 1.5057592378698152e-9*e1 + 1.8182477549678398e-7*e2 + 0.000265942360295654*e3 - 0.012169801366573908*e4 + 0.25838101370941824*e5 - 2.921437946175593*e6 + 16.963191123856383*e7 - 18.356276131836637*e8 - 293.3566651404204*e9 + 1095.410458138189*e10 - 874.716911001333*e11 + 3112.672324311266*e12 - 2020.806516817968*e13 + 460.99064883666415*e14 + 1220.3791075527747*e15 - 747.2397158495869*e16 + 244.50848087109242*e17 - 15.198338794559106*e18 + 7.3838001744904345*e19 + 0.35977103084379103*e20 + 0.0034073946427534752*e21 + 3.176466242678107e-6*e22)/(-2.029332754087875e-12 - 6.659314443251908e-9*e1 + 1.2786688789845978e-6*e2 + 0.0009960480360445558*e3 - 0.0408121389743489*e4 + 0.7038607492465706*e5 - 4.455539775096916*e6 - 25.243468414027458*e7 + 676.8145376717912*e8 - 4910.207748605095*e9 + 15483.360208014285*e10 - 26085.447785506487*e11 + 56744.64239156274*e12 - 57771.290373059644*e13 + 39984.54561622914*e14 - 9413.899973292644*e15 - 181.03406978253508*e16 + 868.729840087911*e17 - 25.779160873597263*e18 + 32.31113412515543*e19 + 1.2344162871494986*e20 + 0.01063222579594702*e21 + 9.683731676774376e-6*e22);
#else
	cs2 = 1./3.;
#endif
	return cs2;
}


precision equation_of_state::effective_temperature(precision conformal_prefactor)
{
	precision T;

#ifdef LATTICE_QCD
	T = (1.0097526475300114e-7 + 0.00021838483878330327*e1 + 0.0527665048468103*e2 + 1.3443038870124433*e3 - 29.734024004370692*e4 + 269.2769967747579*e5 + 1866.2652129862722*e6 - 3760.609426819512*e7 - 7339.004735460014*e8 + 1606.4002928709317*e9 + 22406.963604607237*e10 - 20502.863988381978*e11 - 1681.250215337671*e12 + 31774.024065918482*e13 - 8929.575720414807*e14 + 1060.5047767092692*e15 + 1556.1430619374776*e16 + 79.24641177767323*e17 + 0.8843663758438196*e18 + 0.002262761812575013*e19 + 9.097931466167723e-7*e20)/(6.52335668683017e-7 + 0.0008874335630072429*e1 + 0.1442403491043538*e2 + 1.6787993054371462*e3 - 56.63303674664881*e4 + 699.7663335101611*e5 + 1778.0902926314081*e6 - 4084.075030043697*e7 - 17409.764307553494*e8 + 24711.788909162286*e9 - 2882.4287321466554*e10 + 2951.387363303661*e11 - 20969.597421345836*e12 + 53914.5094221402*e13 - 19175.577404922642*e14 + 4057.4347774133366*e15 + 1550.173009591607*e16 + 54.95207917670983*e17 + 0.43418792515403754*e18 + 0.00074734801321898*e19 + 1.598170725740521e-7*e20);
#else
	T = pow(e1 / conformal_prefactor, 0.25);
#endif
	return T;
}


precision equation_of_state::z_quasi(precision T)
{
	precision z;	// mass / T

#ifdef LATTICE_QCD
	z = (-6.921627474395535e-11 - 2.068820219007291e-7*e1 - 0.00008004794094185476*e2 - 0.004170294272191056*e3 + 0.09994318330775476*e4 - 0.6028049712926495*e5 - 17.32997254829035*e6 + 221.86420008446356*e7 - 1161.528708956738*e8 - 1723.4122240576517*e9 + 339.98693004139375*e10 - 1371.3149826668873*e11 - 333.6959481389614*e12 - 321.321802085032*e13 + 131.94753133513424*e14 - 122.69087685089039*e15 + 5.370697224359077*e16 - 6.052394085388524*e17 - 0.7038207218913564*e18 - 0.01757160010668318*e19 - 0.00012991380242842205*e20 - 2.4594753226751026e-7*e21 - 7.352913098840251e-11*e22)/(-7.977262417393307e-12 - 2.786830987719674e-8*e1 - 0.000011762180354256496*e2 - 0.0006815255296076517*e3 + 0.014251257515341783*e4 - 0.0445565524634568*e5 - 3.3945388533360994*e6 + 37.08500361346862*e7 - 166.1413947461126*e8 - 471.2995122873681*e9 + 116.40828587862126*e10 - 380.4316969931931*e11 + 14.113450070587334*e12 - 186.43236298050886*e13 + 78.29603797399548*e14 - 42.553940433457285*e15 + 1.8975639224075795*e16 - 1.6361357363557714*e17 - 0.3074277311308281*e18 - 0.010250374977627853*e19 - 0.00009012790614650718*e20 - 1.9169861873399123e-7*e21 - 6.394439165661743e-11*e22);
#else
	z = 0;
#endif
	return z;
}


precision equation_of_state::mdmde_quasi()
{
	precision mdmde;

#ifdef LATTICE_QCD
	mdmde = (-1.1123631636297233e-16 + 1.4139243361369686e-14*e1 + 3.721480625854963e-11*e2 - 1.1051034554557073e-8*e3 + 4.293632053721226e-7*e4 + 0.00010102526719429564*e5 - 0.006711775629413175*e6 + 0.22211586478553663*e7 - 4.679885529595576*e8 + 68.90440687134294*e9 - 722.7034703795854*e10 + 5350.068342373652*e11 - 26821.524352996694*e12 + 85241.07873522438*e13 - 165070.78592788626*e14 + 220008.1433721365*e15 - 238821.97449171962*e16 + 192649.01193661397*e17 - 129339.10610817105*e18 + 47362.068861485284*e19 - 9098.976926100564*e20 + 345.1125470989994*e21 - 65.18059783017668*e22 + 0.5336316588688652*e23 + 0.1314236985827908*e24 + 0.0005071653091929672*e25 + 1.2622782105581163e-7*e26)/(-8.470594112097918e-20 - 3.0334858946453655e-16*e1 + 1.9221008294469902e-13*e2 + 3.2970489778770567e-12*e3 - 1.2829998132301256e-8*e4 + 1.3773757095343025e-6*e5 + 0.00002022154904931217*e6 - 0.002930770822779598*e7 + 0.10660208553788184*e8 - 2.2098135701741763*e9 + 30.171364252428436*e10 - 267.5251584133791*e11 + 1413.9660089847882*e12 - 2646.540333937091*e13 - 14290.774304779678*e14 + 87788.200508876*e15 - 167569.32061878*e16 + 240390.70234491615*e17 - 213428.71705193262*e18 + 138010.45687451633*e19 - 51209.10145588782*e20 + 11781.485918505261*e21 - 1217.3899418058113*e22 + 235.38359728669485*e23 + 12.01744346594077*e24 + 0.11795181077625658*e25 + 0.00015208597485394176*e26);
#else
	mdmde = 0;
#endif
	return mdmde;
}


precision equation_of_state::equilibrium_mean_field(precision T)
{
	precision Beq;

#ifdef LATTICE_QCD
	Beq = (-0.00014062752263119843 + 0.2379242517626019*e1 - 7.652070613493643*e2 + 103.74860438908097*e3 - 737.5821559333572*e4 + 3247.263538490574*e5 - 10515.549916941185*e6 + 11733.571827739737*e7 + 4891.077685224742*e8 - 34610.18746768503*e9 + 59089.19275974825*e10 - 51132.605796603144*e11 + 23240.95858982552*e12 - 6351.270469223261*e13 + 1243.3103601685511*e14 - 271.02350624683777*e15 + 77.28890422589751*e16 + 6.861570057022561*e17 - 0.061074992910151324*e18 - 0.0008134689379938509*e19 - 1.1256801223385317e-6*e20)/(-2.3158198765018394 + 68.74089731092687*e1 - 827.2937948361855*e2 + 4956.8913570637615*e3 - 19441.306038112383*e4 + 64613.252524115334*e5 + 18840.864260856848*e6 + 217.74430236461092*e7 - 19791.664416644526*e8 - 29180.14576093904*e9 - 39626.69721247169*e10 - 30613.476453009993*e11 + 62354.45307262448*e12 - 33917.77621528999*e13 + 9343.735868937449*e14 - 602.9029343012371*e15 + 399.2255336058777*e16 + 7.766954875118376*e17 + 0.036231945687140756*e18 + 0.00003768663903333004*e19 + 5.947311885856504e-11*e20);
#else
   Beq = 0;
#endif
   return Beq;
}


precision equation_of_state::beta_shear(precision T, precision conformal_prefactor)
{
	double beta_shear;
	double T1 = (double)T;

#ifdef LATTICE_QCD
	double T2  = T1  * T1;
	double T3  = T2  * T1;
	double T4  = T3  * T1;
	double T5  = T4  * T1;
	double T6  = T5  * T1;
	double T7  = T6  * T1;
	double T8  = T7  * T1;
	double T9  = T8  * T1;
	double T10 = T9  * T1;
	double T11 = T10 * T1;
	double T12 = T11 * T1;
	double T13 = T12 * T1;
	double T14 = T13 * T1;
	double T15 = T14 * T1;
	double T16 = T15 * T1;

	beta_shear = (0.11389523069479185 - 2.705127789711982*T + 28.742274634585783*T2 - 179.70799505135633*T3 + 730.506119388531*T4 - 2007.1240462918665*T5 +
	3721.0091345659325*T6 - 4366.923562236631*T7 + 2347.927819272824*T8 + 1559.8495170149938*T9 - 4257.643429592799*T10 +
	3938.5865394866464*T11 - 1999.1891431300426*T12 + 551.2544143662707*T13 - 64.67817053548931*T14)/
	(-8.065469391396581 + 92.36219194757648*T - 453.6119231622029*T2 + 1282.9299920278522*T3 - 2333.58921706985*T4 + 2843.5221412837436*T5 -
	2291.373282720429*T6 + 1092.688273629689*T7 - 140.80190149751965*T8 - 177.18993543735297*T9 + 126.47676060292035*T10 -
	39.29732109229348*T11 + 6.547795713678674*T12 - 0.6094551149533467*T13 + 0.023633386710670345*T14);
#else
	beta_shear = 4./15. * conformal_prefactor * T1 * T1 * T1 * T1;
#endif
	return beta_shear;
}


precision equation_of_state::beta_bulk(precision T)
{
	double beta_bulk;

#ifdef LATTICE_QCD
	double T1 = (double)T;
	double T2  = T1  * T1;
	double T3  = T2  * T1;
	double T4  = T3  * T1;
	double T5  = T4  * T1;
	double T6  = T5  * T1;
	double T7  = T6  * T1;
	double T8  = T7  * T1;
	double T9  = T8  * T1;
	double T10 = T9  * T1;
	double T11 = T10 * T1;
	double T12 = T11 * T1;

	beta_bulk = (-0.07227796340001295 + 1.3143764840615304*T - 10.618632711354438*T2 + 50.49226038033788*T3 - 157.4895501326704*T4 + 339.56597470170806*T5 -
	518.8479144079654*T6 + 565.401177446674*T7 - 434.90715638011426*T8 + 229.10104947248428*T9 - 77.73510380062434*T10 +
	14.99701288452569*T11 - 1.2012010330086933*T12)/
	(36.286491731381375 - 385.4575151988946*T + 1863.5710590765027*T2 - 5417.869823144005*T3 + 10538.668878148343*T4 - 14431.58178913021*T5 +
	14242.335632907661*T6 - 10182.425674573899*T7 + 5216.5453702051245*T8 - 1858.714375839709*T9 + 434.4504198478326*T10 -
	59.390426532205105*T11 + 3.581827839432567*T12);
#else
   beta_bulk = 0;
#endif
   return beta_bulk;
}


// this is only used in initial conditions (This is the BEST version)

precision equilibrium_energy_density(precision T, precision conformal_prefactor)
{
	precision e;

#ifdef LATTICE_QCD
	double T1 = (double)T;
	double T2 = T1 * T1;
	double T3 = T2 * T1;
	double T4 = T3 * T1;
	double T5 = T4 * T1;
	double T6 = T5 * T1;
	double T7 = T6 * T1;
	double T8 = T7 * T1;
	double T9 = T8 * T1;
	double T10 = T9 * T1;
	double T11 = T10 * T1;
	double T12 = T11 * T1;
	double T13 = T12 * T1;
	double T14 = T13 * T1;
	double T15 = T14 * T1;
	double T16 = T15 * T1;
	double T17 = T16 * T1;
	double T18 = T17 * T1;
	double T19 = T18 * T1;
	double T20 = T19 * T1;

	e = (0.40627518438643867 - 11.011572539859957*T1 + 130.77928985869593*T2 - 877.6121918925887*T3 + 3535.421433912612*T4 - 8078.914267299846*T5 + 6439.154404609683*T6 + 15129.165217897591*T7 - 38768.2682434212*T8 - 8551.62549456636*T9 + 126398.49102427797*T10 - 76579.50897685449*T11 - 311008.5908671992*T12 + 508668.90387515235*T13 + 355111.33690457407*T14 - 2.0632279590489047e6*T15 + 3.1401856418813644e6*T16 - 2.69069250116831e6*T17 + 1.421083090444586e6*T18 - 443129.7467853764*T19 + 64263.40447458266*T20)/(-59.6120292971754 + 717.1241260536837*T1 - 3645.3702819435307*T2 + 9765.42393031536*T3 - 12832.431386201157*T4 + 142.82590432252346*T5 + 25186.86259123896*T6 - 28316.004650857147*T7 - 17930.007137334782*T8 + 77200.73271762159*T9 - 80850.88417786805*T10 + 23824.933947902435*T11 + 25225.154910216275*T12 - 21367.56797902718*T13 - 7567.354921918593*T14 + 20258.528356323113*T15 - 13624.52545866848*T16 + 4563.423116419881*T17 - 761.8274089342926*T18 + 76.44496845688465*T19 - 3.4969057033074002*T20);
#else
	e = conformal_prefactor * T * T * T * T;
#endif
	return e;
}



// This is the HOTQCD + SMASH version
precision equilibrium_energy_density_new(precision T, precision conformal_prefactor)
{
	precision e;

#ifdef LATTICE_QCD
#ifndef CONFORMAL_EOS
	double T1 = T;
	double T2 = T1 * T1;
	double T3 = T2 * T1;
	double T4 = T3 * T1;
	double T5 = T4 * T1;
	double T6 = T5 * T1;
	double T7 = T6 * T1;
	double T8 = T7 * T1;
	double T9 = T8 * T1;
	double T10 = T9 * T1;
	double T11 = T10 * T1;
	double T12 = T11 * T1;
	double T13 = T12 * T1;
	double T14 = T13 * T1;
	double T15 = T14 * T1;
	double T16 = T15 * T1;
	double T17 = T16 * T1;
	double T18 = T17 * T1;

	e = (-0.40698677369016695 + 26.90263255626566*T1 - 451.40931096731606*T2 + 3748.50059581573*T3 - 18498.9525114122*T4 + 57417.46806755477*T5 - 109624.23666681767*T6 + 105527.65145392402*T7 + 24561.692186035198*T8 - 187413.3245788481*T9 + 143250.8494521164*T10 + 143393.6299794811*T11 - 337897.5513221823*T12 + 177408.70085616288*T13 + 137332.23093322886*T14 - 264268.15252307296*T15 + 175507.67395950254*T16 - 58067.60179151776*T17 + 8046.346263452795*T18)/(-173.8835104731062 + 1529.128399904588*T1 - 5768.283494742235*T2 + 11835.156422023792*T3 - 13416.39237101848*T4 + 6707.496209857805*T5 + 437.7325829205757*T6 + 704.6255257249378*T7 - 2851.15588379858*T8 - 5693.046198149581*T9 + 18686.5840958459*T10 - 21262.903672709494*T11 + 12920.1164479426*T12 - 4290.49508258817*T13 + 652.3067905178337*T14 - 18.705607055329114*T15 + 1.824854560118584*T16 - 0.10710494603163746*T17 + 0.002855313723557289*T18);
#endif
#endif

#ifdef CONFORMAL_EOS
#ifndef LATTICE_QCD
	e = conformal_prefactor * T * T * T * T;
#endif
#endif

#ifdef CONFORMAL_EOS
#ifdef LATTICE_QCD
	printf("equilibrium_energy_density error: not eos switch yet\n");
	exit(-1);
#endif
#endif

	return e;
}

