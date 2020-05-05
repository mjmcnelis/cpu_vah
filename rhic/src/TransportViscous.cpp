#include <stdlib.h>
#include "../include/Macros.h"
#include "../include/EquationOfState.h"
#include "../include/TransportViscous.h"
#include "../include/Precision.h"
#include "../include/Parameters.h"

viscous_transport_coefficients::viscous_transport_coefficients(double T_in, double e_in, double p_in, int kinetic_theory_model_in)
{
	T = T_in;
	e = e_in;
	p = p_in;
	kinetic_theory_model = kinetic_theory_model_in;

	T1  = T;
	T2  = T1  * T;
	T3  = T2  * T;
	T4  = T3  * T;
	T5  = T4  * T;
	T6  = T5  * T;
	T7  = T6  * T;
	T8  = T7  * T;
	T9  = T8  * T;
	T10 = T9  * T;
	T11 = T10 * T;
	T12 = T11 * T;
	T13 = T12 * T;
	T14 = T13 * T;
	T15 = T14 * T;
	T16 = T15 * T;
	T17 = T16 * T;
	T18 = T17 * T;
	T19 = T18 * T;
	T20 = T19 * T;
	T21 = T20 * T;
	T22 = T21 * T;
}


viscous_transport_coefficients::~viscous_transport_coefficients()
{

}


void viscous_transport_coefficients::compute_shear_transport_coefficients(precision etas)
{
	if(kinetic_theory_model == 0)				// small fixed m/T << 1
	{
	#ifdef PIMUNU
		taupi_inverse = T / (5. * etas);
		betapi = (e + p) / 5.;
		delta_pipi = 4./3.;
		tau_pipi = 10./7;
	#ifdef PI
		lambda_pibulkPi = 1.2;
	#endif
	#endif
	}
	else										// quasiparticle m/T
	{
	#ifdef PIMUNU
		betapi = (-27.142377054641432 + 462.7152568336012*T1 - 3088.911493836734*T2 + 9294.165656439045*T3 - 5546.362522785967*T4 - 43310.703454899645*T5 + 122769.7671032288*T6 - 103467.21649054765*T7 - 54811.25840471966*T8 + 112983.19635536015*T9 + 33003.2645574495*T10 - 86513.1184180193*T11 - 31137.48855008747*T12 + 50820.04595819615*T13 + 21301.57924579058*T14 - 25062.560360078023*T15 - 3501.9015715618584*T16 + 18229.12958093956*T17 - 16217.694005589445*T18 - 3457.235338268151*T19 + 12393.178531714608*T20 - 6597.871611845416*T21 + 1470.5575446982516*T22)/(1015.5789835794669 - 7321.710661943294*T1 + 16400.56439744739*T2 - 11859.495858498654*T3 - 4578.450202173476*T4 + 6669.392709012055*T5 + 4352.0336723167575*T6 - 2746.3644144823443*T7 - 4767.192568645696*T8 - 975.4366636513871*T9 + 3400.4703589080195*T10 + 3500.6952093601494*T11 - 463.088334916719*T12 - 3332.053030372763*T13 - 1131.0706154419624*T14 + 2028.5128862558045*T15 - 27.60401477286506*T16 - 234.79138256755024*T17 + 3.310604161225659*T18 + 65.90446768041976*T19 - 7.579889386648616*T20 + 0.5007173440741145*T21 - 0.014536235575445516*T22);

		taupi_inverse = T * betapi / ((e + p) * etas);

		delta_pipi = (290.1364770212991 - 2344.2374896737947*T1 + 7705.597405115362*T2 - 12175.058183968922*T3 + 7594.429413531012*T4 - 1763.325991489814*T5 + 14041.482894863113*T6 - 23753.270748274073*T7 - 10363.183651269264*T8 + 40485.272396799504*T9 + 7940.160522257949*T10 - 51551.101906932025*T11 - 9743.077085423502*T12 + 55133.33527602634*T13 + 15502.401708250372*T14 - 51552.20036436098*T15 - 25018.003732821147*T16 + 50232.89849131191*T17 + 40872.26771774201*T18 - 107998.10808613256*T19 + 80447.539511616*T20 - 27840.839471109335*T21 + 3856.8849054519865*T22)/(177.93397574571506 - 1236.040125960039*T1 + 2647.9157253558237*T2 + 1898.1020657025902*T3 - 18880.1467909692*T4 + 32835.631636424376*T5 - 13986.623023051017*T6 - 21369.770041580945*T7 + 15797.10677816289*T8 + 17434.314395717025*T9 - 7335.7302418985155*T10 - 18898.168567709392*T11 - 6804.814642776499*T12 + 22752.307517120906*T13 + 23816.825903146684*T14 - 31004.655800962577*T15 - 37090.379952027004*T16 + 52225.16331977536*T17 + 24335.305742409895*T18 - 79528.37254625904*T19 + 60205.230161727224*T20 - 20882.24058476642*T21 + 2891.10510273954*T22);

		tau_pipi = (145.35316880200702 + 945.8325735330569*T1 - 9432.059748907666*T2 + 22596.42779440415*T3 - 9784.043776272396*T4 - 34181.96035568426*T5 + 34197.702700924594*T6 + 30055.934244274533*T7 - 36561.874733744364*T8 - 28335.292218000726*T9 + 5977.11973697689*T10 + 71234.42447385502*T11 - 19709.508486414004*T12 - 82629.40191181841*T13 + 68245.16973379398*T14 - 25418.083938156746*T15 + 45281.46732287938*T16 - 38396.78562887629*T17 - 15722.536401858844*T18 + 35065.044789858*T19 - 14703.327275131714*T20 + 495.3046612376971*T21 + 635.0944225725251*T22)/(69.9794464389132 + 772.8645690587082*T1 - 6817.641409568438*T2 + 17584.264254853453*T3 - 15010.797317231938*T4 - 7604.751289966955*T5 + 12135.331636783556*T6 + 13007.192278996585*T7 - 6668.68065563448*T8 - 29953.01393477996*T9 + 11951.131959755548*T10 + 31191.392088457585*T11 - 8055.18103666383*T12 - 15957.418747977383*T13 - 31662.729960808883*T14 + 56472.190299920985*T15 - 12650.878253326782*T16 - 7907.762379445128*T17 - 17210.85822297314*T18 + 25963.67908353706*T19 - 10431.250099856581*T20 + 339.8973007777852*T21 + 443.0411392095061*T22);

	#ifdef PI
		lambda_pibulkPi = (9.879822491210334 + 0.8683013755063738*T1 - 412.0168092840492*T2 + 1859.8897057941656*T3 - 3386.9698076878626*T4 + 1989.5209222287162*T5 + 1196.9683409630013*T6 + 1843.4532383308515*T7 - 11669.90892798136*T8 + 13589.591202560308*T9 - 2064.5832531605697*T10 - 5129.924321205942*T11 - 1270.9597932639356*T12 + 4969.691189124399*T13 + 1476.4287689692435*T14 - 2469.6798775791253*T15 - 7366.048915431875*T16 + 13425.678258951335*T17 - 8867.023678280702*T18 + 2136.9070574729917*T19 + 425.1877672917072*T20 - 341.04228416885667*T21 + 54.09326254145119*T22)/(4.937861274931255 + 14.935238446135134*T1 - 343.85758817237814*T2 + 1504.4350271262738*T3 - 3025.9542913799596*T4 + 2671.1134258080274*T5 + 308.62701814901914*T6 - 2464.0300399859802*T7 + 1864.3341840548514*T8 - 1981.1219328572442*T9 + 2504.7735656955665*T10 + 2568.5435792825647*T11 - 9337.769790966759*T12 + 5653.562525933961*T13 + 5005.865164685112*T14 - 6141.7134596114465*T15 - 3401.1197153364037*T16 + 9251.848595994144*T17 - 6124.4610901660335*T18 + 1245.259122596452*T19 + 465.6533817587116*T20 - 287.13463604549423*T21 + 43.27396774217259*T22);
	#endif
	#endif
	}
}


void viscous_transport_coefficients::compute_bulk_transport_coefficients(precision zetas, precision third_cs2)
{
	if(kinetic_theory_model == 0)				// small fixed m/T << 1
	{
	#ifdef PI
		taubulk_inverse	= 15. * third_cs2 * third_cs2 * T / zetas;		// double check this
		betabulk = 15. * third_cs2 * third_cs2 * (e + p);
		delta_bulkPibulkPi = 2./3.;
	#ifdef PIMUNU
		lambda_bulkPipi = 1.6 * third_cs2;
	#endif
	#endif
	}
	else										// quasiparticle m/T
	{
	#ifdef PI
		betabulk = (-0.007759718657876827 + 0.13527701183460097*T1 - 1.007325721445928*T2 + 4.094006434717201*T3 - 9.08657965249736*T4 + 6.256305009441373*T5 + 24.00870073056565*T6 - 82.77151389599645*T7 + 114.99885289024945*T8 - 41.30662650662664*T9 - 123.40789469814824*T10 + 218.9629399579876*T11 - 107.35118678049548*T12 - 110.3502510719402*T13 + 177.84317606715905*T14 - 14.291071345426523*T15 - 189.09921836490983*T16 + 240.92678202232733*T17 - 158.11181677331444*T18 + 62.4799308503998*T19 - 14.596624269558063*T20 + 1.7479763171335152*T21 - 0.06607848782788789*T22)/(14.2447584800057 - 143.0380383405217*T1 + 614.8883332122975*T2 - 1425.7281326691227*T3 + 1753.3236957895263*T4 - 619.9402701118898*T5 - 1166.5869323039597*T6 + 1355.0856774027493*T7 + 234.5120173060178*T8 - 1000.5813430192038*T9 - 164.88968228295698*T10 + 1084.016783336133*T11 - 163.14576407262282*T12 - 1280.9908387234796*T13 + 1608.0736134282931*T14 - 997.9363437538555*T15 + 383.5274234695*T16 - 107.08189506761974*T17 + 28.683056883497215*T18 - 7.769351756848533*T19 + 1.4593392411101056*T20 - 0.1306427916323396*T21 + 0.004536368033247522*T22);

		taubulk_inverse = T * betabulk / ((e + p) * zetas);

		delta_bulkPibulkPi = (-61.61231938958508 + 806.6821677513204*T1 - 4673.811226156834*T2 + 15236.287917055077*T3 - 29027.408534326747*T4 + 28263.742159301823*T5 - 1567.3564057188726*T6 - 22637.647081077546*T7 + 6597.76822674973*T8 + 23681.78826888058*T9 - 12138.484557927837*T10 - 19833.331819860814*T11 + 10973.194698033922*T12 + 17830.31559833379*T13 - 8926.119222102436*T14 - 16918.605763029642*T15 + 8959.17806057012*T16 + 19420.334383673337*T17 - 28189.958971142518*T18 + 16291.389808161735*T19 - 4583.2818996397455*T20 + 480.40397634992524*T21 + 16.53253610602593*T22)/(-27.477910316962866 + 180.06161425136645*T1 - 50.51956245444595*T2 - 3233.13204255847*T3 + 14440.761916675727*T4 - 29419.572729042873*T5 + 25962.73975810756*T6 + 8661.855973989968*T7 - 35065.432174092835*T8 + 3713.272572182705*T9 + 43716.26701628178*T10 - 18901.359998391614*T11 - 44956.13471353282*T12 + 29113.404048071432*T13 + 50609.25347438864*T14 - 59284.29889354307*T15 - 23738.396172313736*T16 + 85806.9132772426*T17 - 70800.84239686087*T18 + 27611.93279001703*T19 - 3948.5694518857335*T20 - 588.9749588157591*T21 + 198.24856308861118*T22);

	#ifdef PIMUNU
		lambda_bulkPipi = (165.9169457564851 - 1733.7255101736666*T1 + 8141.661368121232*T2 - 21552.520797330573*T3 + 32165.478060847825*T4 - 19604.407771507897*T5 - 13917.481148794082*T6 + 28284.235009617634*T7 - 4542.615058278952*T8 - 16626.136171010978*T9 + 8458.702146461757*T10 + 757.0687966250806*T11 - 1876.2933220036368*T12 + 18027.46983725655*T13 - 37290.035630808765*T14 + 33574.7686004331*T15 - 20505.02545193438*T16 + 17628.8331647611*T17 - 17018.27029304184*T18 + 10285.225557917576*T19 - 3205.077703340601*T20 + 365.5856994160321*T21 + 16.643673120140875*T22)/(298.763309725077 - 2116.2916420222905*T1 + 6740.778573027383*T2 - 16934.805964025552*T3 + 47858.83298949663*T4 - 112882.27822256282*T5 + 152535.73896841123*T6 - 55598.9929108521*T7 - 117031.33258443855*T8 + 117562.94685582837*T9 + 44112.68620677705*T10 - 38706.804729682306*T11 - 120403.40655041352*T12 + 55456.081432853956*T13 + 157806.49780950157*T14 - 120196.62423650328*T15 - 25098.931207320944*T16 - 115402.55172559206*T17 + 389664.73550140776*T18 - 420024.12948954443*T19 + 230958.18008319841*T20 - 66732.00430546697*T21 + 8132.911849886708*T22);
	#endif
	#endif
	}
}


