#include <math.h>
#include <cmath>
#include <stdlib.h>
#include "qcd.hpp"


double z_Quasiparticle(double T)
{
	// z = m/T (quasiparticle model; nf = 3 flavors)
  double T1 = T;
  double T2 = T * T;
  double T3 = T2 * T;
  double T4 = T3 * T;
  double T5 = T4 * T;
  double T6 = T5 * T;
  double T7 = T6 * T;
  double T8 = T7 * T;
  double T9 = T8 * T;
  double T10 = T9 * T;
  double T11 = T10 * T;
  double T12 = T11 * T;
  double T13 = T12 * T;
  double T14 = T13 * T;
  double T15 = T14 * T;
  double T16 = T15 * T;
  double T17 = T16 * T;
  double T18 = T17 * T;
  double T19 = T18 * T;
  double T20 = T19 * T;

  // z_quasi fit from hotqcd smash.nb

  // instructions for updating the return value:
  //    1) run the Mathematica notebook hotqcd smash.nb
  //    2) do sh cform.sh qcd_eos (only once!)
  //    3) copy/paste the C++ formula from thermal_function_fits/qcd_eos/z_quasi.txt (it should look like this one)

  return (0.5953313988318101 - 2.2043653777008614*T1 + 6.412298975786382*T2 - 55.999935626852256*T3 + 266.60638354021063*T4 - 613.7443665674955*T5 + 655.4361962409046*T6 - 15.985620813188662*T7 - 722.964445473731*T8 + 559.5675143434362*T9 + 331.3111327115995*T10 - 744.1136744233224*T11 + 232.70107465490355*T12 + 485.8050268195538*T13 - 772.6165293970012*T14 + 655.2619938242386*T15 - 392.69621238742855*T16 + 159.90084483026934*T17 - 34.580773905851075*T18 + 0.3552068133864607*T19 + 0.9529211080685058*T20)/(1.1621492604568943e-6 + 0.8316553006125841*T1 - 6.145185278692966*T2 + 18.593245529858304*T3 - 29.383537074661632*T4 + 29.562837167596356*T5 - 38.44431456235*T6 + 73.50473831542084*T7 - 91.7352784806228*T8 + 47.099875767546756*T9 + 17.377610722099814*T10 - 38.314006662117855*T11 + 19.446036400472437*T12 + 24.055731836948876*T13 - 80.16900358183341*T14 + 90.87090243010566*T15 - 30.23093829725649*T16 - 31.73728103092395*T17 + 38.419610895809136*T18 - 16.148353279684848*T19 + 2.5456531976175*T20);
}
