//Coupled model of photosynthesis-stomatal conductance-energy balance for a maize leaf
// this unit simulates Maize leaf gas-exchange characteristics
// including photosynthesis, traspiration, boundary and stomatal conductances,
// and leaf temperature based on von Caemmerer (2000) C4 model, BWB stomatal
// conductance (1987) and Energy balance model as described in Campbell and
// Norman (1998) photosynthetic parameters were calibrated with PI3733 from
// SPAR experiments at Beltsville, MD in 2002 Stomatal conductance parameters
// were not calibrated
// Ver 1.0, S.Kim, 11/2002, Originally written in Pascal
// Translated into C++ by S.Kim June 11, 2003 following D.Timlin's translation of C3 model
// IMPORTANT: This model must not be released until validated and published
// 2006-2007 modified to use leaf water potentials to adjust photosynthesis for water stress Y. Yang
// modified 2009 to adjust photosynthesis for nitroge stress Y. Yang

#include "stdafx.h"
#include "gas_exchange.h"
#include  <cmath>
#include <stdlib.h>

#define R 8.314  // idealgasconstant
#define maxiter 200
#define epsilon 0.97
#define sbc 5.6697e-8
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
inline double Square(double a) { return a * a; }
inline double Min(double a, double b, double c) {return (__min(__min(a,b),c));}



CGas_exchange::CGas_exchange(std::string sType_in, double n_content)
{
	sType=sType_in;
	lfNContent = n_content; 
}
CGas_exchange::~CGas_exchange()
{
}

void CGas_exchange::getParms()
{
	Parms.EaVp     =        75100;
	Parms.EaVc     =        55900; // Sage (2002) JXB
	Parms.Eaj      =        32800;
	Parms.Hj       =       220000;
	Parms.Sj       =          702.6;
	Parms.Vpm25    =          70.0;
	Parms.Vcm25    =          50.0;
	Parms.Jm25     =         300.0;
	// Kim et al. (2007), Kim et al. (2006)
	//In von Cammerer (2000), Vpm25=120, Vcm25=60,Jm25=400
	//In Soo et al.(2006), under elevated C5O2, Vpm25=91.9, Vcm25=71.6, Jm25=354.2 YY
	// Values in Kim (2006) are for 31C, and the values here are normalized for 25C. SK
	Parms.Rd25     =          2.0;
	Parms.Ear      =        39800;
	Parms.g0 = 0.04;
    Parms.g1 =  4.0;   //in P. J. Sellers, et al.Science 275, 502 (1997), g0 is b, of which the value for c4 plant is 0.04
	                   //and g1 is m, of which the value for c4 plant is about 4 YY)
	Parms.beta_ABA = 1.48e2; //Tardieu-Davies beta, Dewar (2002) Need the references !?
	Parms.delta = -1.0;
	Parms.a_ABA = 1.0e-4;
	Parms.lamda_r = 4.0e-12; // Dewar's email
	Parms.lamda_l = 1.0e-12; 
	Parms.K_max = 6.67e-3; //max. xylem conductance (mol m-2 s-1 MPa-1) from root to leaf, Dewar (2002);
}

	



void CGas_exchange::GasEx_psil(double leafp, double et_supply)
{ //runs
    double Tleaf_old;
    int   iter=1;

    Tleaf = Tair; Tleaf_old = 0;
	Ci = 0.4*CO2;
	gb = gbw();
	gs = gsw(leafp);
    A_net = (CO2-Ci)/(1.57/gs+1.37/gb)*Press/100;
	
	while ((fabs(Tleaf_old -Tleaf)>0.01) && (iter < maxiter))
    {
      Tleaf_old = Tleaf;
      Photosynthesis(leafp);
      EnergyBalance(et_supply);
	  iter2 =++iter; //iter=iter+1, iter2=iter; 
    } 
}






void CGas_exchange::Photosynthesis(double pressure)    //Incident PFD, Air temp in C, CO2 in ppm, RH in percent
{ //keep
    const double    f = 0.15;             //spectral correction
    const double    O = 210;             // gas units are mbar
    const double    theta = 0.5;
    const double    scatt = 0.15;        //leaf reflectance + transmittance
    const int       Kc25 = 650;    //* Michaelis constant of rubisco for CO2 of C4 plants (2.5 times that of tobacco), ubar, Von Caemmerer 2000 */
    const int       Ko25 = 450;    //* Michaelis constant of rubisco for O2 (2.5 times C3), mbar */
	const int		Kp25 = 80;     //* Michaelis constant for PEP caboxylase for CO2 */
    const long      Eac = 59400;  const long       Eao = 36000;   // activation energy values
	const int	    Vpr = 80; //* PEP regeneration limited Vp, value adopted from vC book */
    const double	gbs = 0.003; //* bundle sheath conductance to CO2, mol m-2 s-1 */
    const double    x = 0.4;  //* Partitioning factor of J, yield maximal J at this value */
    const double    alpha = 0.0001; //* fraction of PSII activity in the bundle sheath cell, very low for NADP-ME types  */
    const double    gi = 1.0; //* conductance to CO2 from intercelluar to mesophyle, mol m-2 s-1, assumed */
    const double    beta = 0.99; //* smoothing factor */
    const double    gamma1 = 0.193; //* half the reciprocal of rubisco specificity, to account for O2 dependence of CO2 comp point, note that this become the same as that in C3 model when multiplied by [O2] */

	double Kp, Kc, Ko, Km, Ia, I2, Tk, Vpmax, Jmax, Vcmax,  Om, Rm, J, Ac1, Ac2, Ac, Aj1,
		   Aj2, Aj, Vp1, Vp2, Vp, P,  Ca, Cm, Cm_last, Cm_next,
		    Os, GammaStar, Gamma;
    int iter;
   getParms(); //Reset values changed for N status
//* Light response function parameters */
  Ia = PFD*(1-scatt);    //* absorbed irradiance */
  I2 = Ia*(1-f)/2;    //* useful light absorbed by PSII */

//* other input parameters and constants */
  Tk = Tleaf + 273.0;
  P  = Press/100;
  Ca = CO2*P; //* conversion to partial pressure */
  Om = O;   //* mesophyle O2 partial pressure */
  Kp = Kp25; // T dependence yet to be determined
  Kc = Kc25*exp(Eac*(Tleaf-25)/(298*R*(Tleaf+273)));
  Ko = Ko25*exp(Eao*(Tleaf-25)/(298*R*(Tleaf+273)));
  Km = Kc*(1+Om/Ko); //* effective M-M constant for Kc in the presence of O2 */
  Rd = Parms.Rd25*exp(Parms.Ear*(Tleaf-25)/(298*R*(Tleaf+273)));

  double CriticalNitrogen;
  CriticalNitrogen=__max(0.25,lfNContent);
  Parms.Vcm25 = Parms.Vcm25*(2/(1+exp(-2.9*(CriticalNitrogen-0.25)))-1);
  Parms.Jm25 = Parms.Jm25*(2/  (1+exp(-2.9*(CriticalNitrogen-0.25)))-1);
  Parms.Vpm25 = Parms.Vpm25*(2/(1+exp(-2.9*(CriticalNitrogen-0.25)))-1); //in Sinclair and Horie, 1989 Crop sciences, it is 4 and 0.2; 
                                                                  //In J Vos. et al. Field Crop study, 2005, it is 2.9 and 0.25;
                                                                  //In Lindquist, weed science, 2001, it is 3.689 and 0.5. 
  
  Vpmax = Parms.Vpm25*exp(Parms.EaVp*(Tleaf-25)/(298*R*(Tleaf+273)));
  Vcmax = Parms.Vcm25*exp(Parms.EaVc*(Tleaf-25)/(298*R*(Tleaf+273)));
  Jmax = Parms.Jm25*exp(((Tk-298)*Parms.Eaj)/(R*Tk*298))*(1+exp((Parms.Sj*298-Parms.Hj)/(R*298)))/(1+exp((Parms.Sj*Tk-Parms.Hj)/(R*Tk)));
  Rm = 0.5*Rd;

  Cm = Ci; //* mesophyle CO2 partial pressure, ubar, one may use the same value as Ci assuming infinite mesohpyle conductance */
  Cm_next = Cm;
  iter = 1;
  Cm_last = 0;
//* iteration to obtain Cm from Ci and A, could be re-written using more efficient method like newton-raphson method */
  while ((fabs(Cm-Cm_last) > 0.01) && (iter < maxiter))
  {
      Cm_last = Cm;
      gs = gsw(pressure);
      Cm = __min(__max(0.0, Ca - A_net*(1.6/gs + 1.37/gb)*P), 2*Ca);
      Vp1 = (Cm*Vpmax)/(Cm+Kp); //* PEP carboxylation rate, that is the rate of C4 acid generation */
      Vp2 = Vpr;
      Vp = __max(__min(Vp1, Vp2),0);
      //* Enzyme limited A (Rubisco or PEP carboxylation */
      Ac1 = (Vp+gbs*Cm-Rm); 
	 // Ac1 = __max(0,(Vp+gbs*Cm-Rm)); //prevent Ac1 from being negative Yang 9/26/06
      Ac2 = (Vcmax-Rd);
      Ac = __min(Ac1,Ac2);
      //* Light and electron transport limited  A mediated by J */
      J =  QuadSolnLower(theta, -(I2 + Jmax), I2*(Jmax)) ; //* rate of electron transport */
      Aj1 = (x*J/2-Rm+gbs*Cm);
      Aj2 = (1-x)*J/3-Rd;
      Aj = __min(Aj1,Aj2);
      A_net = ((Ac+Aj) - sqrt(Square(Ac+Aj)-4*beta*Ac*Aj))/(2*beta); //* smooting the transition between Ac and Aj */
      iter++;
  }

//  Convergence = True;
  iter1 = iter;
  Os = alpha*A_net/(0.047*gbs)+Om; //* Bundle sheath O2 partial pressure, mbar */
//       Cbs = Cm + (Vp-A_net-Rm)/gbs; //* Bundle sheath CO2 partial pressure, ubar */
  GammaStar = gamma1*Os;
  Gamma = (Rd*Km + Vcmax*GammaStar)/(Vcmax-Rd);
  Ci = Cm;
  A_gross = __max(0, A_net + Rd); // gets negative when PFD = 0, Rd needs to be examined, 10/25/04, SK
  
}



void CGas_exchange::EnergyBalance(double Jw)
// see Campbell and Norman (1998) pp 224-225
// because Stefan-Boltzman constant is for unit surface area by denifition,
// all terms including sbc are multilplied by 2 (i.e., gr, thermal radiation)
{ 
    const long lamda = 44000;  //KJ mole-1 at 25oC
    const double psc = 6.66e-4;
    const double Cp = 29.3; // thermodynamic psychrometer constant and specific hear of air (J mol-1 °C-1)
	double gha, gv, gr, ghr, psc1, Ea, thermal_air, Ti, Ta;
    Ta = Tair;
    Ti = Tleaf;
    gha = gb*(0.135/0.147);  // heat conductance, gha = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Mol m-2 s-1 ?
    gv = gs*gb/(gs+gb); 
    gr = (4*epsilon*sbc*pow(273+Ta,3)/Cp)*2; // radiative conductance, 2 account for both sides
    ghr = gha + gr;
    thermal_air = epsilon*sbc*pow(Ta+273,4)*2; // emitted thermal radiation
    psc1 = psc*ghr/gv; // apparent psychrometer constant
    this->VPD = Es(Ta)*(1-RH); // vapor pressure deficit
    Ea = Es(Ta)*RH; // ambient vapor pressure
	// debug dt I commented out the changes that yang made for leaf temperature for a test. I don't think they work
	if (Jw==0)
	{
	  Tleaf = Ta + (psc1/(Slope(Ta) + psc1))*((R_abs-thermal_air)/(ghr*Cp)-this->VPD/(psc1*Press)); //eqn 14.6b linearized form using first order approximation of Taylor series
	}
	else
	{
		Tleaf = Ta + (R_abs-thermal_air-lamda*Jw)/(Cp*ghr);
	}
    double Es_leaf = Es(Tleaf);
	double temp = Slope(Ta);
	double temp1 = R_abs-thermal_air;
    ET = __max(0, gv*((Es(Tleaf)-Ea)/Press)/(1-(Es(Tleaf)+Ea)/(Press))); //04/27/2011 dt took out the 1000 everything is moles now
}




double CGas_exchange::gsw(double pressure)  // stomatal conductance for water vapor in mol m-2 s-1 
{ //keep
	double Ds, aa, bb, cc, hs, Cs, Gamma, tmp;
		
	double temp = set_leafpEffect(pressure);
	Gamma = 10.0;
    Cs = CO2 - (1.37*A_net/gb); // surface CO2 in mole fraction
	    if (Cs <= Gamma) Cs = Gamma + 1;
        aa = temp*Parms.g1*A_net/Cs;
        bb = Parms.g0+gb-(temp*Parms.g1*A_net/Cs);
        cc = (-RH*gb)-Parms.g0;
        hs = QuadSolnUpper(aa,bb,cc);
	if (hs > 1) hs = 1; else if (hs <= 0.30) hs = 0.30;    //preventing bifurcation
    Ds = (1-hs)*Es(Tleaf); // VPD at leaf surface 
    tmp = (Parms.g0+Parms.g1*temp*(A_net*hs/Cs));

	// this below is an example of how you can write temporary data to a debug window. It can be copied and 
	// pasted into excel for plotting. Dennis See above where the CString object is created.
	//CString strdbg;
	//strdbg.Format("Source = %s tmp = %f pressure = %f Ds= %f Tleaf = %f Cs = %f Anet = %f hs = %f RH = %f\n",
		// sType, tmp, pressure, Ds, Tleaf, Cs, A_net, hs, RH);
    //OutputDebugString(strdbg);	

	if (tmp < Parms.g0) tmp=Parms.g0;
	return tmp;
}


//void CGas_exchange::set_leafpEffect(double leafp)  //calculating effect of leaf water potential on stomatal conductance
                                                   //model of Tuzet et al. 2003 used  Yang 8/21/06


double CGas_exchange::set_leafpEffect(double pressure)
{ 
  //pressure - leaf water potential MPa...
  double sf, phyf; 
  sf = 2.3;  //sensitivity parameter Tuzet et al. 2003 Yang
  phyf = -1.2; //reference potential Tuzet et al. 2003 Yang
  this->leafpEffect = (1+exp(sf*phyf))/(1+exp(sf*(phyf-pressure)));
  double temp = this->leafpEffect;
  return temp;
}



double CGas_exchange::gbw(void)
{
	const double stomaRatio = 1.0; // maize is an amphistomatous species, assume 1:1 (adaxial:abaxial) ratio.
	double ratio;
	double d;
	ratio = Square(stomaRatio+1)/(Square(stomaRatio)+1);
    d = width*0.72; // characteristic dimension of a leaf, leaf width in m
  //  return 1.42; // total BLC (both sides) for LI6400 leaf chamber
    return (1.4*0.147*sqrt(__max(0.1,wind)/d))*ratio; 
    // return (1.4*1.1*6.62*sqrt(wind/d)*(Press/(R*(273.15+Tair)))); // this is an alternative form including a multiplier for conversion from mm s-1 to mol m-2 s-1
	// 1.1 is the factor to convert from heat conductance to water vapor conductance, an avarage between still air and laminar flow (see Table 3.2, HG Jones 2014)
	// 6.62 is for laminar forced convection of air over flat plates on projected area basis
	// when all conversion is done for each surface it becomes close to 0.147 as given in Norman and Campbell
	// multiply by 1.4 for outdoor condition, Campbell and Norman (1998), p109, also see Jones 2014, pg 59 which suggest using 1.5 as this factor.
	// multiply by ratio to get the effective blc (per projected area basis), licor 6400 manual p 1-9
}

double CGas_exchange::Es (double T) //Campbell and Norman (1998), p 41 Saturation vapor pressure in kPa
{
    return (0.611*exp(17.502*T/(240.97+T)));
}

double CGas_exchange::Slope(double T) // slope of the sat vapor pressure curve: first order derivative of Es with respect to T
{
   double Temp, Temp1, Temp2;
 // units of b and c are  degrees C
   const double b= 17.502; const double c= 240.97;
   Temp1=Es(T);
   Temp2=Square(c+T);
   Temp=(Es(T)*(b*c)/Square(c+T)/Press); // for checking the value
   return (Es(T)*(b*c)/Square(c+T)/Press); 
}


void CGas_exchange::SetVal_psil(double PFD, double Tair, double CO2, double RH, double wind,  double Press, double width, double leafp, double ET_supply)
{
	const  double scatt = 0.15;
	this->PFD = PFD;
    double PAR = (PFD/4.55);
    double NIR = PAR; // If total solar radiation unavailable, assume NIR the same energy as PAR waveband
    this->R_abs = (1-scatt)*PAR + 0.15*NIR + 2*(epsilon*sbc*pow(Tair+273,4)); // times 2 for projected area basis
	// shortwave radiation (PAR (=0.85) + NIR (=0.15) solar radiation absorptivity of leaves: =~ 0.5
    this->CO2 = CO2;
    this->RH = __min(100.0, __max(RH, 10.0))/100;
    this->Tair = Tair;  //C
    this->width = width/100;  //meters
	this->wind = wind;        //meters s-1
    this->Press = Press;      //kPa
	this->leafp = leafp;
	this->leafpEffect = 1; //At first assume there is not drought stress, so assign 1 to leafpEffect. Yang 8/20/06
    getParms();
    GasEx_psil(leafp, ET_supply); //override GasEx() function so as to pass leaf water potential
	                    //into the calculation of gas exchange. Yang 8/20/06
}




double QuadSolnUpper (double a, double b, double c )
{
    if (a==0) return 0;
	else if ((b*b - 4*a*c) < 0) return -b/a;   //imaginary roots
    else  return (-b+sqrt(b*b-4*a*c))/(2*a);
}

double QuadSolnLower (double a, double b, double c )
{
    if (a==0) return 0;
	else if ((b*b - 4*a*c) < 0) return -b/a;   //imaginary roots
    else  return (-b-sqrt(b*b-4*a*c))/(2*a);
}
