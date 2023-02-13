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

// 2006-2007 modified to use leaf water potentials to adjust photosynthesis for water stress Y. Yang
// modified 2009 to adjust photosynthesis for nitroge stress Y. Yang

// there are two ways to do max to make it linux friendly __max(blah, blah) or (std::max)(blah,blah)
// the second form has to be used with std as there is some macro somewhere for max

/*! @file
   Contains code to simulate gas exchange 
   */
/*

class CGasExchange 
   Additional functions and variables
   file GasExchange.CPP Holds class definitions
*/

#include "stdafx.h"
#include "gas_exchange.h"
#include <cmath>
#include <stdlib.h>

	// General fixed parameters
#define R 8.314  //!< \b R idealgasconstant
#define maxiter 200 //!< \b maxiter maximum number of iterations
#define epsilon 0.97   //!< epsilon emissivity See Campbell and Norman, 1998, page 163 (CHECK) 
#define sbc 5.6697e-8  //!< stefan-Boltzmann constant Wm-2 k-4. Actually varies somewhat with temperature
#define    O 205.0     //!<  Oxygen partial pressure gas units are mbar
#define Q10 2.0        //!< Q10 factor



#ifdef _DEBUG
#define new DEBUG_NEW
#endif

	
	//note that '<' indicates member is before the block and not after for Doxygen

	CGasExchange::CGasExchange(std::string sType_in, double n_content, const TGasExSpeciesParam& photoparam)
		//! This is the constructor

		/** Constructor - initialization of some variables is done here.

		*/ 
	{ 
		sType = sType_in;
		lfNContent = n_content;

		isCiConverged=false;
		errTolerance = 0.001;
		eqlTolerance = 1.0e-6;
		sParms = photoparam;
		

	}


	CGasExchange::~CGasExchange()
		//! This is the destructor
		/**Destructor - nothing is done here
		*/
	{
	}

	void CGasExchange::SetVal(double PhotoFluxDensity, double Tair, double CO2, double RH, double wind,  
		                      double Press, 
		                      double width, double leafp, double ET_supply)
		/**Sets environment variables for a single execution of the module 

		* Calls GasEx() to calculate photosynthetic rate and stomatal conductance.

		* @param[in] PhotoFluxDensity	Photosynthetic Flux Density (umol Quanta m-2 s-1) (check)
		* @param[in] Tair	Air Temperature (C)
		* @param[in] CO2	CO2 concentration of the air (umol mol-1)
		* @param[in] RH	    Relative Humidity (%)
		* @param[in] wind	Windspeed at 2.5 m, m s-1
		* @param[in] Press	Atmospheric pressure (kpa m-2)
		* @param[in] ConstantTemperature boolian if true, leaf temperature=air temperature when calculating gas exchange
		\return nothing
		*/


	{
		this->PhotoFluxDensity = PhotoFluxDensity;
		double PAR = (PhotoFluxDensity/4.55); //PAR is watts m-2
		double NIR = PAR; // If total solar radiation unavailable, assume NIR the same energy as PAR waveband
		this->R_abs = (1-sParms.scatt)*PAR + 0.15*NIR + 2*(epsilon*sbc*pow(Tair+273,4)); // times 2 for projected area basis
		// shortwave radiation (PAR (=0.85) + NIR (=0.15) solar radiation absorptivity of leaves: =~ 0.5
		//transfer variables to local scope
		this->CO2 = CO2;
		this->RH = __min(100.0, __max(RH, 20.0))/100; /* made this 20 from 10 for testing DT 2/21*/
		this->Tair = Tair;
		this->wind = wind;  //m s-1
		this->Press = Press; //kPa
		ConstantLeafTemperature=0;
		this->leafp = leafp;
		//TODO do we need GetParms as in the original?
		this->leafpEffect = 1; //At first assume there is not drought stress, so assign 1 to leafpEffect. Yang 8/20/06
		GasEx(leafp,ET_supply);   // Gas exchange calculations here
	}

	void CGasExchange::GasEx(double leafp, double ET_supply)
		/** 
		* carries out calculations for photosynthesis and stomatal conductance.
		* no parameters, returns nothing

		* @see SearchCi(), @see EnergyBalance(), @see CalcStomatalConductance()
		\return nothing
		*/
	{
		double Tleaf_old;  //previous leaf temperture (for iteration)
		int   iter=1;
		iter_total=0;
		Tleaf = Tair; Tleaf_old = 0;
		Ci = sParms.internalCO2Ratio*CO2;
		BoundaryLayerConductance = CalcTurbulentVaporConductance();
		StomatalConductance = CalcStomatalConductance(leafp);
		while ((abs(Tleaf_old -Tleaf)>0.01) && (iter < maxiter))
		{
			Tleaf_old=Tleaf;
			Ci=SearchCi(Ci);
			StomatalConductance=CalcStomatalConductance(leafp);
			EnergyBalance();
			iter2 =++iter; //iter=iter+1, iter2=iter; 
			if (ConstantLeafTemperature) Tleaf=Tair;
		} 

	}
	

	void CGasExchange::Photosynthesis(double Ci)    
		/**
		* Calculates photosynthesis for C4 plants. 
		* Requires Incident PhotoFluxDensity, Air temp in C, CO2 in ppm, RH in percent
		@see SetVal()
		* 
		@param[in] Ci - internal CO2 concentration, umol mol-1

		\return nothing
		*/
	{
		const double    curvature=0.995; //!< \b curvature factor of Av and Aj colimitation
		const long       Eao = 36000;  //*!< \b EAO, activation energy for Ko */
		const int	    Vpr25 = 80;    //*!<   \b Vpr25, PEP regeneration limited Vp at 25C, value adopted from vC book */
		const double    x = 0.4;       //*!< \b x Partitioning factor of J, yield maximal J at this value */
		const double    alpha = 0.0001; //*!< \b alpha, fraction of PSII activity in the bundle sheath cell, very low for NADP-ME types  */
		const double    beta = 0.99;   //*!< \b beta, smoothing factor */
		const double    theta = 0.5;

		double Kp, Kc, Ko, Km;         //!<\b Kp, \b Kc, \b Ko, \b Km, Calculated Michaelis params as a function of temperature
		double Ia, I2;                 // secondary calculated light variables
		double Vpmax, Jmax, Vcmax, Eac, Rm, Om,  J, Ac1, Ac2, Ac, Aj1,
			Aj2, Aj, Vp1, Vp2, Vp, P,  Ca, Cm, Vpr,
			Os, GammaStar, Gamma, a1, b1, c1; //secondary calculated variables

		//* Light response function parameters */
		Ia = PhotoFluxDensity*(1-sParms.scatt);    //* absorbed irradiance */
		I2 = Ia*(1-sParms.f)/2;    //* useful light absorbed by PSII */
		//* other input parameters and constants */
		P  = Press/100;
		Ca = CO2*P; //* conversion to partial pressure Atmospheric partial pressure of CO2, kPa*/
		Om = O;   //* mesophyle O2 partial pressure */
		Eac=sParms.EaVc;




		Kp = sParms.Kp25*pow(Q10,(Tleaf-25.0)/10.0);
		Vpr = Vpr25*pow(Q10,(Tleaf-25.0)/10.0);
		Kc = sParms.Kc25*exp(Eac*(Tleaf-25)/(298*R*(Tleaf+273))); //Kc adjusted for temperature
		Ko = sParms.Ko25*exp(Eao*(Tleaf-25)/(298*R*(Tleaf+273)));
		Km = Kc*(1+Om/Ko); //* effective M-M constant for Kc in the presence of O2 */
		DarkRespiration = sParms.Rd25*exp(sParms.Ear*(Tleaf-25)/(298*R*(Tleaf+273)));
		// The following are Arrhenius Equations for parameter temperature dependencies
		// Vpm25 (PEPC activity rate) , Vcm25  (Rubisco Capacity rate) and Jm25 (Whole chain electron transport rate) are the rates at 25C for Vp, Vc and Jm
		
		double CriticalNitrogen;
		CriticalNitrogen = __max(0.25, lfNContent);
		double Vcm25_L = sParms.Vcm25*(2 / (1 + exp(-2.9*(CriticalNitrogen - 0.25))) - 1);
		double Jm25_L = sParms.Jm25*(2 / (1 + exp(-2.9*(CriticalNitrogen - 0.25))) - 1);
		double Vpm25_L = sParms.Vpm25*(2 / (1 + exp(-2.9*(CriticalNitrogen - 0.25))) - 1); //in Sinclair and Horie, 1989 Crop sciences, it is 4 and 0.2; 


		Vpmax = Vpm25_L*exp(sParms.EaVp*(Tleaf-25)/(298*R*(Tleaf+273)));
		Vcmax = Vcm25_L*exp(sParms.EaVc*(Tleaf-25)/(298*R*(Tleaf+273)));
		Jmax = Jm25_L*exp((((Tleaf+273)-298)*sParms.Eaj)/(R*(Tleaf+273)*298))*(1+exp((sParms.Sj*298-sParms.Hj)/(R*298)))
			/(1+exp((sParms.Sj*(Tleaf+273)-sParms.Hj)/(R*(Tleaf+273.0))));
		Rm = 0.5*DarkRespiration;

		Cm=Ci; //* mesophyle CO2 partial pressure, ubar, one may use the same value as Ci assuming infinite mesohpyle conductance */
		double gs_last=0;

		StomatalConductance = CalcStomatalConductance(leafp);
		Vp1 = (Cm*Vpmax)/(Cm+Kp); //* PEP carboxylation rate, that is the rate of C4 acid generation  Eq 1 in Kim 2007*/
		Vp2 = Vpr;
		Vp = __max(__min(Vp1, Vp2),0);
		//* Enzyme limited A (Rubisco or PEP carboxylation */
		Ac1 = (Vp+sParms.gbs*Cm-Rm);
		Ac2 = (Vcmax-DarkRespiration);
		//* Quadratic expression to solve for Ac */
		a1 = 1-(alpha/0.047)*(Kc/Ko);
		b1 = -(Ac1 + Ac2 + sParms.gbs*Km + (alpha/0.047)*(sParms.gamma1*Vcmax + DarkRespiration*Kc/Ko));
		c1 = Ac1*Ac2-(Vcmax*sParms.gbs*sParms.gamma1*Om+DarkRespiration* sParms.gbs*Km);
		//Ac = QuadSolnLower(a1,b1,c1);
		Ac = __min(Ac1,Ac2);
		//* Light and electron transport limited  A mediated by J */
		//J=minh(I2,Jmax,theta);  //* rate of electron transport */
		J = QuadSolnLower(theta, -(I2 + Jmax), I2*(Jmax)); //* rate of electron transport */

		Aj1 = (x*J/2-Rm+ sParms.gbs*Cm);  // Eq 4 in Kim, 2007
		Aj2 = (1-x)*J/3-DarkRespiration;       //Eq 4 in Kim, 2007
		Aj = __min(Aj1,Aj2);      //Eq 4 in Kim, 2007
		AssimilationNet = ((Ac+Aj) - sqrt(Square(Ac+Aj)-4*beta*Ac*Aj))/(2*beta); //* smooting the transition between Ac and Aj */
		//AssimilationNet=minh(Ac,Aj, curvature);
		gs_last=StomatalConductance;
		Os = alpha*AssimilationNet/(0.047*sParms.gbs)+Om; //* Bundle sheath O2 partial pressure, mbar */
		GammaStar = sParms.gamma1*Os;
		Gamma = (DarkRespiration*Km + Vcmax*GammaStar)/(Vcmax-DarkRespiration);
		AssimilationGross = __max(0, AssimilationNet + DarkRespiration); 

	}


	void CGasExchange::EnergyBalance()
		/** 
		Calculates Transpiration rate (T) and leaf temperature (Tleaf). Iterates by recalculating
		photosynthesis until leaf temperatures converge

		See Campbell and Norman (1998) pp 224-225
        
		Because Stefan-Boltzman constant is for unit surface area by denifition,
		all terms including sbc are multilplied by 2 (i.e., RadiativeConductance, thermal radiation)
		
		Does not have input 

		\return nothing but calculates transpiration (T) and leaf temperature (Tleaf)

		
		*/

	{
		const long Lambda = 44000; //latent heat of vaporization of water J mol-1 - not used in this implementation
		const double Cp = 29.3; // thermodynamic psychrometer constant and specific hear of air, J mol-1 C-1
		const double psc = 6.66e-4; //psycrometric constant units are C-1
		//psc=Cp/Lambda = 29.3/44000 See Campbell and Norman, pg 232, after eq 14.11

		//The following are secondary variables used in the energy balance
		double HeatConductance,  //heat conductance J m-2 s-1
			VaporConductance, //vapor conductance ratio of stomatal and heat conductance mol m-2 s-1
			RadiativeConductance, //radiative conductance J m-2 s-1
			RadiativeAndHeatConductance, //radiative+heat conductance
			psc1,  // apparent psychrometer constant Campbell and Norman, page 232 after eq 14.11
			Ea,   //ambient vapor pressure kPa
			thermal_air; // emitted thermal radiation Watts  m-2
		double lastTi, newTi;
		int    iter;
		bool badval = false;

		HeatConductance = BoundaryLayerConductance*(0.135/0.147);  // heat conductance, HeatConductance = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Boundary Layer Conductance to Heat
		// Since BoundaryLayerConductance is .147*sqrt(u/d) this scales to 0.135*sqrt(u/d) - HeatConductance on page 109 of Campbell and Norman, 1998
		// Wind was accounted for in BoundaryLayerConductance already  as BoundaryLayerConductance (turbulent vapor transfer) was calculated from CalcTurbulentVaporConductance() in GasEx. 
		// units are J m-2 s-1 
		VaporConductance = StomatalConductance*BoundaryLayerConductance/(StomatalConductance+BoundaryLayerConductance);      //vapor conductance, StomatalConductance is stomatal conductance and is given as gvs in Campbell and Norman.      
		                                                                                                                    // note units are moles m-2 s-1. 
		RadiativeConductance = (4*epsilon*sbc*pow(273+Tair,3)/Cp)*2; // radiative conductance, *2 account for both sides
		RadiativeAndHeatConductance = HeatConductance + RadiativeConductance;
		thermal_air = epsilon*sbc*pow(Tair+273,4)*2; //Multiply by 2 for both surfaces
		psc1 = psc*RadiativeAndHeatConductance/VaporConductance; 
		this->VPD = Es(Tair)*(1-RH); // vapor pressure deficit Es is saturation vapor pressure at air temperature
		// iterative version
		newTi=-10;
		iter=0;
		lastTi=Tleaf;
		double Res, dRes; //temporary variables
		double thermal_leaf;
		Ea = Es(Tair)*RH; // ambient vapor pressure
		while ((abs(lastTi-newTi)>0.001) && (iter <maxiter)) 
		{
			lastTi=newTi;
			
				 // the commented line is from the original but doesn't make much difference for dry matter
            Tleaf= Tair + (R_abs- thermal_air-Lambda*VaporConductance*this->VPD/Press)/(Cp*RadiativeAndHeatConductance+Lambda*Slope(Tair)*VaporConductance); // eqn 14.6a
			if (std::isnan(Tleaf)) badval = true;
			thermal_leaf=epsilon*sbc*pow(Tleaf+273,4)*2;
			Res = R_abs - thermal_leaf - Cp*HeatConductance*(Tleaf - Tair) - Lambda*VaporConductance*0.5*(Es(Tleaf)-Ea)/Press; // Residual function: f(Ti), KT Paw (1987)
			dRes= -4*epsilon*sbc*pow(273+Tleaf,3)*2-Cp*HeatConductance*Tleaf-Lambda*VaporConductance*Slope(Tleaf); // derivative of residual: f'(Ti)
			newTi = Tleaf + Res/dRes; // newton-rhapson iteration
			iter++;
		}
		Tleaf=newTi;

		Transpiration =__max(0,VaporConductance*((Es(Tleaf) - Ea) / Press) / (1 - (Es(Tleaf) + Ea) / (Press))); //Don't need Lambda - cancels out see eq 14.10 in Campbell and Norman, 1998
		// mol m-2 s-1. everything in moles, units of VaporConductance are moles. 
		// the above was divided by (1-Es(Tleaf)+Ea)/(press))
	}



	double CGasExchange::CalcStomatalConductance(double pressure)  
		/**
		* calculates and returns stomatal conductance for CO2 in umol CO2 m-2 s-1 
		* Uses Ball-Berry model.
		* @see Es()
		\return stomatal conductance for CO2 umol m-2 s-1
		*
		*
		*/
		//*! \fn CalcStomatalConductance()
	{
		//** \code
		double Ds, //! \b Ds, VPD at leaf surface 
			aa,    //! \b aa, a value in quadratic equation 
			bb,    //! \b bb, b value in quadratic equation 
			cc,    //! \b cc, calcuation variable (x) in quadratic equation
			hs,    //! \b hs, solution for relative humidity
			Cs,    //! \b Cs, estimate of mole fraction of CO2 at the leaf surface
			Gamma, //! \b Gamma, CO2 compensation point in the absence of mitochondirial respiration, in ubar
			StomatalConductance;    //! \b StomatalConductance, temporary variable to hold stomatal conductance
		Gamma = 10.0; 
		//** \endcode

		double temp = set_leafpEffect(pressure);
		double P=Press/100;  
		Cs = (CO2 - (1.37*AssimilationNet/BoundaryLayerConductance))*P; // surface CO2 in mole fraction
		if (Cs == Gamma) Cs = Gamma + 1;
		if (Cs <= Gamma) Cs = Gamma + 1;
		// Quadratic equation to obtain hs by combining StomatalConductance with diffusion equation
		aa = temp*sParms.g1*AssimilationNet/Cs;
		bb = sParms.g0+BoundaryLayerConductance-(temp*sParms.g1*AssimilationNet/Cs);
		cc = (-RH*BoundaryLayerConductance)-sParms.g0;
		hs = QuadSolnUpper(aa,bb,cc);
		if (hs > 1) hs = 1;  else if (hs <= 0.30) hs = 0.30; //preventing bifurcation
		if (hs<0) hs = 0;
		Ds = (1-hs)*Es(Tleaf); // VPD at leaf surface
		StomatalConductance = (sParms.g0+temp*sParms.g1*(AssimilationNet*hs/Cs));
		if (StomatalConductance < sParms.g0) StomatalConductance=sParms.g0; //Limit StomatalConductance to mesophyll conductance 
	// this below is an example of how you can write temporary data to a debug window. It can be copied and 
	// pasted into excel for plotting. Dennis See above where the CString object is created.
	//CString strdbg;
	//strdbg.Format("Source = %s tmp = %f pressure = %f Ds= %f Tleaf = %f Cs = %f Anet = %f hs = %f RH = %f\n",
		// sType, tmp, pressure, Ds, Tleaf, Cs, A_net, hs, RH);
	//OutputDebugString(strdbg);
		return StomatalConductance;  // moles m-2 s-1
	}





	double CGasExchange::CalcTurbulentVaporConductance(void)
	{
		/**
		* calculates conductance for turbulant vapor transfer in air - forced convection
		\return conductance for turbulent vapor transfer in air (mol m-2 s-1)
		*/

		double ratio; /*!< temporary holding variable for stomatal ratio calculations*/
		double Char_Dim; /*!<  characteristic dimension of leaf */
		ratio = Square(sParms.stomaRatio+1)/(Square(sParms.stomaRatio)+1);
		Char_Dim = sParms.LfWidth*sParms.widthPara; // characteristic dimension of a leaf, leaf width in m
		// wind is in m per second
		return (1.4*0.147*sqrt(__max(0.1,wind)/Char_Dim))*ratio; 
		// multiply by 1.4 for outdoor condition, Campbell and Norman (1998), p109, gva
		// multiply by ratio to get the effective blc (per projected area basis), licor 6400 manual p 1-9
	}

	double CGasExchange::Es(double Temperature) 
	{
		/**
		* calculates and returns Saturation vapor pressure (kPa). Campbell and Norman (1998), p 41. 
		@param[in] Temperature
		\return saturated vapor pressure (kPa)
		*/

		double result;
		// a=0.611 kPa, b=17.502 C and c=240.97 C 
		//Units of Es are kPa
		result=(0.611*exp(17.502*Temperature/(240.97+Temperature)));
		return result;
	}

	double CGasExchange::Slope(double Temperature) 
		/**
		Calculates the slope of the sat vapor pressure curve: 
		first order derivative of Es with respect to T

		@param[in] Temperature (C)
		@see Es()
		\return slope of the vapor pressure curve kPa T-1
		*/

	{
		double VPSlope;
		// units of b and c are  degrees C
		const double b= 17.502; const double c= 240.97;
		VPSlope=(Es(Temperature)*(b*c)/Square(c+Temperature)/Press);
		return VPSlope; 
	}
	
	double CGasExchange::set_leafpEffect(double pressure)   //calculating effect of leaf water potential on stomatal conductance
	{                                                        //model of Tuzet et al. 2003 used  Yang 8/21/06
		//pressure - leaf water potential MPa...
		// from Tuzet et al., 2003
		this->leafpEffect = (1 + exp(sParms.sf*sParms.phyf)) / (1 + exp(sParms.sf*(sParms.phyf - pressure)));
		double temp = this->leafpEffect;
		return temp;
	}
	double CGasExchange::SearchCi(double CO2i)
	{
		/**
		* does a secant search to find the optimal internal CO2 concentration (ci)
		* Calls:
		* @see EvalCi()
		@param[in] CO2i - internal CO2 concentration, (umol mol-1)
		\return Ci (umol mol-1)
		*/

		int iter;
		double fprime, Ci1, Ci2, Ci_low, Ci_hi, Ci_m;
		double temp;
		Ci1 = CO2i;
		Ci2 = CO2i + 1.0;
		Ci_m = (Ci1+Ci2)/2.0;
		iter_Ci = 0;
		iter = 0;
		isCiConverged = true;

		do 
		{
			iter++;
			//Secant search method
			if (abs(Ci1-Ci2) <= errTolerance) {break;}
			if (iter >= maxiter) 
			{
				isCiConverged = false;
				break;
			}
			fprime = (EvalCi(Ci2)-EvalCi(Ci1))/(Ci2-Ci1);  // f'(Ci)
			if (fprime != 0.0) 
			{
				Ci_m = (std::max)(errTolerance, Ci1-EvalCi(Ci1)/fprime); // use (std::max) because somewhere max is defined as a macro
			}
			else
				Ci_m = Ci1;
			Ci1 = Ci2;
			Ci2 = Ci_m;
			temp=EvalCi(Ci_m);
			double temp2=maxiter;
		} while ((abs(EvalCi(Ci_m)) >= errTolerance) || (iter < maxiter));




		// C4 photosynthesis fails to converge at low soil water potentials using secant search, 6/8/05 SK
		// Bisectional type search is slower but more secure
		//Bisectional search
		if (iter > maxiter)
		{
			Ci_low = 0.0;
			Ci_hi = 2.0*CO2;
			isCiConverged = false;

			while (abs(Ci_hi-Ci_low) <= errTolerance || iter > (maxiter*2))
			{
				Ci_m = (Ci_low + Ci_hi)/2;
				if (abs(EvalCi(Ci_low)*EvalCi(Ci_m)) <= eqlTolerance) break;
				else if (EvalCi(Ci_low)*EvalCi(Ci_m) < 0.0) {Ci_hi = (std::max)(Ci_m, errTolerance);}
				else if (EvalCi(Ci_m)*EvalCi(Ci_hi) < 0.0)  {Ci_low =(std::max)(Ci_m, errTolerance);}
				else {isCiConverged = false; break;}
			}

		}

		CO2i = Ci_m;
		Ci_Ca = CO2i/CO2;
		iter_Ci = iter_Ci + iter;
		iter_total = iter_total + iter;
		return CO2i;

	}
	double CGasExchange::EvalCi(double Ci)
	{
		/**
		* Called by SearchCi() to calculate a new value of Ci for the current values of photosynthesis and stomatal conductance
		* determined using parameters from a previous step where the energy balance was solved.
		@see SearchCi()

		@param[in] Ci **, estimate of internal CO2 concentration, umol mol-1

		\return the difference between the passed value of Ci (old)and the new one. 

		*/
		double newCi;

		Photosynthesis(Ci);
		if (abs(StomatalConductance) > eqlTolerance) 
		{
			newCi = max(1.0,CO2 - AssimilationNet*(sParms.SC_param/StomatalConductance+ sParms.BLC_param/BoundaryLayerConductance)*(Press/100.0));
		}
		else
			newCi = max(1.0,CO2 - AssimilationNet*(sParms.SC_param/eqlTolerance+ sParms.BLC_param/BoundaryLayerConductance)*(Press/100.0));
		return (newCi-Ci);
	}

	//These two functions solve the quadratic equation.
	double CGasExchange::QuadSolnUpper (double a, double b, double c )
	{

		/** solves the uppper part of the quadratic equation ax2+bx2=c

		@param[in] a
		@param[in] b
		@param[in] c 

		\return lower portion of x		*/
		if (a==0) return 0;
		else if ((b*b - 4*a*c) < 0) return -b/a;   //imaginary roots
		else  return (-b+sqrt(b*b-4*a*c))/(2*a);
	}

	double CGasExchange::QuadSolnLower (double a, double b, double c )
	{
		/** solves the lower part of the quadratic equation ax2+bx=c

		@param[in] a
		@param[in] b
		@param[in] c 
		\return lower portion of x
		*/
		if (a==0) return 0;
		else if ((b*b - 4*a*c) < 0) return -b/a;   //imaginary roots
		else  return (-b-sqrt(b*b-4*a*c))/(2*a);
	}

	//*! hyperbolic min
	double CGasExchange::minh(double fn1,double fn2,double theta2)
	{

		/**  
		@param [in] fn1 first value to be compared for min
		@param [in] fn2 second value to be compared for min
		@param [in] theta2  curvature factor

		\return hyperbolic minimum
		*/
		double x, res;

		x = ((fn1+fn2)*(fn1+fn2)-4*theta2*fn1*fn2);
		if (x<0)
		{
			res = __min(fn1,fn2); 
			return res;
		}
		if (theta2==0.0)
		{
			res= fn1*fn2/(fn1+fn2);
			return res;
		}
		else
		{
			res = ((fn1+ fn2) - sqrt(x))/(2*theta2); // hyperbolic minimum
			return res;
		}
	}


