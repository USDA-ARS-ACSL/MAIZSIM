#include "stdafx.h"
#include "leaf.h"
#include "weather.h"
#include "initinfo.h"
#include <cmath>
#include <algorithm>
#define MINUTESPERDAY (24*60);
using namespace std;


CLeaf::CLeaf(int n, CDevelopment * dv): COrgan()
{
	rank = n;
	length=width = area = PotentialArea = plastochrons = 0.0;
	SLA = 200.0; // temporary for now - it should vary by age. Value comes from some of Soo's work
	                
	PotentialAreaIncrease = 0;
	RelativeAreaIncrease=0;
	FullyExpandedArea=0;
	totalLeaves = dv->get_totalLeaves(); //todo: should be a plant parameter not leaf
	TotalDroppedLeaves=0;
	actualArea = actualgreenArea=greenArea = senescentArea = droppedArea = 0;
	actualLength=actualwidth=0;
	initiated = growing = prolific = aging = terminated = dropped = false;
	phase1Delay = phase2Duration=ptnLength=ptnWidth = 0;
	elongRate = 0.564;  //0.564 cm dd-1 Fournier and Andrieu 1998 Pg239. This is the "potential" elongation rate with no water stress Yang
	elongAge = 0.0;
	old_leaf = 0; //record leaf area in the last time step
    N_content = 3.0; //no N stress
	LeafCalibTemperature=dv->get_CalibTemperature();   
	

		
}

void CLeaf::initialize (CDevelopment * dv) 
// set potential leaf area of the current rank based on generic leaf no, see Fournier and Andrieu (1998), Birch et al., (1998)
// Routine below to calculate potentialArea will repeat in CLeaf::update to get actual leaf area after the inductive phase to any additional leaves developed before tassel initiation, SK 1-20-12
{

	COrgan::initialize();
	const double LM_min = 115, k = 24.0;
	double L_max = (sqrt(LM_min*LM_min + k*(totalLeaves - dv->get_initInfo().genericLeafNo)));
	// If this routine runs before TI, totalLeaves = genericLeafNo, SK
	double n_m, a, b;

	n_m = 5.93 + 0.33*totalLeaves; // the rank of the largest leaf. YY
	a = -10.61 + 0.25*totalLeaves;
	b = -5.99 + 0.27*totalLeaves;
    //equation 7 in Fournier and Andrieu (1998). YY
	ptnLength = L_max*exp(a/2*pow(rank/n_m-1,2)+b/2*pow(rank/n_m-1,3)); //equa 8(b)(Actually eqn 6? - eqn 8 deals with leaf age - DT)
	                                                                   //in Fournier and Andrieu(1998). YY
	phase2Duration = ptnLength/elongRate; // linear phase duration in physiological time equa 8(a)Fournier and Andrieu(1998) YY
    phase1Delay = __max(0.0, -5.16+1.94*rank); //Fournier's value : -5.16+1.94*rank;equa 11 Fournier and Andrieu(1998) YY
	//phase1Delay = __max(0.0, -5.16+2.14*rank);
	double W_max = L_max*0.106;//Fournier and Andrieu(1998) Pg242 YY
	double LA_max = L_max*W_max*0.75; // daughtry and hollinger (1984) Fournier and Andrieu(1998) Pg242 YY
    double lfno_effect = max(0.5, min(1.0, exp(-1.17+0.047*totalLeaves))); // Fig 4 of Birch et al. (1998) 
	PotentialArea=lfno_effect*LA_max*exp(a*pow(rank/n_m-1,2)+b*pow(rank/n_m-1,3)); //equa 6. Fournier and Andrieu(1998) multiplied by Birch et al. (1998) leaf no effect
	                                                               //LA_max the area of the largest leaf
	                                                               //PotentialArea potential final area of a leaf with rank "n". YY

	initiated = true;
	growing=true; // DT
	first=true;
	std::cout << " GDDay " << dv->get_GDDsum()  << " leaf " << dv->get_LvsInitiated() << " " << rank << " totalLeaves " << totalLeaves 
              <<  " LAmax " << LA_max << " lfno_eff " << lfno_effect << " Area " << PotentialArea << std::endl;
}


void CLeaf::update(CDevelopment * dv, double predawnlwp)
{ 
	COrgan::set_temperature(dv->get_Tcur());
	COrgan::update();
	{
	// The routine within this bracket is the same as in initialize but repeats here to account for updated totalLeaves, SK
	// Can become a separate function to be called twice in both processes, SK
	//const double phyllochron = 39.60; //PI3733, (Topt-Tbase)/Rmax = (31-8)/0.581
	                                  // should this be a variety parameter?
	                                  // dennis - trying a lower value for Colorado data
	                                  //using rmax =0.394 (31-8)/0.394
	const double LM_min = 115, k = 24.0;
	                                     //LM_min is a length characteristic of the longest leaf,in Fournier and Andrieu 1998, it was 90 cm
	                                     //LA_max is a fn of leaf no (Birch et al, 1998 fig 4) with largest reported value near 1000cm2. This is implemented as lfno_effect below, SK
	                                     //LM_min of 115cm gives LA of largest leaf 1050cm2 when totalLeaves are 25 and Nt=Ng, SK 1-20-12
	                                     // Without lfno_effect, it can be set to 97cm for the largest leaf area to be at 750 cm2 with Nt ~= Ng (Lmax*Wmax*0.75) based on Muchow, Sinclair, & Bennet (1990), SK 1-18-2012
	                                     // Eventually, this needs to be a cultivar parameter and included in input file, SK 1-18-12
	                                     // the unit of k is cm^2 (Fournier and Andrieu 1998 Pg239). YY
	                                     // T_peak is the temperature at which potential leaf size is maximized based on Kim et al (2007) EEB data Fig 1a and 1b,
	                                     // T_peak doesn't seem to be cultivar specific, 1/18/12 --SK
	/* 
	  T_peak, formerly T_opt, has been re-examined and now determined based on combined and normalized data from Tollennar (1989), Bos (2000), Hesketh (1989), Fournier (1998) and Kim (2007)
	  T_peak is determined based on a base temperature (Tb) of 8.0 deg C (Birch et al., 2003, EJA).
	  Calibration with above data confirms Tb = 8.03 C.
	  Todo: Probably better to put these constants in the header for the leaf class
	   -- 1/17/2012, SK
	*/

	double L_max = (sqrt(LM_min*LM_min + k*(totalLeaves - dv->get_initInfo().genericLeafNo)));
// L_max is the length of the largest leaf when grown at T_peak. Here we assume LM_min is determined at growing Topt with minmal (generic) leaf no, SK 8/2011
// Revised to remove this adjustment assuming LM_min was determined near T_peak 1/17/12, SK
// L_max is multiplied by leaf rank fn, T response, etc. 
//	double L_max = (sqrt(LM_min*LM_min + k*(totalLeaves - dv->get_initInfo().genericLeafNo)));
          //equation 8(c)  in Fournier and Andrieu(1998) the "genericLeafNo" is the number of
	       // leaf induced at the start of the inductive period. YY
	double n_m, a, b;

	n_m = 5.93 + 0.33*totalLeaves; // the rank of the largest leaf. YY
	a = -10.61 + 0.25*totalLeaves;
	b = -5.99 + 0.27*totalLeaves;
    //equation 7 in Fournier and Andrieu (1998). YY
	//Some of Fournier and Andrieu's variables and equations are not used at this time
	ptnLength = L_max*exp(a/2*pow(rank/n_m-1,2)+b/2*pow(rank/n_m-1,3)); //equa 8(b)(Actually eqn 6? - eqn 8 deals with leaf age - DT)
	                                                                   //in Fournier and Andrieu(1998). YY
	phase2Duration = ptnLength/elongRate; // linear phase duration in physiological time equa 8(a)Fournier and Andrieu(1998) YY /*Not used*/
    phase1Delay = __max(0.0, -5.16+1.94*rank); //Fournier's value : -5.16+1.94*rank;equa 11 Fournier and Andrieu(1998) YY  /*Not used*/
	//phase1Delay = __max(0.0, -5.16+2.14*rank);
	double W_max = L_max*0.106;//Fournier and Andrieu(1998) Pg242 YY
	double LA_max = L_max*W_max*0.75; // daughtry and hollinger (1984) Fournier and Andrieu(1998) Pg242 YY
    double lfno_effect = max(0.5, min(1.0, exp(-1.17+0.047*totalLeaves))); // Fig 4 of Birch et al. (1998) 
	PotentialArea=lfno_effect*LA_max*exp(a*pow(rank/n_m-1,2)+b*pow(rank/n_m-1,3)); //equa 6. Fournier and Andrieu(1998) multiplied by Birch et al. (1998) leaf no effect
	                                                               //LA_max the area of the largest leaf
	                                                               //PotentialArea potential final area of a leaf with rank "n". YY

	}

	calcLongevity(predawnlwp);
	if (!dropped)
		{
		  Expand(dv, predawnlwp);
	      senescence(dv, predawnlwp);
		}
   
	double effect=LWPeffect(predawnlwp);
	if (area >= 0.99*PotentialArea*effect) prolific = true; else prolific = false; //assuming the lenght of the expanding
	                                                                        //period for each leaf is not affected by
	                                                                       //draught stress. Then the total leaf area
	                                                                       //of each leaf is affected by leaf water potential
	                                                                       //the same way as the expansion rate? YY
	if (rank == 0) greenArea = 0.0; // coleoptile doesn't photosynthesize
	else 
	{
		if (terminated && !dropped) // make sure we only enter this loop once for each leaf
	    {   
		    droppedArea=area;
			area = greenArea = 0;
			dropped=true; 

		}
			 // dropped
	else
	{
		greenArea = max(0.,area-senescentArea);  //To do- I don't know why we need two variables for this (DT)
		actualgreenArea = greenArea;
		FullyExpandedArea=max(area,FullyExpandedArea);
	}
	}


}


void CLeaf::Expand(CDevelopment * dv, double predawnlwp)
//introduing the effect of predawn leaf water potential on 
//the potential leaf area growth based on Lizaso et al. (2003) YY
// called after emergence
{
	const double T_peak = 18.7, Tb = 8.0, phyllochron = (dv->get_T_Opt()- Tb)/(dv->get_Rmax_LTAR()); // this corresponds to PHY in Lizaso (2003); phyllochron needed for next leaf appearance in degree days (GDD8) - 08/16/11, Soo. 
	                                     // T_opt is the temperature at which potential leaf size is maximized based on Kim et al (2007) EEB data Fig 1a and 1b,
	                                     // T_opt can be cultivar specific, 8/16/11 --SK
	                                     // Changed T_opt to T_peak -- 1/17/12 SK
		const double T_opt_ke  = 18.99, T_max_ke = 42.25; //Beta fn parms for ke, determined using SPAR data from Kim et al. (2007) EEB, SK 08-18-2011      
	double Wk = totalLeaves/8.18;
	const double tt2 = phyllochron; //GDD sum from emegence to second tip appearance, guessed
	double tt = (rank-2)*phyllochron + tt2;
//	double temp;
//	temp=dv->get_GDDsum();
	double T_cur = dv->get_Tcur(); // get current temperature;
    double Ke_Topt, Ke, te;
	{
		Ke_Topt = 0.026+0.174*exp(-(pow((double)rank-1,2)/(2*pow(Wk,2)))); // Eqn 8 in Lizaso et al. 
														   // note that ko=.02  (from figure 2) and kx=0.2
														   // rank is leaf number. This is a 
														   //slope parameter
														   // ke_Topt: Assume this value represent intrinsic growth rate at T_opt -- SK 1-17-12
														   // values corrected (k0=0.026 and kx= 0.174) to reflect the figure 2 -- SK 1-17-12

		Ke = Ke_Topt *(dv->beta_fn(T_cur, 1.0, T_opt_ke, T_max_ke)); 
		// Assuming Ke values in Lizaso (2003) paper were determined at 20C of mean growing season temperature
		// Normalize Ke values in reference to te at 20C (Ke20), SK, 08-18-2011; Revised to remove this normalization assuming Lizaso values represent optimal values -- SK 1-17-12
		// Ke=Ke_Topt;

		/*
		  Turned off T effects on ke. This line needs to be commented out to turn on temperature effect. SK
          Turned off because
		  1) we lack estimates for ke_Topt, and
		  2) T dependence determined from whole-plant leaf area responses could be quite different from individual leaf growth responses per review 2's comments
		   -- SK 1-17-12
		*/
	}
	{
		double T_effect_te;
		if (rank > 2)
		{
			const double Q10_te = 1.546; //Q10 value for te, determined using SPAR data from Kim et al. (2007) EEB, SK 08-18-2011; Revised value, 1-17-11,SK
			const double T_opt_te = T_opt_ke; //Assume T_opt_ke represents temperature condition for te_Top determination, SK 1-17-12

			T_effect_te = pow(Q10_te, (T_cur-T_opt_te)/10.0);
		// Assuming te values in Lizaso (2003) paper were determined at Topt of mean growing season temperature
		// Normalize te values with respect to te at 20C (te20), SK, 08-18-2011; Revised to set the effect = 1.0 at T_opt, 1-17-11 SK
		// With this, extreme high T will make it to take longer to reach 50% expansion, SK */
	       
			// T_effect_te = 1.0; //comment it out to turn off the effect, SK
			te = tt + (2.197/Ke_Topt)*T_effect_te; 
											   //Eq 6 in Lizaso, te is time to 50% expansion
											   //tt is time to appearance of leaf tip in thermal time.
											   //Leaf appears when the leaf blade
											   // is approx 10% expanded
											   // Adjust te with respect to T_cur. T_effect = 1 at optimal T, SK 1-18-12
											   // This accounts for longer growth duration at high T and applies only to expansion portion (not tt), SK
											   // T effect on tt is taken care of in phyllochron determination using beta fn, SK
			                                   // Ke_Topt (not Ke) is used as denominator because te values in Lizaso are assumed to come from optimal growth conditions, SK 1-18-12
			                                   
		                                       
		}            
		else te = 25.0*rank;
	}

	double dL=0.0;
	double actualDL = 0.0;
	double water_effect=LWPeffect(predawnlwp);
	if (first)
	{
		first=false;
	}

	double T_effect_size = max(0.0, (T_cur-Tb)/(T_peak-Tb)*exp(1.0-(T_cur-Tb)/(T_peak-Tb))); // adjust T effect on leaf size, 8/2011, SK; revised it to adjust with respect to Tb, SK 1-17-12
	//T_effect = water_effect = CarbonRatio = 1.0; for debugging only - assumes no stress.

    PotentialAreaIncrease =max(0.0,PotentialArea*Ke*exp(-Ke*(dv->get_GDDsum()-te))/pow(1+exp(-Ke*(dv->get_GDDsum()-te)),2)*dv->get_dGDD());
	                                           //potential leaf area increase without carbon limitation YY
	dL = min(water_effect,T_effect_size)*PotentialAreaIncrease;
	//	PotentialArea*Ke*exp(-Ke*(dv->get_GDDsum()-te))/pow(1+exp(-Ke*(dv->get_GDDsum()-te)),2)*dv->get_dGDD();
	// water effect and T effect should probably multiplicative, SK 1-18-12
	                         
	actualDL = dL;
	area +=dL;
	actualArea += actualDL;
	if (dv->get_GDDsum() < 2*te) growing = true; else growing = false;
	//set_PotentialCarboIncrement(PotentialAreaIncrease/SLA); //SK 8/20/2010: Doesn't this have to be set before CarbonRatio has been read?
	                                                  // DT 1/11/2012: this is for the next time step.
 

	return;
}



void CLeaf::senescence(CDevelopment * dv, double predawnlwp)
//potential leaf area growth based on Lizaso et al. (2003)
 //ToDo - need to add water stress
{
    const double T_opt_ks  = 18.99, T_max_ks = 42.25; //Beta fn parms for ks, determined using SPAR data from Kim et al. (2007) EEB, SK 08-18-2011      
	// not ks is senescent form of ke

	const double T_peak = 18.7;
	double water_effect=LWPeffect(predawnlwp);
  
	const double Tb = 8.0, phyllochron = (dv->get_T_Opt()- Tb)/(dv->get_Rmax_LTAR()); // *2.0phyllochron for senesence, set longer (twice) than appearance
	//DT - Nov 20, 2012 modified senescence to correspond with expansion - use same optimum and actual Ks 
	//DT took at *2 as a test
	double PotentialAreaDecrease=0;
	double Ks, Ks_Topt, ts;
	double Wk = totalLeaves/8.18;
	const double tt2 = phyllochron; //GDD sum from emegence to second tip appearance, guessed
	double tt = (rank-2)*phyllochron + tt2;
	double T_cur = dv->get_Tcur(); // get current temperature;
	double GDD_cur=dv->get_GDDsum(); //get current GDD DT- I think we need to make this relative
    double N_index = (2/(1+exp(-2.9*(N_content-0.25)))-1); //SK 8/20/10: as in Sinclair and Horie, 1989 Crop sciences, N availability index scaled between 0 and 1 based on 
	                                                       // This assumes 0.25mg/m2 minimum N required
	double N_stress = 0.1*(1.0-N_index); //SK 8/20/10: scaled between 0 and 0.1 to adjust the slope (Kx) bewteen 0.1 and 0.2
	double T_effect_size = max(0.0, (T_cur-Tb)/(T_peak-Tb)*exp(1.0-(T_cur-Tb)/(T_peak-Tb)));

	Ks_Topt = 0.01+(0.08+N_stress)*exp(-(pow((double)rank-1,2)/(2*pow(Wk,2)))); //this slope is slower than expansion by half under no N stress (0.1 compared to 0.2)
                                                       // the same as Eqn 8 in Lizaso et al. for growth 
	                                                   // note that ko=.02 and kx=0.2 (for expansion)
	                                                   // rank is leaf number. This is a 
	                                                   //slope parameter

	/*
	SK 8/22/10: The rate of senescence is regulated by N stress via N_index. When senescence is accelerated due to N stress, 
	this in turn reduces greenLeafArea and implicit redistribution of N takes place in plant.cpp when updating to calculate avg leaf N content which is a member
	of plant class currently. Eventually gas-exchange and all other balancing need to be done at individual nodal unit instead of the lumped method at the plant level
	as is being done currently.
	*/
	Ks = Ks_Topt; // taking out temperture for now *(dv->beta_fn(T_cur, 1.0, T_opt_ks, T_max_ks));  //See expand routine for explanation

    // ts is the senescent form of te used in expansion
	double T_effect_ts = 0.0, T_opt_ts;
	if (rank > 2) 
		{
		   const double Q10_ts=1.546;
		   T_opt_ts=T_opt_ks;
		   T_effect_ts= pow(Q10_ts, (T_cur-T_opt_ts)/10.0);
		   ts = tt + (2.197/Ks_Topt)*T_effect_ts; 
			
		}
    	else ts = 25*rank;
	double ll=get_longevity(); //need to add water and temperature effects here? will affect longevity?
	ts += get_longevity();
	PotentialAreaDecrease =max(0.0,FullyExpandedArea*Ks*exp(-Ks*(GDD_cur-ts))/pow(1+exp(-Ks*(GDD_cur-ts)),2)*dv->get_dGDD());
	                                           //potential leaf area increase without carbon limitation YY
	double dL=PotentialAreaDecrease;
	senescentArea+=dL;   
	senescentArea = min(FullyExpandedArea,senescentArea);
	if (senescentArea >= 0.01) aging = true;
	if (GDD_cur > (ts+ll/2.0)) terminated = true; // assume a leaf drops when GDD = ts + half longevity
	                //notice we increased longivity above
	//if (rank==4) cout << "senescent area for leaf 4->"<< senescentArea << " GA->" <<greenArea << " Area->" << area << "Fully" << FullyExpandedArea << endl;
	return;
}



void CLeaf::calcLongevity(double predawnlwp)
// see Lizaso et al. (2003)
{
	const double L0 = 150.0;
	const double Lx = 850.0;
	double LLy;  //longevity from previous time step
	double water_effect=LWPeffect(predawnlwp);
	double LN_l = 3.59 + 0.498*totalLeaves; // nodal position of most longevous leaf
	double Wi = (1.0/3.0)*totalLeaves; // Width function of the bell shape
	double LL =L0 + Lx * exp(-pow(rank-LN_l,2)/(2*pow(Wi,2)));
	LL=LL-L0*(1-water_effect);
	set_longevity(LL);
}

double CLeaf::GTI(double T_avg)
// general thermal index
// improved calculation of GDD, Steward et al. (1998)
{
	double b1 = 0.0432;
	double T_opt = 32.2;
	return b1*T_avg*T_avg*(1-0.6667*T_avg/T_opt);
}
CLeaf::~CLeaf() {}

//create a function which simulates the reducing in leaf expansion rate
//when predawn leaf water potential decreases. Parameterization of rf_psil
//and rf_sensitivity are done with the data from Boyer (1970) and Tanguilig et al (1987) YY
double CLeaf::LWPeffect(double predawn_psil)
{
    //DT Oct 10, 2012 changed this so it was not as sensitive to stress near -0.5 lwp
	double rf_psil=-1.0;
	double rf_sensitivity=0.5;
	double effect;
	effect=(1+exp(rf_psil*rf_sensitivity))/(1+exp(rf_sensitivity*(rf_psil-(predawn_psil+0.5))));
	if (effect >1 ) effect=1;
	return effect;
}


