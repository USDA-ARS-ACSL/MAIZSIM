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
	length=width = area = PotentialArea = plastochrons = GDD2mature = 0.0;
	SLA = 200.0; // temporary for now - it should vary by age. Value comes from some of Soo's work
	                
	PotentialAreaIncrease = 0;
	RelativeAreaIncrease=0;
	FullyExpandedArea=0;
	TotalDroppedLeaves=0;
	actualArea = actualgreenArea=greenArea = senescentArea = droppedArea = unstressedArea=0;
	actualLength=actualwidth=0;
	initiated = appeared = growing = mature = aging = dead = dropped = false;
	phase1Delay = growthDuration=stayGreenDuration = seneDuration = 0.0;
	ptnLength=ptnWidth = 0.0;
	elongRate = 0.564;  //0.564 cm dd-1 Fournier and Andrieu 1998 Pg239. This is the "potential" elongation rate with no water stress Yang
	maxElongRate = 12.0; //max elongation rate (cm per day) at optipmal temperature (Topt: 31C with Tbase = 9.8C using 0.564 cm/dd rate from Fournier 1998 paper above
	elongAge = 0.0; //physiological age during expansion phase
	seneAge = 0.0; //age during senscence phase
	activeAge = 0.0; // age during active phase after fully expansion before senescence
	old_leaf = 0; //record leaf area in the last time step
    N_content = 3.0; //no N stress 
	WLRATIO = 0.106; // leaf lamina width to length ratio
	A_LW = 0.75; // leaf area coeff with respect to L*W
	stayGreen = 4.5; // staygreen trait of the hybrid
                                  // stay green for this value times growth period after peaking before senescence begins
                                  // An analogy for this is that with no other stresses involved, it takes 15 years to grow up, stays active for 60 years, and age the last 15 year if it were for a 90 year life span creature.
	                              //Once fully grown, the clock works differently so that the hotter it is quicker it ages
}

CLeaf::~CLeaf() {}


void CLeaf::initialize (CDevelopment * dv) 
// set potential leaf area of the current rank based on generic leaf no, see Fournier and Andrieu (1998), Birch et al., (1998)
// Routine below to calculate potentialArea will repeat in CLeaf::update to get actual leaf area after the inductive phase to any additional leaves developed before tassel initiation, SK 1-20-12
{
		COrgan::initialize();
		calc_dimensions(dv);
		initiated = true;
	//	growing=true; // DT
		first=true;
		stayGreenDuration = stayGreen*growthDuration;
		seneDuration = growthDuration; 
	// uncomment for debugging
	#if _DEBUG
			  std::cout << " GDDay " << dv->get_GDDsum()  << " leaf " << dv->get_LvsInitiated() << " " << rank << " totalLeaves " << dv->get_totalLeaves() 
				  <<  " potential area " << PotentialArea << std::endl;
	#endif
}

void CLeaf::calc_dimensions(CDevelopment *dv)
{
	const double LM_min = 115, k = 24.0;
	double totalLeaves = dv->get_totalLeaves(); //todo: should be a plant parameter not leaf
	double L_max = (sqrt(LM_min*LM_min + k*(totalLeaves - dv->get_initInfo().genericLeafNo)));
	                                     //LM_min is a length characteristic of the longest leaf,in Fournier and Andrieu 1998, it was 90 cm
	                                     //LA_max is a fn of leaf no (Birch et al, 1998 fig 4) with largest reported value near 1000cm2. This is implemented as lfno_effect below, SK
	                                     //LM_min of 115cm gives LA of largest leaf 1050cm2 when totalLeaves are 25 and Nt=Ng, SK 1-20-12
	                                     // Without lfno_effect, it can be set to 97cm for the largest leaf area to be at 750 cm2 with Nt ~= Ng (Lmax*Wmax*0.75) based on Muchow, Sinclair, & Bennet (1990), SK 1-18-2012
	                                     // Eventually, this needs to be a cultivar parameter and included in input file, SK 1-18-12
	                                     // the unit of k is cm^2 (Fournier and Andrieu 1998 Pg239). YY
  										 // L_max is the length of the largest leaf when grown at T_peak. Here we assume LM_min is determined at growing Topt with minmal (generic) leaf no, SK 8/2011
                                         // If this routine runs before TI, totalLeaves = genericLeafNo, and needs to be run with each update until TI and total leaves are finalized, SK
	double n_m, a, b;

	n_m = 5.93 + 0.33*totalLeaves; // the rank of the largest leaf. YY (need to adjust here for the maximum leaf size as a function of plant pop -DT 8/12/2015)
	a = -10.61 + 0.25*totalLeaves;
	b = -5.99 + 0.27*totalLeaves;
    //equation 7 in Fournier and Andrieu (1998). YY

	ptnLength = L_max*exp(a/2*pow(rank/n_m-1,2)+b/2*pow(rank/n_m-1,3)); //*dv->get_shadeEffect()
	                                                                    //equa 8(b)(Actually eqn 6? - eqn 8 deals with leaf age - DT)
	                                                                   //in Fournier and Andrieu(1998). YY
	growthDuration = ptnLength/maxElongRate; // shortest possible linear phase duration in physiological time (days instead of GDD) modeified form of equa 8(a)Fournier and Andrieu(1998)
    phase1Delay = __max(0.0, -5.16+1.94*rank); //not used in MAIZSIM because LTAR is used to initiate leaf growth. Fournier's value : -5.16+1.94*rank;equa 11 Fournier and Andrieu(1998) YY, This is in plastochron unit

	double W_max = L_max*WLRATIO;//Fournier and Andrieu(1998) Pg242 YY
	double LA_max = L_max*W_max*A_LW; // daughtry and hollinger (1984) Fournier and Andrieu(1998) Pg242 YY
    double lfno_effect = max(0.5, min(1.0, exp(-1.17+0.047*totalLeaves))); // Fig 4 of Birch et al. (1998) 

	{
		PotentialArea=lfno_effect*LA_max*exp(a*pow(rank/n_m-1,2)+b*pow(rank/n_m-1,3)); //equa 6. Fournier and Andrieu(1998) multiplied by Birch et al. (1998) leaf no effect
	                                                               //LA_max the area of the largest leaf
	}                                                          //PotentialArea potential final area of a leaf with rank "n". YY

}

void CLeaf::update(CDevelopment * dv, double PredawnLWP)
{ 
	COrgan::set_temperature(dv->get_Tcur());
	COrgan::update();
//	calc_dimensions(dv);
	expand(dv, PredawnLWP);
	senescence(dv, PredawnLWP);
    greenArea = max(0.0, area-senescentArea);
}

void CLeaf::expand(CDevelopment * dv, double PredawnLWP)
//leaf expansiopn rate based on a determinate sigmoid function by Yin et al. (2003)
{
    double CriticalNitrogen;
	CriticalNitrogen= __max(0.25,this->N_content);
	double N_effect=(2/(1+exp(-2.9*(CriticalNitrogen-0.25)))-1);
	const double T_peak = 18.7, Tb = 8.0, phyllochron = (dv->get_T_Opt()- Tb)/(dv->get_Rmax_LTAR()); 
    // T_peak is the optimal growth temperature at which the potential leaf size determined in calc_mophology achieved. Similar concept to fig 3 of Fournier and Andreiu (1998) 
	// phyllochron corresponds to PHY in Lizaso (2003); phyllochron needed for next leaf appearance in degree days (GDD8) - 08/16/11, SK. 

	double T = dv->get_Tcur();
	double T_gro = dv->get_Tgrow();
	double T_effect_size = max(0.0, (T_gro-Tb)/(T_peak-Tb)*exp(1.0-(T_gro-Tb)/(T_peak-Tb))); 
	//final leaf size is adjusted by growth temperature determining cell size during elongation
    //final leaf 
	// See Kim et al. (2012) Agro J. for more information on how this relationship has been derermined basned o multiple studies and is applicable across environments
    const double psi_threshold_bars = -0.8657; 
	double water_effect=LWPeffect(PredawnLWP, psi_threshold_bars);
	double dD = dv->get_initInfo().timeStep/MINUTESPERDAY; // time step as day fraction

	if (dv->get_LvsAppeared() >= rank && !appeared) 
	{
		appeared = true;
	}

	double t_e = growthDuration; // end of growth period, time to maturity
	double t_m = t_e/2; // max. growth rate assumed to be half of t_e
  //+(1-C2_effect) proposed C2_effect for elongAge, add to 1.0
	if (appeared && !mature)
	{
		elongAge += dv->beta_fn(T, 1.0,dv->get_T_Opt(), dv->get_T_ceil())*dD; // Todo: implement Parent and Tardieu (2011, 2012) approach for leaf elongation in response to T and VPD, and normalized at 20C, SK, Nov 2012
		// elongAge indicates where it is now along the elongation stage or duration. duration is determined by totallengh/maxElongRate which gives the shortest duration to reach full elongation in the unit of days.
		elongAge = __min(t_e, elongAge);

//		area = __max(0.0, water_effect*T_effect_size*PotentialArea*(1.0 + (t_e-elongAge)/(t_e-t_m))*pow(elongAge/t_e, (t_e/(t_e-t_m))));
		double maxExpansionRate = T_effect_size*PotentialArea*(2*t_e-t_m)/(t_e*(t_e-t_m))*pow(t_m/t_e,t_m/(t_e-t_m));
		PotentialAreaIncrease =__max(0.0,maxExpansionRate*__max(0.0, (t_e-elongAge)/(t_e-t_m)*pow(elongAge/t_m,t_m/(t_e-t_m)))*dD);
	                                           //potential leaf area increase without any limitations
      //C2_effect = 1.0; // place holder
		double dA = PotentialAreaIncrease; // growth temperature effect is included in determining potential area
		water_effect=1.0;
		N_effect=1.0;
		area += dA*__min(water_effect,N_effect)*dv->get_shadeEffect();
		unstressedArea+=dA;
		if (unstressedArea >= PotentialArea || elongAge >= t_e) 
		{
			mature = true;
			set_GDD2mature (get_physAge());
			growing = false;
		}
		else growing = true;

	}
	return;
}


void CLeaf::senescence(CDevelopment * dv, double PredawnLWP)
{
	double dD = dv->get_initInfo().timeStep/MINUTESPERDAY;
	double T = dv->get_Tcur();
	double T_opt = dv->get_T_Opt();
	double T_grow = dv->get_Tgrow();
    double N_index = __max(0.0, (2/(1+exp(-2.9*(N_content-0.25)))-1)); //SK 8/20/10: as in Sinclair and Horie, 1989 Crop sciences, N availability index scaled between 0 and 1 based on 
	// This assumes 0.25mg/m2 minimum N required, and below this the value is 0.0.
	const double psi_threshold_bars = -4.0; //threshold predawn leaf water potential (in bars) below which water stress triggers senescence, needs to be substantiated with lit or exp evidence, SK
	// This is the water potential at which considerable reduction in leaf growth takes place in corn, sunflower, and soybean in Boyear (1970)
	double water_effect=LWPeffect(PredawnLWP, psi_threshold_bars);

	double t_e, t_m; // max. growth rate assumed to be half of t_e

	double Q10 = 2.0;
	double q10fn = pow(Q10,(T - T_opt)/10);
        // Assumes physiological time for senescence is the same as that for growth though this may be adjusted by stayGreen trait
		// a peaked fn like beta fn not used here because aging should accelerate with increasing T not slowing down at very high T like growth,
		// instead a q10 fn normalized to be 1 at T_opt is used, this means above Top aging accelerates. 

	double scale = 1.0; // scale for reduction in leaf lifespan and aging rate

	if (!mature && !aging && !dead)
	{
		stayGreenDuration = stayGreen*growthDuration;
		seneDuration = growthDuration; // end of growth period, time to maturity
	}
	else if (mature && !aging && !dead)
	{
		activeAge += q10fn*dD;
		scale = 0.5;
		stayGreenDuration = __max(0.0, stayGreenDuration - scale*(1.0 - water_effect)*dD);
		//One day of cumulative severe water stress (i.e., water_effect = 0.0 around -4MPa) would result in a reduction of leaf lifespan in relation staygreeness and growthDuration, SK
        //if scale is 1.0, one day of severe water stress shortens one day of stayGreenDuration
		if (activeAge >= stayGreenDuration)
		{
			activeAge = stayGreenDuration;
			aging = true;
		}
	}	
	else if (aging && !dead)
	{
  		seneAge += q10fn*dD; 
		scale = 0.5;
        //if scale is 0.5, one day of severe water stress at predawn shortens one half day of agingDuration
		seneDuration = __max(0.0, seneDuration - scale*(1.0 - water_effect)*dD);
		t_e = seneDuration;
		t_m = t_e/2;
		seneAge = __min(t_e, seneAge);
//		double maxRate = area*(2*t_e-t_m)/(t_e*(t_e-t_m))*pow(t_m/t_e,t_m/(t_e-t_m));
//		double nonstressAgingRate = __max(0.0,maxRate*__max(0.0, (t_e-seneAge)/(t_e-t_m)*pow(seneAge/t_m,t_m/(t_e-t_m))));
// 		double dA = nonstressAgingRate*dD;
//Leaf senescence accelerates with drought and heat. see http://www.agry.purdue.edu/ext/corn/news/timeless/TopLeafDeath.html
//		senescentArea += dA;
		senescentArea = __max(0.0, this->area*(1.0 + (t_e-seneAge)/(t_e-t_m))*pow(seneAge/t_e, (t_e/(t_e-t_m))));
		if (senescentArea >= area || seneAge >= t_e) 
		{
			senescentArea = area; dead = true;
		}
		else dead = false; 
	}
	else if (dead && get_physAge() >= get_GDD2mature())
	{

		dropped = true;
	}

	return;
}


//create a function which simulates the reducing in leaf expansion rate
//when predawn leaf water potential decreases. Parameterization of rf_psil
//and rf_sensitivity are done with the data from Boyer (1970) and Tanguilig et al (1987) YY
double CLeaf::LWPeffect(double predawn_psi_bars, double threshold)
{
    //DT Oct 10, 2012 changed this so it was not as sensitive to stress near -0.5 lwp
	//SK Sept 16, 2014 recalibrated/rescaled parameter estimates in Yang's paper. The scale of Boyer data wasn't set correctly
	//sensitivity = 1.92, LeafWPhalf = -1.86, the sensitivity parameter may be raised by 0.3 to 0.5 to make it less sensitivy at high LWP, SK

	double psi_f=-1.4251; // -1.0;
	double s_f=0.4258; // 0.5;
	double psi_th = threshold; // threshold wp below which stress effect shows up
	double effect;
	effect=__min(1.0, (1+exp(psi_f*s_f))/(1+exp(s_f*(psi_f-(predawn_psi_bars-psi_th)))));
	if (effect >1 ) effect=1;
	return effect;
}


