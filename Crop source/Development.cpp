#include "stdafx.h"
#include <cmath>
#include "development.h"
#include <iostream>
#include <string>
#include <algorithm>
//#using <mscorlib.dll>
using namespace std;



CDevelopment::CDevelopment(const TInitInfo& info)
{    // DT -8/16/2016 reorganized some of these variables to make the initialization more consistent.
	 // removed SetParms method since it duplicated a much of what is done in the constructor. I moved
	 // the code from SetParms to here.

	/*note that info.timestep should have been updated - it seems to be pulling the 
		default value */
	dt = info.timeStep / MINUTESPERDAY; //converting minute to day decimal, 1= a day
	LvsAppeared = LvsExpanded = progressToAnthesis = 0; progressToTasselEmerg = 0;
	GerminationRate = EmergenceRate = LvsInitiated = 0;
	GDDsum = GDDgrain = dGDD = phyllochronsFromTI = 0;
	GTIsum = dGTI = 0;
	emerge_gdd = 0;
	stayGreen = info.stayGreen;
	LM_min = info.LM_min;
	Rmax_LIR = info.Rmax_LIR;
	Rmax_LTAR = info.Rmax_LTAR;
	DayLengthSensitive = info.DayLengthSensitive;
	PhyllochronsToSilk = info.PhyllochronsToSilk;
	PhyllochronsToTassel = info.PhyllochronsToTassel;
	totLeafNo = juvLeafNo = info.genericLeafNo;
	GDD_rating = info.GDD_rating;
	totLeafNo = juvLeafNo = info.genericLeafNo;
	Rmax_Germination = Rmax_Emergence = 0;

	Q10MR = info.Q10MR;
	Q10LeafSenescence = info.Q10LeafSenescense;
	WLRATIO = info.WLRATIO;
	A_LW = info.A_LW;
	leafNumberFactor_a1 = info.leafNumberFactor_a1;
	leafNumberFactor_a2 = info.leafNumberFactor_a2;
	leafNumberFactor_b1 = info.leafNumberFactor_b1;
	leafNumberFactor_b2 = info.leafNumberFactor_b2;




	initLeafNo = youngestLeaf = 5;
	curLeafNo = 1;
	LvsAtTI = 1;
	LvsInitiated = initLeafNo;

	LvsToInduce = 0.0;
	inductionPeriod = 0.0;
	inductions = 0;

	T_grow_sum = steps = 0.0;
	T_grow = T_ind = -99;

	shadeEffect = 1.0;

	Rmax_Germination = 0.45*dt; // max rate of germination per day, assume it takes two day at 31 C, needs to be replaced
								// Itabari et al., 1993 Expl Agric. 29:351-364 has info on this
	Rmax_Emergence = 0.2388*dt;
	Rmax_LTAR = Rmax_LTAR*dt; // Kim et al. (2007); Kim and Reddy (2004), 0.581 from Yan and Hunt (1999), equivalent phyllochron in CERES
							  //cdt changed from 0.524 to test for colorado data used 0.374
	Rmax_LIR = Rmax_LIR*dt; // best fit of K and W (1983), Kim and Reddy (2004) originally 0.978 (for colorado tried 0.558
	T_base = info.T_base;
	T_opt =  info.T_opt; // These Topt and Tceil values from Kim et al. (2007), also see Kim and Reddy (2004), Yan and Hunt (1999), SK
	T_ceil = info.T_ceil;
	T_opt_GDD = info.T_opt_GDD;
	

	P2 = 0.5;
	// Save for later
	initInfo = info;

}

CDevelopment::~CDevelopment(void)
{
}


int CDevelopment::update(const TWeather& wthr)
{
	double Jday = wthr.jday;
	T_cur = max(0., wthr.airT);
	T_air = wthr.airT;
	if (LvsAppeared < 9) T_cur = max(0., wthr.soilT);
	double addLeafPhotoPeriod, addLeafTemperature, addLeafTotal;
	//	double dt = initInfo.timeStep/(24*60); //converting minute to day decimal, 1= a day

	if (!germination.done)
	{
		//TODO: implement germination rate model of temperature.
		// for now assume it germinates immidiately after sowing
		GerminationRate += beta_fn(T_cur, Rmax_Germination, T_opt, T_ceil);
		if (GerminationRate >= 0.5)
		{
			germination.done = true;
			germination.daytime = wthr.daytime;
			// initialize T_grow
			T_grow = T_cur;
			T_grow_sum = T_cur;

			cout << " Germinated: GDDsum " << GDDsum << " time step (min): " << dt*(24 * 60) << endl;
		}
	}
	else // if (germination.done)
	{
		T_grow_sum += T_cur;
		steps++;
		T_grow = T_grow_sum / steps; // mean growing season temperature since germination, SK 1-19-12

		if (!emergence.done)
		{
			EmergenceRate += beta_fn(T_cur, Rmax_Emergence, T_opt, T_ceil);
			if (EmergenceRate >= 1.0)
			{
				emergence.done = true;
				emergence.daytime = wthr.daytime;
				
				cout << " Emergence: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
				GDDsum = 0.0; //reset GDDsum from emergernce, SK
				emerge_gdd = GDDsum; //gdd at emergence YY 4/2/09
				//	corn->LvsAppeared = 1.0;
				//LvsAppeared = 2;
			}
		}
		if (!tasselInitiation.done &&  emergence.done)
		{
			double temp1 = LvsInitiated;
			LvsInitiated += beta_fn(T_cur, Rmax_LIR, T_opt, T_ceil);
			if (LvsInitiated < temp1)
			{
				int iii = 1;
			}
			curLeafNo = (int)LvsInitiated;
			if (LvsInitiated >= juvLeafNo)
				// inductive phase begins after juvenile stage and ends with tassel initiation
			{
				//Equation 4 in Grant 1989 Ag. J. (81)
			   //dt 12/11/2012 broke the equation in two to separate temperature and daylenght effects
				if (T_ind == -99) { T_ind = T_grow; } //mean temperature during induction period
				addLeafTemperature = __max(0.0, (13.6 - 1.89*T_ind + 0.081*T_ind*T_ind - 0.001*T_ind*T_ind*T_ind));
				addLeafPhotoPeriod = 0.0;
				if (DayLengthSensitive)
				{
					addLeafPhotoPeriod = __max(0.0, 0.1*(juvLeafNo - 10.0)*(wthr.dayLength - 12.5));
				}
				addLeafTotal = (addLeafTemperature + addLeafPhotoPeriod);

				//addLeaf = __max(0, 0.1*(juvLeafNo-10.0)*(wthr.dayLength-12.5) + (13.9-1.89*T_cur+0.0795*T_cur*T_cur - 0.001*T_cur*T_cur*T_cur)); //Equation 4 in Grant 1989 Ag. J. (81)

				// effect of photoperiod and temperature on leaf no. used as Grant (1989)
				// Added back the temperature effect on leaf number and revised the algorithm to accumulate addLeafNo to totLeafNo.
				// Changed to respond to mean growing season temperature upto this point. 
				// This has little mechanistic basis. Needs improvements. SK 1-19-12
				LvsToInduce = (LvsToInduce*inductions + addLeafTotal) / (inductions + 1);
				T_ind = (T_ind*inductions + T_cur) / (inductions + 1);
				inductions++;
				inductionPeriod += dt;
				//	totLeafNo = juvLeafNo + addedLvs/inductionPeriod; //get a mean value over this period
				//	LvsAtTI = LvsInitiated; //Should be LvsInitiated. Already confirmed with Soo. 7/27/2006
					// uncomment the following for debugging
				 //   cout << " Inductive phase: " << LvsInitiated << " " << totLeafNo << " " << juvLeafNo << " " << addedLvs/inductionPeriod << endl;
				double actualAddedLvs = LvsInitiated - juvLeafNo;
				if (actualAddedLvs >= LvsToInduce)
				{
					youngestLeaf = totLeafNo = (int)LvsInitiated;
					curLeafNo = youngestLeaf;
					tasselInitiation.done = true;
					tasselInitiation.daytime = wthr.daytime;
					LvsInitiated = youngestLeaf;
					LvsAtTI = LvsAppeared;
					cout << " Tassel initiation: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
				}
			}
		}
		else if (tasselInitiation.done)
		{
			phyllochronsFromTI += beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil); // to be used for C partitoining time scaling, see Plant.cpp
		}

		if ((LvsAppeared < (int)LvsInitiated))
		{
			LvsAppeared += beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil);
			if (LvsAppeared >= (int)LvsInitiated)
			{
				LvsAppeared = LvsInitiated;
			}
		}
		//DT Sep 21, 2016 added a variable to indicate tasseling is done at 1 phyllocrhon from the last leaf appearing
		if (tasselInitiation.done && (LvsAppeared >= (int)LvsInitiated))
			//todo should move !silking.done to the if statement for silking
//		 if (((tasselInitiation.done) && (!silking.done) && !DayLengthSensitive)
//			 ||  ((LvsAppeared >= (int) LvsInitiated) && (!silking.done) && DayLengthSensitive))
		{
			if (!tasselFull.done)
				progressToTasselEmerg += beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil); // Assume full tassel emergence occurs at total tip appeared + 1.5 phyllochrons

			if ((progressToTasselEmerg >= PhyllochronsToTassel) && (!tasselFull.done)) //was 3
			{
				tasselFull.done = true;
				tasselFull.daytime = wthr.daytime;
				cout << " Tassel fully emerged: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
			}
			if (!silking.done && tasselFull.done)
				progressToAnthesis += beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil); // Assume 75% Silking occurs at total tip appeared + 3 phyllochrons

			if ((progressToAnthesis >= PhyllochronsToSilk) && !silking.done) //was 3
			{
				silking.done = true;
				silking.daytime = wthr.daytime;
				cout << " Silking: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
			}
		}
		if (silking.done)
		{
			GDDgrain += calcGDD(T_cur, T_opt_GDD)*dt;
			if (GDDgrain >= 170 && (!beginGrainFill.done)) // where is this number '170' from? SK
				//Todo: GTI was found more accurate for grain filling stage, See Thijs phenolog paper (2014)
			{
				beginGrainFill.done = true;
				beginGrainFill.daytime = wthr.daytime;
				cout << " Grain filling begins: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
			}
		}

	}

	//	if (!maturity.done)
	{
		dGTI = (calcGTI(T_cur, silking.done)*dt);
		dGDD = calcGDD(T_cur, T_opt_GDD)*dt;
		GDDsum += dGDD;
		GTIsum += dGTI;
//		if (GDDsum >= GDD_rating && (!maturity.done))
//		{
//			maturity.done = true;
//			maturity.daytime = wthr.daytime;
//			cout << " Matured: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
//		}

	}

	return 0;
}


double CDevelopment::beta_fn(double t, double R_max, double t_o, double t_c)
{
	// beta function, See Yin et al. (1995), Ag For Meteorol., Yan and Hunt (1999) AnnBot, SK
	double f, g, alpha;
	const double t_b = 0.0, beta = 1.0;
	if (t <= t_b || t >= t_c) return 0.0;
	if (t_c <= t_o || t_o <= t_b) return 0.0;

	f = (t - t_b) / (t_o - t_b);
	g = (t_c - t) / (t_c - t_o);
	alpha = beta*(t_o - t_b) / (t_c - t_o);

	return R_max*pow(f, alpha)*pow(g, beta);
}


double CDevelopment::calcGTI(double T_avg, bool Silked)
// General Thermal Index, Stewart et al. (1998)
 //Phenological temperature response of maize. Agron. J. 90: 73–79.
{
	//	double b1 = 0.0432;
	double b1 = 0.011178;
	double T_opt = 32.2;
	//return b1*T_avg*T_avg*(1-0.6667*T_avg/T_opt);
	//if (Silked = false) return b1*T_avg*T_avg*(1-0.6667*T_avg/T_opt);
	//else 
	return 5.358 + 0.011178*T_avg*T_avg;
}

double CDevelopment::calcGDD(double T_avg, double T_opt)
// GDD model with base 8. See Birch et al. (2003) Eu J Agron
{
	//double const T_base = 8.0;
	double const T_base = 8.0;
	//double const T_opt = 34.0;
	//double const T_opt = 30.0;
	return min(T_avg, T_opt) - T_base;
	//	if (Silked = false)  return b1*T_avg*T_avg*(1-0.6667*T_avg/T_opt);
	//	else return 5.358 + 0.011178*T_avg*T_avg;
}

//create a function which simulates the reducing in leaf expansion rate
//when predawn leaf water potential decreases. Parameterization of rf_psil
//and rf_sensitivity are done with the data from Boyer (1970) and Tanguilig et al (1987) YY
double CDevelopment::LWPeffect(double predawn_psi_bars, double threshold)
{
	//DT Oct 10, 2012 changed this so it was not as sensitive to stress near -0.5 lwp
	//SK Sept 16, 2014 recalibrated/rescaled parameter estimates in Yang's paper. The scale of Boyer data wasn't set correctly
	//sensitivity = 1.92, LeafWPhalf = -1.86, the sensitivity parameter may be raised by 0.3 to 0.5 to make it less sensitivy at high LWP, SK

	double psi_f = -1.4251; // -1.0, was -1.4251   later changed to -2.3;
	double s_f = 0.4258; // was 0.4258 0.5;
	double psi_th = threshold; // threshold wp below which stress effect shows up
	double effect;
	effect = __min(1.0, (1 + exp(psi_f*s_f)) / (1 + exp(s_f*(psi_f - (predawn_psi_bars - psi_th)))));
	if (effect >1) effect = 1;
	return effect;
}
double CDevelopment::LeafN_effect(double CriticalN)
{
	// calculate N effect on leaf growth
	double N_effect;
	N_effect = __max(0.0, (2 / (1 + exp(-2.9 * (CriticalN - 0.25))) - 1));
	N_effect = __min(1, N_effect);
	return N_effect;
}


