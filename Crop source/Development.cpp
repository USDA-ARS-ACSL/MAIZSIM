#include "StdAfx.h"
#include <cmath>
#include "development.h"
#include "radiation.h"
#include <iostream>
#include <string>
#include <algorithm>
//#using <mscorlib.dll>
using namespace std;



CDevelopment::CDevelopment(const TInitInfo& info)
{
	LvsInitiated = LvsAppeared = LvsExpanded = Anthesis=0;
    GerminationRate = EmergenceRate = LvsInitiated = LvsAppeared = LvsExpanded = Anthesis =0;
	GDDsum = GDDgrain = dGDD= P2 = phyllochronsFromTI =0;
	GTIsum = dGTI = 0;
	GDD_rating = 1000;
	emerge_gdd = 0;
	Rmax_LIR =info.Rmax_LIR;
    Rmax_LTAR = info.Rmax_LTAR;
	DayLengthSensitive=info.DayLengthSensitive;
	Rmax_Germination = Rmax_Emergence =0;
	T_base = 8.0;  T_opt = 30.0; T_ceil = 40.0; 
	totLeafNo = juvLeafNo = info.genericLeafNo; 
	initLeafNo =  youngestLeaf = 5;
	curLeafNo =1; 
	LvsAtTI = 1;
	initInfo = info;
	LvsToInduce = 0.0;
	inductionPeriod = 0.0;
	inductions = 0;
	dt = initInfo.timeStep/MINUTESPERDAY; //converting minute to day decimal, 1= a day
	T_grow_sum = steps = 0.0;
    T_grow =  T_ind = -99;
	PhyllochronsToSilk=info.PhyllochronsToSilk;
	shadeEffect=1.0;
	setParms();
}

CDevelopment::~CDevelopment(void)
{
}

void CDevelopment::setParms() // dt in days
{
	initLeafNo = 5;
	totLeafNo = juvLeafNo=initInfo.genericLeafNo;
	Rmax_Germination = 0.45*dt; // max rate of germination per day, assume it takes two day at 31 C, needs to be replaced
	Rmax_Emergence = 0.2388*dt ;
	Rmax_LTAR = Rmax_LTAR*dt; // Kim et al. (2007); Kim and Reddy (2004), 0.581 from Yan and Hunt (1999), equivalent phyllochron in CERES
	                //cdt changed from 0.524 to test for colorado data used 0.374
	Rmax_LIR = Rmax_LIR*dt; // best fit of K and W (1983), Kim and Reddy (2004) originally 0.978 (for colorado tried 0.558
	T_base = 8.0;
	T_opt  = 32.1; // These Topt and Tceil values from Kim et al. (2007), also see Kim and Reddy (2004), Yan and Hunt (1999), SK
	T_ceil = 43.7;
	LvsInitiated = initLeafNo;
	GDD_rating = initInfo.GDD_rating;
	P2 = 0.5;
}


int CDevelopment::update(const TWeather& wthr)
{
	double Jday = wthr.jday;
    T_cur = max(0.,wthr.airT);
	if (LvsAppeared < 9) T_cur = max(0.,wthr.soilT); 
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
			cout << "* Germinated: GDDsum " << GDDsum << " time step (min): " << dt*(24*60) << endl;
		}
	}
	else // if (germination.done)
	{
		T_grow_sum += T_cur; 
		steps ++;
		T_grow = T_grow_sum/steps; // mean growing season temperature since germination, SK 1-19-12

		if(!emergence.done)
		{
		    EmergenceRate += beta_fn(T_cur, Rmax_Emergence, T_opt, T_ceil);
			// corn->LvsAppeared = corn->EmergenceRate;
			if(EmergenceRate >= 1.0)
			{
				emergence.done = true;
				emergence.daytime = wthr.daytime;
				GDDsum = 0.0; //reset GDDsum from emergernce, SK
				cout << "* Emergence: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
				emerge_gdd = GDDsum; //gdd at emergence YY 4/2/09
				//	corn->LvsAppeared = 1.0;
			}
		}
		if (!tasselInitiation.done)
		{
			LvsInitiated += beta_fn(T_cur, Rmax_LIR, T_opt, T_ceil);
			curLeafNo = (int) LvsInitiated;
			if (LvsInitiated >= juvLeafNo)
			// inductive phase begins after juvenile stage and ends with tassel initiation
			{
                  //Equation 4 in Grant 1989 Ag. J. (81)
			     //dt 12/11/2012 broke the equation in two to separate temperature and daylenght effects
				if (T_ind == -99) {T_ind = T_grow;} //mean temperature during induction period
				addLeafTemperature = __max(0.0,(13.6-1.89*T_ind+0.081*T_ind*T_ind - 0.001*T_ind*T_ind*T_ind)); 
				addLeafPhotoPeriod=0.0;
				if (DayLengthSensitive)
				{
				  addLeafPhotoPeriod = __max(0.0, 0.1*(juvLeafNo-10.0)*(wthr.dayLength-12.5)); 
				}
				addLeafTotal=(addLeafTemperature + addLeafPhotoPeriod);
			
				//addLeaf = __max(0, 0.1*(juvLeafNo-10.0)*(wthr.dayLength-12.5) + (13.9-1.89*T_cur+0.0795*T_cur*T_cur - 0.001*T_cur*T_cur*T_cur)); //Equation 4 in Grant 1989 Ag. J. (81)
               
                // effect of photoperiod and temperature on leaf no. used as Grant (1989)
				// Added back the temperature effect on leaf number and revised the algorithm to accumulate addLeafNo to totLeafNo.
				// Changed to respond to mean growing season temperature upto this point. 
				// This has little mechanistic basis. Needs improvements. SK 1-19-12
				LvsToInduce = (LvsToInduce*inductions + addLeafTotal)/(inductions+1);
				T_ind = (T_ind*inductions + T_cur)/(inductions +1);
 				inductions ++;
				inductionPeriod += dt;
 			//	totLeafNo = juvLeafNo + addedLvs/inductionPeriod; //get a mean value over this period
			//	LvsAtTI = LvsInitiated; //Should be LvsInitiated. Already confirmed with Soo. 7/27/2006
				// uncomment the following for debugging
			 //   cout << "* Inductive phase: " << LvsInitiated << " " << totLeafNo << " " << juvLeafNo << " " << addedLvs/inductionPeriod << endl;
			    double actualAddedLvs = LvsInitiated - juvLeafNo;    
				if ( actualAddedLvs >= LvsToInduce)
				{
					youngestLeaf = totLeafNo = (int) LvsInitiated;
					curLeafNo = youngestLeaf;
					tasselInitiation.done =true;
					tasselInitiation.daytime = wthr.daytime;
					LvsInitiated = youngestLeaf;
					LvsAtTI = LvsAppeared;
					cout << "* Tassel initiation: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
				}
			}
		}
		else
		{
			phyllochronsFromTI += beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil); // to be used for C partitoining time scaling, see Plant.cpp
		}
			
		if ((LvsAppeared < (int) LvsInitiated))
		{ 
			LvsAppeared += beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil);
			if (LvsAppeared >= (int) LvsInitiated)
			{
                LvsAppeared = LvsInitiated;
			}
		}

		 if (tasselInitiation.done && (LvsAppeared >= (int) LvsInitiated) && !silking.done)
//		 if (((tasselInitiation.done) && (!silking.done) && !DayLengthSensitive)
//			 ||  ((LvsAppeared >= (int) LvsInitiated) && (!silking.done) && DayLengthSensitive))
		{
		    Anthesis += beta_fn(T_cur, Rmax_LTAR, T_opt, T_ceil); // Assume 75% Silking occurs at total tip appeared + 3 phyllochrons
			if (Anthesis >= PhyllochronsToSilk) //was 3
			{
                silking.done = true;
			    silking.daytime = wthr.daytime;
				cout << "* Silking: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
			}
		}
		if (silking.done)
		{
			GDDgrain += calcGDD(T_cur)*dt;
			if (GDDgrain >= 170 && (!beginGrainFill.done)) // where is this number '170' from? SK
				//Todo: GTI was found more accurate for grain filling stage, See Thijs phenolog paper (2014)
			{
				beginGrainFill.done = true;
			    beginGrainFill.daytime = wthr.daytime;
				cout << "* Grain filling begins: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
			}
		}

	}

//	if (!maturity.done)
	{
		dGTI = (calcGTI(T_cur, silking.done)*dt);
		dGDD = calcGDD(T_cur)*dt;
		GDDsum += dGDD;
		GTIsum += dGTI;
		if (GDDsum >= GDD_rating && (!maturity.done))
		{
			maturity.done = true;
			maturity.daytime = wthr.daytime;
			cout << "* Matured: GDDsum " << GDDsum << " Growing season T " << T_grow << endl;
		}

	}

	return 0;
}


double CDevelopment::beta_fn(double t, double R_max, double t_o, double t_c)
{
// beta function, See Yin et al. (1995), Ag For Meteorol., Yan and Hunt (1999) AnnBot, SK
	double f, g, alpha;
	const double t_b = 0.0, beta=1.0;
    if (t <= t_b || t >= t_c) return 0.0;
	if (t_c <= t_o || t_o <= t_b) return 0.0;

	    f = (t-t_b)/(t_o-t_b);
		g = (t_c-t)/ (t_c-t_o);
		alpha = beta*(t_o-t_b)/(t_c-t_o);

    return R_max*pow(f, alpha)*pow(g,beta);
}


double CDevelopment::calcGTI (double T_avg, bool Silked)
// General Thermal Index, Stewart et al. (1998)
 //Phenological temperature response of maize. Agron. J. 90: 73–79.
{
//	double b1 = 0.0432;
	double b1=0.011178;
	double T_opt = 32.2;
	//return b1*T_avg*T_avg*(1-0.6667*T_avg/T_opt);
	//if (Silked = false) return b1*T_avg*T_avg*(1-0.6667*T_avg/T_opt);
	//else 
	return 5.358 + 0.011178*T_avg*T_avg;
}

double CDevelopment::calcGDD (double T_avg)
// GDD model with base 8. See Birch et al. (2003) Eu J Agron
{
	//double const T_base = 8.0;
    double const T_base = 8.0;
	double const T_opt = 34.0;
	//double const T_opt = 30.0;
	return min(T_avg, T_opt)-T_base;
//	if (Silked = false)  return b1*T_avg*T_avg*(1-0.6667*T_avg/T_opt);
//	else return 5.358 + 0.011178*T_avg*T_avg;
}