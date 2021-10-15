#pragma once
#ifndef _DEVELOPMENT_H_
#define _DEVELOPMENT_H_
#include "weather.h"
#include "initinfo.h"
#include <iostream>
#include <string>
using namespace std;

enum EPhase
  {
      Seed, Juvenile, Inductive, preSilking, Reproductive, Maturity
  };

struct TEvent
{
public:
	TEvent() {daytime = 0.0;  done=false;}
	double daytime;
	bool done;
};

class CDevelopment
{
public:
	CDevelopment(const TInitInfo&);
	~CDevelopment();
	void set_shadeEffect(double x) {shadeEffect=x;}
	double beta_fn(double t, double R_max, double t_m, double t_e);
	double calcGTI(double, bool);
	double calcGDD(double temperature,double T_opt);
	int update(const TWeather&); 
	double LWPeffect(double predawn_psi_bars, double threshold);
	double LeafN_effect(double CriticalN);
	TInitInfo get_initInfo() {return initInfo;}
	int get_youngestLeaf() {return youngestLeaf;}
	int get_totalLeaves() {return (int) totLeafNo;} // removed +1 that was previously here to count only those leaves fully initiated, SK 1-19-12
	//double get_totalLeaves() {return totLeafNo;} // take the value as double , SK 1-19-12
	double get_phyllochronsFromTI() {return phyllochronsFromTI;}
	double get_phyllochronsToSilk() { return PhyllochronsToSilk; }
	double get_progressToAnthesis() { return progressToAnthesis; }
	double get_PhyllochronsToTassel() { return PhyllochronsToTassel; }
	double get_stayGreen() { return stayGreen; }
	double get_LM_min() { return LM_min; }
	double get_LvsAtTI() {return LvsAtTI;}
	double get_LvsInitiated(){return LvsInitiated;} 
	double get_LvsAppeared(){return LvsAppeared;} // tip appreance rate is most conservative across cultivars and temperature regimes, persoanl comm with Dr. Tollenaar
	double get_GDDsum() {return GDDsum;}
	double get_EmergeGdd() {return emerge_gdd;}
	double get_dGDD() {return dGDD;}
	double get_dt() {return dt;} // added 08-16-11, SK
	double get_Tcur() {return T_cur;}
	double get_Tavg() {return T_avg;}
	double get_Rmax_LIR() {return Rmax_LIR/dt;}
	double get_Rmax_LTAR() {return Rmax_LTAR/dt;}
	double get_T_Opt() {return T_opt;}
	double get_Tbase() {return T_base;}
	double get_T_ceil() {return T_ceil;}
	double get_Tgrow()  {return T_grow;}
    double get_shadeEffect()  {return shadeEffect;}
	double get_T_Air() { return T_air; }
	double get_Q10LeafSenescence() { return Q10LeafSenescence; }
	double get_Q10MR() { return Q10MR; }
	double get_WLRATIO() { return WLRATIO; }
	double get_A_WL() { return A_LW; }
	double get_leafNumberFactor_a1() { return leafNumberFactor_a1; }
	double get_leafNumberFactor_a2() { return leafNumberFactor_a2; }
	double get_leafNumberFactor_b1() { return leafNumberFactor_b1; }
	double get_leafNumberFactor_b2() { return leafNumberFactor_b2; }



	//don't think we need the set functions
	void set_Q10LeafSenescence(double x) { Q10LeafSenescence = x; }
	void set_Q10MR(double x) { Q10MR = x; }
	void set_WLRATIO(double x) { WLRATIO = x; }
	void set_A_LW(double x) { A_LW = x; }
	void set_leafNumberFactor_a1(double x) {leafNumberFactor_a1 = x; }
	void set_leafNumberFactor_a2(double x) { leafNumberFactor_a2 = x; }
	void set_leafNumberFactor_b1(double x) { leafNumberFactor_b1 = x; }
	void set_leafNumberFactor_b2(double x) { leafNumberFactor_b2 = x; }


	bool Germinated() {return germination.done;}
	bool Emerged() {return emergence.done;}
	bool TasselInitiated() {return tasselInitiation.done;}
	bool Tasseled() { return tasselFull.done;}
	bool Flowered() {return anthesis.done;}
	bool Silked() {return silking.done;}
	bool GrainFillBegan() {return beginGrainFill.done;}
	bool Matured() {return maturity.done;}
	bool Dead() {return death.done;} // when all leaves are senescend and dead, the whole-plant is pronounced dead. SK

	string getNote() {return note;}
	TEvent germination;
	TEvent emergence;
	TEvent tasselInitiation;
	TEvent tasselFull;
	TEvent anthesis;
	TEvent silking;
	TEvent beginGrainFill;
	TEvent maturity;
	TEvent death;

private:
	CDevelopment(const CDevelopment&); // supressing shallow copy constructor
	int GDD_rating;
	double dGDD; // delta GDD
	double T_opt_GDD; //optimum temperature for GDD calcs
	double dt, steps; // timestep, step counter, SK
	double GDDsum; // cumulative GDD from sowing with Tbase= 8.0 and Topt = 34 as in CERES
	double emerge_gdd; //thermal time need for a corn plant to emerge YY 4/2/09
	double GDDgrain; // cumulative GDD from silking
	double dGTI; // delta GTI
	double GTIsum; // cumulative GTI, Stewart et al (1998), equivalent to GDD of Tbase = 10 and Topt 30
	double Rmax_LIR, Rmax_LTAR, Rmax_Germination, Rmax_Emergence, stayGreen, LM_min;
	bool DayLengthSensitive; //True if Day Length Sensitive
	double T_base, T_opt, T_ceil, T_cur, T_avg, T_grow, T_grow_sum, T_ind, T_air; // T_grow: mean temperature of the growing season from day 1 up to now, SK
	double totLeafNo, LvsToInduce, juvLeafNo, LvsAtTI, phyllochronsFromTI; //number of total, juvenile (genetic coeff) lvs, and lvs appeared at tassel initiation
	double P2; //photoperiod sensitivity as used in CERES-Maize
	double GerminationRate, EmergenceRate, LvsInitiated, LvsAppeared, LvsExpanded;
	double progressToAnthesis, progressToTasselEmerg, inductionPeriod;
	int initLeafNo,  youngestLeaf, curLeafNo, inductions; 
	double PhyllochronsToSilk; // number of phyllochrons past tassel initiation when silking takes place
	double PhyllochronsToTassel; // number of phyllochrons past tassel initiation when tassels are fully emerged
	double shadeEffect; //effect of shade with respect to R/FR, see JXB (2014) Zhu, 8-13-2015, SK
	double Q10LeafSenescence, Q10MR; //Q10 for respiration and leaf sensecense
    double WLRATIO; // leaf lamina width to length ratio
	double A_LW; // leaf area coeff with respect to L*W
	// coefficients in equation to calculate area of leaf based on largest leaf from Birch et al.,
	double leafNumberFactor_a1, leafNumberFactor_a2, leafNumberFactor_b1, leafNumberFactor_b2;


	string note;
	TInitInfo initInfo;
	
};
#endif