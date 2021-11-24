#pragma once
#ifndef _LEAF_H_
#define _LEAF_H_
#include "organ.h"
#include "weather.h"
#include "development.h"

class CLeaf: public COrgan
//leaf blade, call it "leaf" for convenience
{
	friend class CSheath;
public:
	CLeaf(int rank, CDevelopment * dv); // take current leaf rank and total leaf number to calculate potentialArea
	~CLeaf();

	bool isInitiated() {return initiated;}
	bool isAppeared() {return appeared;}
	bool isGrowing() {return growing;}
	bool isMature() {return mature;}
	bool isAging() {return aging;}
	bool isDead() {return dead;}
	bool isDropped() {return dropped;}
	double get_area(){return area;}
	double get_greenArea() {return greenArea;}
	double get_senescentArea() {return senescentArea;}
	double get_potentialArea(){return PotentialArea;}
	double get_length() {return length;}
	double get_SLA() {return SLA;}
	double get_N_content() {return N_content;}
	double get_N_Effect() { return N_effect; }
	double get_CriticalN() { return CriticalNitrogen; }
	double get_potentialAreaIncrease () {return PotentialAreaIncrease;}
	double get_RelativeAreaIncrease() {return RelativeAreaIncrease;}
	double get_actualgreenArea() {return actualgreenArea;}
	double get_droppedArea () {return droppedArea;}
	double get_GDD2mature() {return GDD2mature;}
	double get_Elongation_Age() { return elongAge; }
	double get_psi_threshold_bars() { return psi_threshold_bars; }
	int    get_TotalGrowingLeaves() {return TotalGrowingLeaves;}
	int    get_TotalDroppedLeaves() {return TotalDroppedLeaves;}
	int    get_TotalMatureLeaves() { return TotalMatureLeaves; }
	int    get_Rank() {return rank;}
	
	

	void initialize(CDevelopment * dv);
	void set_area(double x) {area=x;}
	void set_length(double x) {length=x;}
	void set_SLA(double x) {SLA=x;}
	void set_TotalGrowingLeaves(int x) {TotalGrowingLeaves=x;}
	void set_TotalDroppedLeaves(int x) {TotalDroppedLeaves=x;}
	void  set_TotalMatureLeaves(int x) {TotalMatureLeaves=x; }

	
	

	double GTI(double);
	void set_RelativeAreaIncrease(double x) {RelativeAreaIncrease=x;}
	void update(CDevelopment *, double pdlwp);
	void expand(CDevelopment *, double pdlwp);
	void senescence(CDevelopment *, double pdlwp);
	void set_N_content(double x) {N_content=x;}
	void set_GDD2mature(double x) {GDD2mature=x;}

	void calc_dimensions(CDevelopment * dv);
	void calcLongevity(double pdlwp);
	double LWPeffect(double predawn_psil, double threshold);

private:
	CLeaf(const CLeaf&);
	bool initiated, appeared, growing, mature, aging, dead, dropped;
	int rank;
	int TotalGrowingLeaves, TotalDroppedLeaves, TotalMatureLeaves;
	double PotentialArea; // potential leaf area
	double FullyExpandedArea;
	double PotentialAreaIncrease, old_leaf;//potential leaf area increase without temperature or water limitation YY
	double RelativeAreaIncrease; //Area increase of this leaf relative to all other leaves. Used to partition carbon
	double area, unstressedArea; // actual leaf area and accumulated potential area
	double greenArea, actualgreenArea;//actualgreenArea is the green area of leaf growing under carbon limitation
	//SK 8/22/10: There appears to be no distinction between these two variables in the code.
	double senescentArea, droppedArea;
	double length;
	double width;
	double SLA;
	double plastochrons, GDD2mature; // Fournier and Andrieu (1998), used for delay between initiation and elongation
	double phase1Delay, growthDuration, elongAge, elongRate, maxElongRate, seneAge, activeAge, seneDuration;
	double stayGreen, stayGreenDuration;
	double ptnLength, ptnWidth,LM_min;
	double actualArea;
	double actualLength,actualwidth; //actual length and width of leaf under both drought stress and carbon limitation
	bool   first;  //indicates if this is the first time we call the elongate method;
	double N_content; 
	double N_effect;
	double WLRATIO, A_LW;
	double T_peak, Tb_Leaf;
	double Q10LeafSenescence;
	double leafNumberFactor_a1, leafNumberFactor_a2,
		   leafNumberFactor_b1, leafNumberFactor_b2;

	const double psi_threshold_bars = -0.8657;
	double CriticalNitrogen = __max(0.25,this->N_content);

	/*
	SK 8/22/10: Leaf N content in mg/m2
	*/

};
#endif