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
	bool isGrowing() {return growing;}
	bool isProlific() {return prolific;}
	bool isAging() {return aging;}
	bool isTerminated() {return terminated;}
	bool isDropped() {return dropped;}
	double get_area(){return area;}
	double get_greenArea() {return greenArea;}
	double get_senescentArea() {return senescentArea;}
	double get_potentialArea(){return PotentialArea;}
	double get_length() {return length;}
	double get_SLA() {return SLA;}
	double get_N_content() {return N_content;}
	double get_potentialAreaIncrease () {return PotentialAreaIncrease;}
	double get_RelativeAreaIncrease() {return RelativeAreaIncrease;}
	double get_actualgreenArea() {return actualgreenArea;}
	double get_droppedArea () {return droppedArea;}
	int    get_TotalGrowingLeaves() {return TotalGrowingLeaves;}
	int    get_TotalDroppedLeaves() {return TotalDroppedLeaves;}
	int    get_Rank() {return rank;}
	

	void initialize(CDevelopment * dv);
	void set_area(double x) {area=x;}
	void set_length(double x) {length=x;}
	void set_SLA(double x) {SLA=x;}
	void set_TotalGrowingLeaves(int x) {TotalGrowingLeaves=x;}
	void set_TotalDroppedLeaves(int x) {TotalDroppedLeaves=x;}
	double GTI(double);
	void set_RelativeAreaIncrease(double x) {RelativeAreaIncrease=x;}
	void update(CDevelopment *, double pdlwp);
	void Expand(CDevelopment *, double pdlwp);
	void senescence(CDevelopment *, double pdlwp);
	void set_N_content(double x) {N_content=x;}

	void calcLongevity(double pdlwp);
	 double LWPeffect(double predawn_psil);

private:
	CLeaf(const CLeaf&);
	bool initiated, growing, prolific, aging, terminated, dropped;
	int rank;
	int totalLeaves; // potential total leaf number?
	int TotalGrowingLeaves, TotalDroppedLeaves;
	double PotentialArea; // potential leaf area
	double FullyExpandedArea;
	double PotentialAreaIncrease, old_leaf;//potential leaf area increase without temperature or water limitation YY
	double RelativeAreaIncrease; //Area increase of this leaf relative to all other leaves. Used to partition carbon
	double area; // actual leaf area
	double greenArea, actualgreenArea;//actualgreenArea is the green area of leaf growing under carbon limitation
	//SK 8/22/10: There appears to be no distinction between these two variables in the code.
	double senescentArea, droppedArea;
	double length;
	double width;
	double SLA;
	double plastochrons; // Fournier and Andrieu (1998), used for delay between initiation and elongation
	double phase1Delay, phase2Duration, elongAge, elongRate;
	double ptnLength, ptnWidth;
	double LeafCalibTemperature; //temperature at which experiments run where parameters for leaf expansion were determined
	double actualArea;
	double actualLength,actualwidth; //actual length and width of leaf under both drought stress and carbon limitation
	bool   first;  //indicates if this is the first time we call the elongate method;
	double N_content; 
	
	/*
	SK 8/22/10: Leaf N content in mg/m2
	*/

};
#endif