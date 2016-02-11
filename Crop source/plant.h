#pragma once
#ifndef _PLANT_H_
#define _PLANT_H_
#include "organ.h"
#include "nodalUnit.h"
#include "development.h"
#include "roots.h"
#include "ear.h"
#include "gas_exchange.h"
#include <iostream>
#include <string>

struct TStage
{
public:
	TStage() { V = 0.0; R = 0.0;}
	double V, R;
};

class CPlant
{
public:
	CPlant(const TInitInfo&);
	~CPlant();

	CNodalUnit* get_nodalUnit() {return nodalUnit;}
	CEar * get_ear() {return ear;}
	CRoots * get_roots() {return roots;}
	CDevelopment * get_develop() {return develop;}
	CGas_exchange * get_sunlit() {return this->sunlit;}
	CGas_exchange * get_shaded() {return this->shaded;}  //get access to the pointers that point to the sunlit/shade leaves Yang 8/29/06

	int get_nodeNumber() {return nodeNumber;}
	int get_finalNodeNumber() {return finalNodeNumber;}

	double get_mass() {return mass;}
	double get_age() {return age;}
	double get_CH2O() {return CH2O;}
	double get_N() {return TotalNitrogen;} //Total N in plant mg plant-1
	double get_Pg() {return photosynthesis_gross;}
	double get_Pn() {return photosynthesis_net;}
	double get_assimilate() {return assimilate;}
	double get_ET() {return transpiration;}
	double get_ET_Old() {return transpirationOld;}
	double get_tmpr() {return temperature;}
	double get_C_pool() {return C_pool;}
	double get_C_pool_root() {return C_pool_root;}
	double get_C_reserve() {return C_reserve;}
	double get_MaintenanceRespiration() {return maintRespiration;}
	double get_stemMass() {return stemMass;}
	double get_leafMass() {return leafMass;}
	double get_earMass() {return earMass;}
	double get_shootMass() {return shootMass;}
	double get_rootMass() {return rootMass;}
	double get_shootPart() {return shootPart;}
	double get_leafPart()  {return leafPart;}
	double get_rootPart() {return rootPart;}
	double get_DroppedLeafMass() {return droppedLeafmass;}
	double get_activeLeafMass() {return activeLeafMass;}
	double get_conductance() {return conductance;}
	double get_VPD() {return VPD;}
	double get_LeafArea() {return leafArea;}
	double get_LeafN(){return leaf_N;}
	double get_LeafNFraction () {return leaf_NFraction;}
	double get_HourlyNitrogenDemand () {return HourlyNitrogenDemand;}
	double get_CumulativeNitrogenDemand () {return CumulativeNitrogenDemand;}
	double get_HourlyNitrogenSoilUptake() {return HourlyNitrogenSoilUptake;}
	double get_CumulativeNitrogenSoilUptake() {return CumulativeNitrogenSoilUptake;}
	double get_droppedLfArea() { return currentDroppedLfArea;}
	double get_ptnLfIncrease() { return potentialLeafAreaIncrease;}
	double getSunlitLAI() {return SunlitLAI;}
	double getShadedLAI() {return ShadedLAI;}
	double getC2_effect() {return C2_effect;}
	double getSunlitRatio() {return SunlitRatio;}
	string getNote() {return note;}

	TStage get_stage() {return stage;}

	void setMass();
	void set_age(double x) {age=x;}
	void set_CH2O();
	void set_N(double x) {TotalNitrogen= x;} // Units are grams. Was scaled from mg in crop.cpp
	void set_HourlyNitrogenDemand (double x) {HourlyNitrogenDemand=x;}
	void set_CumulativeNitrogenDemand (double x) {CumulativeNitrogenDemand=x;}
	void set_HourlyNitrogenSoilUptake(double x) {HourlyNitrogenSoilUptake=x;}
	void set_CumulativeNitrogenSoilUptake(double x) {CumulativeNitrogenSoilUptake=x;}
	void set_C_pool_root(double x) {C_pool_root=x;}
	void set_NitrogenRatio(double x) {NitrogenRatio=x;}

	

	void update(const TWeather &, double PredawnLWP);
	void calcGasExchange(const TWeather & weather);
	void calcMaintRespiration(const TWeather&);

	double calcLeafArea();
	double calcTotalLeafMass();
	double calcActiveLeafMass();
	double calcDroppedLeafMass();
	double calcGreenLeafArea();
	double calcActualGreenLeafArea();
	double calcPotentialLeafArea();
	double calcPotentialLeafAreaIncrease(void);//calculate potential leaf area increase without carbon limitation YY
	void calcPerLeafRelativeAreaIncrease();
	double calcSenescentLeafArea();
	double calcDroppedLeafArea();
	double calcGreenLeafArea2(); // empirical fit of plant green leaf area from SPAR 02 field exp
    double calcPotentialCarbondemand(); //calculate potential carbon demand for potential leaf growth YY
	void grow();
	void C_allocation(const TWeather&);
	void calcRed_FRedRatio(const TWeather&);
	void writeNote(const TWeather &);
    

private:
	TInitInfo initInfo;
	CNodalUnit * nodalUnit;  //nodal Unit object
	CEar * ear;      //ear object
	CRoots * roots;  //root object
	CDevelopment * develop;  //phenology object
	CGas_exchange * sunlit;  //sunlit photosynthesis object
	CGas_exchange * shaded; //declare two pointers that point to a sunlit and a shaded leaf Yang 8/29/06 
	string note;
	int finalNodeNumber; //final number of nodes
	int nodeNumber; // currently initiated number of nodes 

	// Note C variables are as carbohydrate (same as dry matter)
	// masses are also carbohydrate (dry matter)
	double C_pool; // shorterm C pool, g(CH2O)
	double C_avail; //Availabel carbon from long term reserved carbon pool YY
	double C_pool_root; //storage for carbon allocated to roots but not used
	                    //in the previous time step
	double C_reserve; // longterm C pool
	double C_content;
	double C_demand;
	double C_supply;
	double C_ReserveLeaf;  //holds extra C in leaf - allows SLA to change
	double mass, seedMass,stemMass, leafMass, shootMass, rootMass, seedRootMass, earMass, activeLeafMass, droppedLeafmass; // this is redundant, but for convenience of access
	double maintRespiration;
	double sowingDay;
	double age;
	double CH2O; // carbohydrate, also dry matter, g
	double N; //Need to delete this - replaced by TotalNitrogen
	double N_pool; //SK 8/20/10: Short-term N pool for remobilization, this should come mostly from senescing leaves and can be purged daily to active leaves, not implemented at the moment
	
	double leafArea, droppedLfArea;
	double currentDroppedLfArea;
	double previousDroppedlfArea;
	double greenLeafArea,actualGreenLeafArea;
	double senescentLeafArea;
	double potentialLeafArea;
	double potentialLeafAreaIncrease; //Increase in leaf area without carbon limitation YY
	double PotentialLeafCarbonDemand; //Carbon demand for potential leaf growth without carbon limitation YY
	
	double photosynthesis_gross; // gros photosynthesis, umolCO2 m-2 s-1
	double photosynthesis_net; // gros photosynthesis, umolCO2 m-2 s-1
	double assimilate; //assimilation flux, g CO2 per plant per timestep
	double transpiration, transpirationOld; //current and previous values of transpiration - g per plant per hr
	double VPD;
	double conductance;

	double temperature;
	double shootPart; //g per plant Carbohydrate partitioined to shoot
	double rootPart;  // g per plant Carbohydrate partitioned to root
	double shootPart_old; //g per plant carbohydrate partitioned to root in the previous time step
    double rootPart_old;
	double leafPart;   //g per plant carbohydrate partitioned to leaf

	double TotalNitrogen;  //This is the total nitrogen content of the plant in g plant-1
	double HourlyNitrogenDemand;  // Nitrogen demand in g N plant-1
	double CumulativeNitrogenDemand;  // cumulativeNitrogen demand in g N plant-1
	double CumulativeNitrogenSoilUptake; // Nitrogen uptake from the soil g N plant-1
	double HourlyNitrogenSoilUptake;
	double leaf_NFraction; //records the fraction of nitrogen in leaves YY
	double leaf_N; //total nitrogen in the leaves of a plant YY (grams N per plant)
	double leaf_N_content; //leaf nitrogen content (per unit square meter) of a plant YY
	double OptimalLeafN;       //g N holds leaf N content that is optimal 
	double NitrogenRatio;     //optimal N ratio according to N Dilution ratio


	
	double emerge_gdd;//records thermal time needed for plant to emergy YY
	double SunlitLAI, ShadedLAI; // sunlit and shaded LAI values
	double SunlitRatio;  // daily ratio of sunlit leaf area to use for scaling leaf expansion due to carbon stress
    double C2_effect;
	TStage stage;



};
#endif

