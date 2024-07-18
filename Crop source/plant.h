#pragma once
#ifndef _PLANT_H_
#define _PLANT_H_
#include "organ.h"
#include "nodalunit.h"
#include "development.h"
#include "roots.h"
#include "ear.h"
#include "gas_exchange.h"
//#include "gas_ex_species_param.h"
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
	CPlant(const TInitInfo&, TGasExSpeciesParam&);
	~CPlant();

	CNodalUnit* get_nodalUnit() { return nodalUnit; }
	CEar * get_ear() { return ear; }
	CRoots * get_roots() { return roots; }
	CDevelopment * get_develop() { return develop; }
	CGasExchange * get_sunlit() { return this->sunlit; }
	CGasExchange * get_shaded() { return this->shaded; }  //get access to the pointers that point to the sunlit/shade leaves Yang 8/29/06

	int get_nodeNumber() { return nodeNumber; }
	int get_finalNodeNumber() { return finalNodeNumber; }

	double get_mass() { return mass; }
	double get_age() { return age; }
	double get_CH2O() { return CH2O; }
	double get_TotalN() { return TotalNitrogen; } //Total N in plant mg plant-1
	double get_Pg() { return photosynthesis_gross; }
	double get_Pn() { return photosynthesis_net; }
	double get_assimilate() { return assimilate; }
	double get_ET() { return transpiration; }
	double get_ET_Old() { return transpirationOld; }
	double get_tmpr() { return temperature; }
	double get_C_pool() { return C_pool; }
	double get_C_pool_root() { return C_pool_root; }
	double get_C_reserve() { return C_reserve; }
	double get_MaintenanceRespiration() { return maintRespiration; }
	double get_stemMass() { return stemMass; }
	double get_leafMass() { return leafMass; }
	double get_earMass() { return earMass; }
	double get_cobMass() { return cobMass; }
	double get_sheathMass() { return sheathMass; }
	double get_grainMass() { return grainMass; }
	double get_shootMass() { return shootMass; }
	double get_rootMass() { return rootMass; }
	double get_shootPart() { return shootPart; }
	double get_leafPart() { return leafPart; }
	double get_rootPart() { return rootPart; }
	double get_DroppedLeafMass() { return droppedLeafmass; }
	double get_activeLeafMass() { return activeLeafMass; }
	double get_conductance() { return conductance; }
	double get_VPD() { return VPD; }
	double get_LeafArea() { return leafArea; }
	double get_LeafN() { return leaf_N; }
	double get_LeafNFraction() { return leaf_NFraction; }
	double get_Leaf_N_Content() { return leaf_N_content; }
	double get_HourlyNitrogenDemand() { return HourlyNitrogenDemand; }
	double get_CumulativeNitrogenDemand() { return CumulativeNitrogenDemand; }
	double get_HourlyNitrogenSoilUptake() { return HourlyNitrogenSoilUptake; }
	double get_CumulativeNitrogenSoilUptake() { return CumulativeNitrogenSoilUptake; }
	double get_droppedLfArea() { return currentDroppedLfArea; }
	double get_ptnLfIncrease() { return potentialLeafAreaIncrease; }
	double get_sunlit_LAI() { return sunlit_LAI; }
	double get_shaded_LAI() { return shaded_LAI; }
	double get_sunlit_PFD() { return sunlit_PFD; }
	double get_shaded_PFD() { return shaded_PFD; }
	double get_sunlit_A_net() { return sunlit_A_net; }
	double get_shaded_A_net() { return shaded_A_net; }
	double get_sunlit_A_gross() { return sunlit_A_gross; }
	double get_shaded_A_gross() { return shaded_A_gross; }
	double get_sunlit_gs() { return sunlit_gs; }
	double get_shaded_gs() { return shaded_gs; }
	double getC2_effect() {return C2_effect; }
	double getSunlitRatio() {return SunlitRatio; }
	string getNote() { return note; }

	TStage get_stage() {return stage;}

	void setMass();
	void set_age(double x) {age=x;}
	void set_CH2O();
	void set_Q10MR(double x) { Q10MR = x; }
	void set_TotalN(double x) {TotalNitrogen= x;} // Units are grams per plant. Was scaled from mg in crop.cpp
	void set_HourlyNitrogenDemand (double x) {HourlyNitrogenDemand=x;}
	void set_CumulativeNitrogenDemand (double x) {CumulativeNitrogenDemand=x;}
	void set_HourlyNitrogenSoilUptake(double x) {HourlyNitrogenSoilUptake=x;}
	void set_CumulativeNitrogenSoilUptake(double x) {CumulativeNitrogenSoilUptake=x;}
	void set_C_pool_root(double x) {C_pool_root=x;}
	void set_NitrogenRatio(double x) {NitrogenRatio=x;}

	

	void update(const TWeather&);
	void calcGasExchange(const TWeather & weather, const TGasExSpeciesParam& photoparam);
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
	double calcPotentialCarbondemand(); //calculate potential carbon demand for potential leaf growth YY
	void grow();
	void C_allocation(const TWeather&);
	void calcRed_FRedRatio(const TWeather&);
	void writeNote(const TWeather &);
	void calcLeafN_Content();
    
	

private:
	TInitInfo initInfo;
	TGasExSpeciesParam gasExparam;
	CNodalUnit * nodalUnit;  //nodal Unit object
	CEar * ear;      //ear object
	CRoots * roots;  //root object
	CDevelopment * develop;  //phenology object
	CGasExchange * sunlit;  //sunlit photosynthesis object
	CGasExchange * shaded; //declare two pointers that point to a sunlit and a shaded leaf Yang 8/29/06 
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
	double N_effectOnKernal;  //N effect on kernel development - to slow it down
	double C_ReserveLeaf;  //holds extra C in leaf - allows SLA to change
	double mass, seedMass,stemMass, leafMass, shootMass, rootMass, seedRootMass, earMass, activeLeafMass, droppedLeafmass, cobMass, sheathMass, grainMass; // this is redundant, but for convenience of access
	double maintRespiration;
	double sowingDay;
	double age;
	double CH2O; // carbohydrate, also dry matter, g
	double N_pool; //SK 8/20/10: Short-term N pool for remobilization, this should come mostly from senescing leaves and can be purged daily to active leaves, not implemented at the moment
	double Q10MR; //Q10 for maintenance respiration
	
	double leafArea, droppedLfArea; // leaf area, m2 per plant
	double currentDroppedLfArea; //leaf area dropped in the current time step
	double previousDroppedlfArea;
	double greenLeafArea,actualGreenLeafArea; //green leaf area, m2 per plant
	double senescentLeafArea; // senescent leaf area, m2 per plant
	double potentialLeafArea; // potential leaf area, m2 per plant with no stresses
	double potentialLeafAreaIncrease; //Increase in leaf area without carbon limitation YY m-2 per plant
	double PotentialLeafCarbonDemand; //Carbon demand for potential leaf growth without carbon limitation m-2 YY
	
	double photosynthesis_gross; // gross photosynthesis, umolCO2 m-2 s-1
	double photosynthesis_net; // gross photosynthesis, umolCO2 m-2 s-1
	double assimilate; //assimilation flux, g CO2 per plant per timestep
	double transpiration, transpirationOld; //current and previous values of transpiration - g per plant per hr
	double VPD;
	double conductance;
	double LAF; // leaf angle factor for corn leaves, Campbell and Norman (1998)


	double temperature;
	double shootPart; //g per plant Carbohydrate partitioined to shoot
	double rootPart;  // g per plant Carbohydrate partitioned to root
	double shootPart_old; //g per plant carbohydrate partitioned to root in the previous time step
    double rootPart_old;
	double leafPart;   //g per plant carbohydrate partitioned to leaf

	double TotalNitrogen;  //This is the total nitrogen content of the plant in g plant-1
	
	double CumulativeNitrogenDemand;  // cumulativeNitrogen demand in g N plant-1
	double CumulativeNitrogenSoilUptake; // cumulative Nitrogen uptake from the soil g N plant-1
	double HourlyNitrogenSoilUptake;   // Nitrogen uptake from the soil g N plant-1 hour-1
	double HourlyNitrogenDemand;  // Nitrogen demand in g N plant-1 hour-1
	double leaf_NFraction; //records the fraction of nitrogen in leaves YY
	double leaf_N; //total nitrogen in the leaves of a plant YY (grams N per plant)
	double leaf_N_content; //leaf nitrogen content (per unit square meter) of a plant YY
	double OptimalLeafN;       //g N holds leaf N content that is optimal not used now
	double NitrogenRatio;     //optimal N ratio according to N Dilution ratio g n per g biomass

	
	double emerge_gdd;//records thermal time needed for plant to emergy YY
	double sunlit_LAI, shaded_LAI; // sunlit and shaded LAI values
	double sunlit_PFD, shaded_PFD;
	double 	sunlit_A_net, shaded_A_net,
		sunlit_A_gross, shaded_A_gross,
		sunlit_gs, shaded_gs;

	double SunlitRatio;  // daily ratio of sunlit leaf area to use for scaling leaf expansion due to carbon stress
    double C2_effect;
	TStage stage;



};
#endif

