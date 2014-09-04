#pragma once
#ifndef _ORGAN_H_
#define _ORGAN_H_
#include "weather.h"
#include "ThermalTime.h"
#include "initinfo.h"
#include "development.h"
#define CO2_MW 44.0098
#define C_MW 12.011
#define CH2O_MW 30.03


struct TElement
{
	TElement() {CHO = 0; water = 0; nitrogen=0;}
	double CHO;
	double water;
	double nitrogen;
};

class COrgan
{
public:
	COrgan();
	COrgan(const TInitInfo&);
	virtual ~COrgan();
	virtual void import_CH2O(double);// import CHO from the reserve, virtual common metabolic reserve 
	virtual void import_N(double); // import N from the reserve
//	virtual double export_CHO(); //export CHO to the reserve, virtual common metabolic reserve 
//	virtual double export_N();
//	virtual void grow(double);
	virtual void respire();
	//dt added for temporary transfer of carbon to 2DSOIL
	double get_PotentialCarboIncrement() {return PotentialCarboIncrement;}
	double get_ActualCarboIncrement() {return ActualCarboIncrement;}

	double get_age() {return age;}
	double get_physAge() {return physAge;}
	double get_mass() {return mass;}
	double get_CH2O() {return CH2O;}
	double get_N() {return N;}
	double get_Nconc() {return N/mass;}
//	TElement * get_element() {return element;}
	double get_temperature() {return temperature;}
	double get_longevity() {return longevity;}
	double get_growthDuration() {return growthDuration;}
	TInitInfo get_initInfo() {return initInfo;}
	CThermalTime * get_GDD() {return GDD;}

	//dt added for temporary 2DSOIL compatability
// dt added this function to return carbon for roots
	// this will change later as we move the root stuff into the plant model
	// from 2DSOIL
	virtual void set_PotentialCarboIncrement(double x) {PotentialCarboIncrement = x;}
	virtual void set_ActualCarboIncrement(double x) {ActualCarboIncrement = x;}
	virtual void set_age(double x) {age = x;}
	virtual void set_physAge(double x) {physAge=x;}
	virtual void set_mass(double x) {mass=x;set_CH2O(x*(C_conc/C_MW*CH2O_MW));}
	virtual void set_CH2O(double x) {CH2O=x;}
	virtual void set_N(double x) {N=x;}
//	virtual void set_element(TElement * x) {element=x;}
	virtual void set_temperature(double x) {temperature=x;}
	void set_longevity(double x) {longevity=x;}
	virtual void set_growthDuration(double x) {growthDuration=x;}
	virtual void initialize();
	virtual void update();

private:
	COrgan(const COrgan&);
	TInitInfo initInfo;
	TElement * element;
	CThermalTime * GDD;
	double age; // chronological age of an organ, days
	double physAge; // physiological age accouting for temperature effect (in reference to endGrowth and lifeSpan, days)
	double mass; // biomass, g
	double CH2O; //glucose, MW = 180.18 / 6 = 30.03 g
	double C_conc, N; //nitrogen content, mg
	double temperature; // organ temperature, C
	double longevity; // life expectancy of an organ in days at optimal temperature (fastest growing temp), days
	double growthDuration; // physiological days to reach the end of growth (both cell division and expansion) at optimal temperature, days
	double PotentialCarboIncrement; //Carbon Allocation to roots or leaves for time increment 
	double ActualCarboIncrement; //Carbon Allocation to roots or leaves for time increment  gr C for roots, gr carbo dt-1
	
};
#endif
