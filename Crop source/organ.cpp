#include "stdafx.h"
#include "organ.h"


COrgan::COrgan()
{
	age = physAge = mass = 0;
	CH2O=N=0;
	C_conc = 0.40;
	temperature = 25;
	growthDuration=10;
	longevity=50;
	GDD = NULL;
	GDD = new CThermalTime();
	//dt added for temporary transfer of carbon to 2DSOIL
    PotentialCarboIncrement = 0;
	ActualCarboIncrement=0;
	
}

COrgan::COrgan(const TInitInfo& info)
{
	initInfo = info;
	temperature=25.0;
	CH2O=N=0;
	C_conc = 0.40;
	age = physAge = mass = 0;
	growthDuration=10;
	longevity=50;
	GDD = NULL;
	GDD = new CThermalTime();
	PotentialCarboIncrement = 0; //Amount of potential carbo needed for new growth
	ActualCarboIncrement=0;      // Amount of actual carbo available for growth
}

COrgan::~COrgan()
{
	if (GDD != NULL) delete GDD;
}

void COrgan::initialize()
{
	if (GDD == NULL) GDD = new CThermalTime();
	GDD->initialize(initInfo.timeStep);
}
void COrgan::update()
{
	GDD->add(temperature);
	age = GDD->get_actualAge();
	physAge=GDD->get_sum();
	mass = CH2O/CH2O_MW*C_MW/C_conc; // C content as carbohydrate
}

void COrgan::import_CH2O(double dCH2O)
{
	CH2O += dCH2O;
	mass = CH2O/CH2O_MW*C_MW/C_conc; // C content as carbohydrate
}

void COrgan::import_N(double dN)
{
	N += dN;
}

void COrgan::respire()
// this needs to be worked on
// currently not used at all
{
	double Rm = 0.02; //maintenance respiration
	double Ka = 0.1; //growth respiration
	CH2O -= Ka*CH2O + Rm*CH2O;
}




