#include "stdafx.h"
#include "nodalunit.h"
#include "weather.h"

CNodalUnit::CNodalUnit()
{
	rank = 0; // coleoptile
	leaf = NULL;
	stem = NULL;
//	sheath = NULL;
//	internode = NULL;
	initiated = appeared = growing = aging =  prolific =terminated = false;
}
void CNodalUnit::initialize(int n, CDevelopment * dv)
{
	rank = n;
	leaf = new CLeaf(n, dv);
	stem = new CStem(n);
//	sheath = new CSheath();
//	internode = new CInternode();
	leaf->initialize(dv);
	stem->initialize();
	mass = leaf->get_mass() + stem->get_mass(); // has no mass here
	initiated = true;
}

CNodalUnit::~CNodalUnit()
{
	if (leaf !=NULL) delete leaf;
	if (stem !=NULL) delete stem;
//	delete sheath;
//	delete internode;
}


void CNodalUnit::update(CDevelopment * dv, double PredawnLWP)
{
	leaf->update(dv, PredawnLWP);
	stem->update(dv);
//	stem->grow(weather); // mass is not updated yet for leaf area increase should be able to add actual carbon increment here
	mass = leaf->get_mass() + stem->get_mass(); // has no mass here
}


