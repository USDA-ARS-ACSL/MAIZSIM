#include "stdafx.h"
#include "roots.h"

CRoots::CRoots(void)
{
	Initialized=false;
}

// AD when changed to Initialized=true, the program never reaches the loop
// if (!pSC->getPlant()->get_roots()->GetInitialized()) in Crop.cpp

CRoots::~CRoots(void)
{
}
