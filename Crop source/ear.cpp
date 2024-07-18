#include "stdafx.h"
#include "ear.h"
//#using <mscorlib.dll>

CEar::CEar(void)
	:COrgan(), totalKernels(0), seedWeight(0), sheathWeight(0), cobWeight(0), grainWeight(0)
{
}

CEar::~CEar(void)
{
}

// values passed into import statements are gr Carbon. 
// need to convert to gr carbohydrate
void CEar::import_sheathWeight(double dsheathCarbon)
{
	sheathWeight += dsheathCarbon / CH2O_MW * C_MW / C_conc;
}

void CEar::import_cobWeight(double dcobCarbon)
{
	cobWeight += dcobCarbon / CH2O_MW * C_MW / C_conc;
}

void CEar::import_grainWeight(double dgrainCarbon)
{
	grainWeight += dgrainCarbon / CH2O_MW * C_MW / C_conc;
}
