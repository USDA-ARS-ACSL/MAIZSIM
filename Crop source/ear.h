#pragma once
#include "organ.h"

class CEar :
	public COrgan
{
public:
	CEar(void);
	~CEar(void);
private:
	unsigned int totalKernels;
	double seedWeight;
};
