#pragma once
#include "organ.h"

class CEar :
	public COrgan
{
public:
	CEar(void);
	~CEar(void);
	void import_sheathWeight(double);
	void import_cobWeight(double); 
	void import_grainWeight(double);
	double get_sheathMass() { return sheathWeight; }
	double get_cobMass() { return cobWeight; }
private:
	unsigned int totalKernels;
	double seedWeight;
	double sheathWeight;
	double cobWeight;
	double grainWeight;
};
