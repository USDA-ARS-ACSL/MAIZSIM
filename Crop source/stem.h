#pragma once
#ifndef _STEM_H_
#define _STEM_H_
#include "organ.h"

class CStem: public COrgan
{
public:
	CStem();
	CStem(int);
	~CStem();

	double get_length() {return length;}
	double get_diameter() {return diameter;}

	void set_length(double x) {length=x;}
	void set_diameter(double x) {diameter=x;}
	void update(CDevelopment * dv);

private:
	int rank;
	double length;
	double diameter;
};
#endif