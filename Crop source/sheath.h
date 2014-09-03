#pragma once
#ifndef _SHEATH_H_
#define _SHEATH_H_
#include "organ.h"


class CSheath :
	public COrgan
{
	friend class CLeaf;
public:
//	CSheath(void);
	CSheath(int);
	virtual ~CSheath(void);
	void initialize();
	void update(CDevelopment *);
	void elongate(CDevelopment *);
	double get_length(){return length;}
	double get_ptnLength(){return ptnLength;}
private:
	CSheath(const CSheath&);
	bool initiated;
	int rank;
	double length;
	double area;
	double ptnLength, PotentialArea;
	double elongRate, elongAge;
	double phase1Delay, growthDuration;
	double plastochrons;

};
#endif
