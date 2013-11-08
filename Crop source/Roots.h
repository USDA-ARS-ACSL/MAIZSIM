#pragma once
#include "organ.h"

class CRoots :
	public COrgan
{
public:
	CRoots(void);
	virtual ~CRoots(void);
	bool GetInitialized() {return Initialized;}
	void SetInitialized() {Initialized=true;}
private:
	bool Initialized;
};
