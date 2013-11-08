#include "stdafx.h"
#include "stem.h"

CStem::CStem()
:COrgan(), length(0.0), diameter(0.0) {}

CStem::CStem(int n): COrgan()
{
	rank = n;
	length = diameter = 0.0;
}
CStem::~CStem() {}

void CStem::update(CDevelopment * dv)
{
	COrgan::set_temperature(dv->get_Tcur());
	COrgan::update();
}