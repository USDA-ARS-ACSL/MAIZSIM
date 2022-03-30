#ifndef __NODALUNIT_H__
#define __NODALUNIT_H__ 

#include "organ.h"
#include "leaf.h"
#include "stem.h"
#include "development.h"

class CNodalUnit
{
public:
	CNodalUnit();
	~CNodalUnit();
	int get_rank() {return rank;}
	bool isInitiated() {return initiated;}
	bool isAppeared() {return appeared;}
	bool isGrowing() {return growing;}
	bool isProlific() {return prolific;}
	bool isAging() {return aging;}
	bool isTerminated() {return terminated;}
    CLeaf * get_leaf() {return leaf;}
	CStem * get_stem() {return stem;}
//	CSheath * get_sheath() {return sheath;}
//	CInternode * get_internode() {return internode;}
	void set_leaf(CLeaf * x) {leaf=x;}
	void set_stem(CStem * x) {stem=x;}
//	void set_sheath(CSheath * x) {sheath=x;}
//	void set_internode(CInternode * x) {internode=x;}
	void update(CDevelopment *, double PredawnLWP);
	void initialize(int, CDevelopment * dv);
//	double get_leafLength(int rank);
private:
	int rank; 
	bool initiated, appeared, growing, prolific, aging, terminated;
	CLeaf * leaf;
	CStem * stem;
	double mass;
//	CSheath * sheath;
//	CInternode * internode;
};
#endif

