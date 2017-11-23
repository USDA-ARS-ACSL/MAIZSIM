#include "solar.h"
class CRadTrans
{
private:

	enum CLeafAngle {Spherical, Horizontal, Vertical, Diaheliotropic, Empirical, Ellipsoidal, Corn};
	//type TCover = (Glass, Acrylic, polyethyl, doublepoly, whitewashed, NoCover);
	double      absorp;   //leaf absorptivity for PAR
	double      clump ;
	double      rho_soil ; // soil reflectivity for PAR band
	double IrradianceDirect, IrradianceDiffuse, LAI, Elev, LeafAngleFactor, KbVal, KdVal;
	CLeafAngle LeafAngle; // Solar elevation of a beam and cumulative LAI at the layer, diffused fraction (fdf)
	bool IsLeafAngleFactorUsed;

	double Qbt(double L);
	double Qb(double L);
	double Qd(double L);
	double Qsoil();
	void   Kb(double theta) ;
	void   Kd(double LA);
public:
	CRadTrans(void); //sdf: diffused fraction of solar radiation
	~CRadTrans(void);

	void SetVal(CSolar Irradiance, double LAI, double leafAngleFactor);
	double Qsc();
	double Qsc(double L);
	double Qtot(double L); //total irradiance (dir + dif) at depth L, simple empirical approach
	double Irradiancetot(); // total PAR at top of the canopy
	double Qsl();
	double Qsh();
	double Qsl(double L);
	double Qsh(double L);
	double Qdm() ;
	double Qsoilm();
	double GetZenith();
	double Reflect();
	double LAIsl();
	double LAIsh();
	double Fsl(double L);
	double Fsh(double L);
	double GetKb(){return KbVal;} ;    // extiction coefficient assuming spherical leaf dist

	double GetKd(){return KdVal;}


};
