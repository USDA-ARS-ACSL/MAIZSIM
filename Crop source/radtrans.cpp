/*unit CanopyRadTrans;
{Basic canopy architecture parameters, 10/10/00 S.Kim
 modified to represent heterogeneous canopies
 Uniform continuouse canopy: Width2 = 0
 Hedgerow canopy : Width1 = row width, Width2 = interrow, Height1 = hedge height, Height2 = 0
 Intercropping canopy: Height1, Width1, LA1 for Crop1, and so on
 Rose bent canopy: Height1=Upright canopy, Height2 = bent portion height, 10/16/02 S.Kim}

 This is the interface into the Solar Class but calculates transmission coefficients. To use it you must
 initialize the Solar object. See Solar.cpp for more information.

 This routine adds calculations that require use of LAI and LAF (Leaf Angle Factor)
 You must pass parameters to the Solar Class to initialize that class and use the methods therin.
 The LeafAngleFactor can be numeric in which case it is used directly in one equation to calculate
 Kd as a function of zenith angle (eqn 15.4 in Campbell and Norman). In this case use 'IsLeafAngleFactorUsed' as true
 Otherwise you can use shapes as an enumeration.
*/
// dateutils ;

#include <cmath>  // need for math functions
#include <algorithm>  // need for max min functions
#include "radtrans.h"

using namespace std;

inline double cot(double a) { return 1/tan(a); }
inline double sqr(double a) {return (a * a);}


// include Math unit for Real mathmatical doubles

CRadTrans::CRadTrans()
{
	LeafAngle = Spherical;
	LeafAngleFactor = 1.0; // for spherical leaves
	IsLeafAngleFactorUsed = false;
	absorp = 0.85;   //leaf absorptivity for PAR
	clump = 1.0;
	rho_soil = 0.10; // soil reflectivity for PAR band

}
 CRadTrans::~CRadTrans()
 {
 }



 void CRadTrans::SetVal(CSolar Irradiance, double SLAI, double LeafAngleFactorIn)
 {

	 IrradianceDirect = Irradiance.GetPFDDirect();
	 IrradianceDiffuse = Irradiance.GetPFDDiffuse();
 // the transmittance values obtained from Day and Bailey (1999), chap 3, ecosystems of the world 20: greenhouse ecosystem, page 76
     LAI = SLAI;
     Elev =  Irradiance.GetSolarElevation();
     LeafAngleFactor = LeafAngleFactorIn;
     IsLeafAngleFactorUsed = true;
	Kb(GetZenith()); // Calculate  KbVal
    Kd(LAI);       // Calculate KdVal

//  GDiffuse = S;
}

double CRadTrans::GetZenith()  // need to move this to CSolar
{
	double zenith;
	// need to constrain zenith to not quite reach PI/2 when elevation is zero
	// i.e., the sun is near the horizon.
	zenith=fabs(PI/2.0-Elev);
	zenith=min(zenith,1.56);
  return zenith;
}

double CRadTrans::Reflect()
 {
	 return (1-sqrt(absorp))/(1+sqrt(absorp))*(2*GetKb()/(GetKb()+GetKd()));
 }



void CRadTrans::Kb(double theta)  // Campbell, p 253, Ratio of projected area to hemi-surface area for an ellisoid
// x is a leaf angle distribution parameter
{
	double x, tmp;

	tmp = 0.5;

	switch (LeafAngle)
	{
	case Spherical:
		tmp = 0.5/sin(Elev); // (0.5/sin_elev); When Lt accounts for total path length, division by sin(elev) isn't necessary
		break;
	case Horizontal:
		tmp = 1/sin(Elev); //1
		break;
	case Vertical:
		tmp = (2/cot(Elev))/PI;
		break;
	case Diaheliotropic:
		tmp = 1; // 1/sin_elev;
		break;
	case Empirical:
		tmp = 0.667; //(0.5)*KDiffuse/(0.8*sqrt(1-scatt)); //(0.5/sin_elev)*KDiffuse/(0.8*sqrt(1-scatt));
		break;
	default:
		tmp= 0.5/sin(Elev);
	}

	if (IsLeafAngleFactorUsed == true) {
		x = LeafAngleFactor; }
	else
	{
		switch (LeafAngle)
		{
		case Spherical:
			x = 1; //
			break;
		case Horizontal:
			x = 10; //1
			break;
		case Vertical:
			x = 0;
			break;
		case Corn:
			x = 1.37;
			break;
		default:
			x =1;
		}
	}
	//   if Sin(Elev) > sin(5) then
	//   tmp =  0.5/sin(Elev)

    tmp = sqrt(sqr(x)+sqr(tan(theta)))/(x+1.774*pow(x+1.182,-0.733));

	//   else tmp = 0.5/sin(5);

	KbVal=tmp*clump;
}


void CRadTrans::Kd(double LA)
{
const double gauss3[3] = {-0.774597,0,0.774597}; // abscissas
const double weight3[3] = {0.555556,0.888889,0.555556};

double K, FDiffuse, tmp, angle, x;
int i;
if (IsLeafAngleFactorUsed == true) {
	x = LeafAngleFactor;}
   else
   {
	     switch (LeafAngle)
		 {
		 case Spherical:
			 x = 1; //
			 break;
		 case Horizontal:
			 x = 10; //1
			 break;
		 case Vertical:
			 x = 0;
			 break;
		 case Corn:
			 x = 1.37;
			 break;
		 default:
			  x = 1;
		 }
   }	//end else

     FDiffuse = 0;
     for (i = 0; i<3; i++)  //diffused light ratio to ambient, itegrated over all incident angles from -90 to 90
       {
         angle = (PI/2)/2*(gauss3[i]) + (PI/2)/2;
         tmp = sqrt(sqr(x)+sqr(tan(angle)))/(x+1.774*pow(x+1.182,-0.733));
         FDiffuse = FDiffuse + (PI/2)/2*(2*exp(-tmp*LA)*sin(angle)*cos(angle))*weight3[i];
       }
	   if (LA <= 0.0) {
		   K = 0.0;
	   }
	   else  {K = -log(FDiffuse)/LA;}

   KdVal=K*clump;
  }


double CRadTrans::Irradiancetot() //total irradiance at the top of the canopy, passed over from either observed PAR or TSolar or TIrradiance
{
 return (IrradianceDirect + IrradianceDiffuse);  //IrradianceDirect: beam radiation at top of canopy, IrradianceDiffuse: diffuse radiation at top.
}

double CRadTrans::Qtot(double L) //total irradiance (Direct + Diffuse) at depth L, simple empirical approach
{
	double result;
	result=Irradiancetot()*exp(-sqrt(absorp)*(GetKd()+GetKb())/2*L); //;  //
	return result;
}

double CRadTrans::Qbt(double L) // total beam radiation at depth L
 {
  return (IrradianceDirect*exp(-sqrt(absorp)*GetKb()*L));
 }

double CRadTrans::Qd(double L) // net diffuse flux at depth of L within canopy
 {
  return IrradianceDiffuse*exp(-sqrt(absorp)*GetKd()*L);
 }

double CRadTrans::Qdm() // weighted average absorved diffuse flux over depth of L within canopy accounting for exponential decay
 {
	 if (LAI <= 0) {
		 return 0;}
	 else {
           return IrradianceDiffuse*(1-exp(-sqrt(absorp)*GetKd()*LAI))/(sqrt(absorp)*GetKd()*LAI); //  Integral Qd / Integral L
	 }
 }

double CRadTrans::Qb(double L) // unintercepted beam (direct beam) flux at depth of L within canopy
 {
  return IrradianceDirect*exp(-GetKb()*L);
 }

double CRadTrans::Qsl() // mean flux density on sunlit leaves
 {
  return GetKb()*IrradianceDirect + Qsh();
 }

double CRadTrans::Qsl(double L) // flux density on sunlit leaves at delpth L
 {
  return GetKb()*IrradianceDirect + Qsh(L);
 }

double CRadTrans::Qsh() // mean flux density on shaded leaves over LAI
 {
  return (Qdm() + Qsc() + Qsoilm());   // include soil reflection
 }

double CRadTrans::Qsh(double L) // diffuse flux density on shaded leaves at depth L
 {
  return (Qd(L) + Qsc(L) + Qsoilm());   // include soil reflection
 }

double CRadTrans::Qsoilm() // weighted average of Soil reflectance over canopy accounting for exponential decay
 {
	 if (LAI <= 0) {
		 return 0;}
	 else {
            return Qsoil()*rho_soil*(1-exp(-sqrt(absorp)*GetKd()*LAI))/(sqrt(absorp)*GetKd()*LAI); //  Integral Qd / Integral L
	 }
 }


double CRadTrans::Qsc() // weighted average scattered radiation within canopy
{
	double totBeam, nonscatt;
	 if (LAI == 0) {
		 return 0;}
  else
    {
      totBeam = IrradianceDirect*(1-exp(-sqrt(absorp)*GetKb()*LAI))/(sqrt(absorp)*GetKb()); // total beam including scattered absorbed by canopy
      nonscatt= IrradianceDirect*(1-exp(-GetKb()*LAI))/(GetKb());   //non scattered beam absorbed by canopy
      return (totBeam-nonscatt)/LAI; //mean scattered flux density
	}
 //  return (Qbt(LAI)-Qb(LAI));
 }

double CRadTrans::Qsc(double L) // scattered radiation at depth L in the canopy
{
    return Qbt(L)-Qb(L); // total beam - nonscattered beam at depth L
}

double CRadTrans::Qsoil() // total PFD at the soil sufrace under the canopy
 {
  return Qtot(LAI);
 }


double CRadTrans::LAIsl() // sunlit LAI assuming closed canopy; thus not accurate for row or isolated canopy
 {
	 if (Elev <= 0.01) {
		 return 0;}
	 else {
		 return (1-exp(-GetKb()*LAI))/GetKb();}
 }

 double CRadTrans::LAIsh()// shaded LAI assuming closed canopy
{
    return LAI - LAIsl();
 }


double CRadTrans::Fsl(double L) // sunlit fraction of current layer
  {
	  if (Elev <= 0.01) {
		  return 0;}
	  else {
		  return exp(-GetKb()*L);}
  }

double CRadTrans::Fsh(double L)
  {
	  return 1 - Fsl(L);
  }
