/*
Unit to calculate solar geometry including solar elevation, declination,
 azimuth etc using TSolar class. Data are hidden. 03/15/00 SK
- 1st Revision 10/10/00: Changed to dealing only upto the top of the canopy. Radiation transter within the canopy is now a separate module.
- added functins to calculate global radiation, atmospheric emissivity, etc as in Spitters et al. (1986), 3/18/01, SK
24Dec03, SK
- added a function to calculate day length based on Campbell and Norman (1998) p 170,
- added overloads to SetVal
2Aug04, SK
- Translated to C++ from Delphi
- revised some functions according to "An introduction to solar radiaiton" by Iqbal (1983)
- (To Do) Add algorithms for instantaneous diffuse and direct radiation predictions from daily global solar radiation for given time
- (To Do) This can be done by first applying sinusoidal model to the daily data to simulate hourly global solar radiation
- (To Do) Then the division model of diffuse and direct radiations was applied
*/
#include "stdafx.h"
#include <stdlib.h>
#include "radiation.h"
#define _USE_MATH_DEFINES
#include <math.h>

//#using <mscorlib.dll>
//using namespace System;


#define PAR2PFD 4.6  // conversion factor from W/m2 to PFD (umol m-2 s) for PAR waveband (median 550 nm of 400-700 nm) of solar radiation, see Campbell and Norman (1994) p 149
	// 4.55 is a conversion factor from W to photons for solar radiation, Goudriaan and van Laar (1994)
	// some use 4.6 i.e., Amthor 1994, McCree 1981, Challa 1995.
#define sc 1370.0   // solar constant

CSolar::CSolar(void)
{
  jday = 0;
  time = 0.5;
  latitude = 38.0;
  longitude = 75.0;
  stdMeridian = 75.0; // assume EST, 75deg W, 5 hours behind GMT, PST is 120W and 8 hours behind GMT
  tau = 0.65;
  solarRad = -99;
  PAR = -99;
  NIR = -99;
  isInputPFD = false; 
  isInputDaily = false;
  isInputObs = false;
  isEast = false; // default : longitude is west of GMT
  halfday = 6.0;
}

CSolar::~CSolar(void)
{
}

void CSolar::SetVal(int day, double lat, double longi) //Basic form: No observed PAR or solarRad available
{
    jday = day;
    time = 12.00;
    latitude = lat*M_PI/180; //convert to radians
    longitude = longi; // leave it as in degrees, used only once for solar noon calculation
    tau = 0.5;
    altitude = 50;  // Assume altitude is 50 m from the sea level
}

void CSolar::SetVal(int day, double tm, double lat, double longi, double t) //Basic form: No observed PAR or solarRad available
{
    jday = day;
    time = tm*24;
    latitude = lat*M_PI/180; //convert to radians
    longitude = longi; // leave it as in degrees, used only once for solar noon calculation
    tau = t;
    altitude = 50;  // Assume altitude is 50 m from the sea level
}

void CSolar::SetVal(int day, double tm, double lat, double longi, double t, double I0) // Observed PAR available
{
    SetVal(day,tm,lat,longi,t);
    if (I0 != -99) 
	{
		PAR = I0;
		isInputPFD = true;
	}
}
/*
void CSolar::SetVal(DateTime day, DateTime tm, double lat, double longi, double t, double I0) // TDateTime format input
{
    int dd = day.DayOfYear;
	double tt = tm.TimeOfDay.TotalHours;
    SetVal(dd, tt, lat, longi, t, I0);
}
*/

void CSolar::SetVal(int day, double tm, double lat, double longi, double alti, double t, double I0) // Altitude available
{
    SetVal(day, tm, lat, longi, t, I0);
    altitude = alti;
}
/*

void CSolar::SetRadiation(int day, double tm, double solarRad_obs, double conversion, bool isDaily) // only daily global solar radiation (MJ m-2 d-1) is input from the weather file
{
	jday = day;
	time = tm*24;
	double dailySolarRad = solarRad_obs*conversion; // conversion is a conversion factor from the unit of input solar radiation to MJ m-2 d-1

	double tm=0.5, t=0.5; // assume noon // assume partly cloudy
	SetVal(day, tm, lat, longi, t);
}

void CSolar::setRadiation(int day, double tm, double PFD_obs, double solarRad_obs) // Observed PAR and solar radiation at top of the canopy available
{
	jday = day;
	time = tm*24;
    if (PFD_obs != -99)  // if data not available input -99
	{
		PFD = PFD_obs;
		isPFDInput = true;
	}
	if (solarRad_obs != -99)
	{
		solarRad = solarRad_obs;
		isSolarRadInput = true;
	}
}


void CSolar::setLocation(double lat, double longi, double alti, bool isEastOfGMT)
{
	latitude = lat;
	longitude = longi;
	altitude = alti;
	isEast = isEastofGMT;
}

*/

double CSolar::solarnoon()
{
	const double PI = M_PI;
	double LC, EqTime, Epsil; // LC is longitude correction for Light noon, Wohlfart et al, 2000; Campbell & Norman 1998
    LC = (stdMeridian-longitude)*1/15;  // standard meridian for pacific time zone is 120 W, Eastern Time zone : 75W
	if (isEast) LC = -LC; // LC is positive if local meridian is east of standard meridian, i.e., 76E is east of 75E
    Epsil = PI*(279.575 + 0.9856*jday)/180; // convert degrees to radians
    EqTime = (-104.7*sin(Epsil)+596.2*sin(2*PI/180*Epsil)+4.3*sin(3*PI/180*Epsil)
            -12.7*sin(4*PI/180*Epsil) - 429.3*cos(Epsil)-2.0*cos(2*PI/180*Epsil)
            + 19.3*cos(3*PI/180*Epsil))/3600;  // Calculating Equation of Time
    return 12 - LC - EqTime;
}

double CSolar::daylength()
{
	//from Iqbal (1983) p 16
	double denom = cos(latitude)*cos(decl()); // this value should never become negative because -90 <= latidute <= 90 and -23.45 < decl < 23.45
	denom = __max(denom, 0.0001); // preventing division by zero for N and S poles
	double cos_omega_s = -sin(latitude)*sin(decl())/denom; // sunrise hour angle
	if (cos_omega_s > 1) halfday = 0; // in the polar region during the winter, sun does not rise
	else if (cos_omega_s < -1) halfday = 12; // white nights during the summer in the polar region
	else
	{
		//double boundLat = 85.0*M_PI/180; // latitude limit to prevent tan(90)
		try
		{
			halfday =180/M_PI*(acos(-tan(latitude)*tan(decl())))/15; // solve sin_elev for elev= 0, in degree, 90 deg = 6 hrs
			//  halfday =arccos((Sin_Elev-sin(Latitude)*sin(Decl))/(cos(Latitude)*cos(Decl))))/15; // in degree, 90 deg = 6 hrs
			/*
			if (abs(latitude) > boundLat) 
			{
				halfday = 180/M_PI*(acos(-tan(boundLat)*tan(decl())))/15;
			}
			else
			{
				halfday =180/M_PI*(acos(-tan(latitude)*tan(decl())))/15; // solve sin_elev for elev= 0, in degree, 90 deg = 6 hrs
			}
			*/
		}
		catch(int code)
		{
//			Console::WriteLine(pe->ToString());
		}
	}
    if (halfday >=0) return halfday*2;
    else return 0;
}

double CSolar::sunrise()
{
    return solarnoon() - halfday;
}

double CSolar::sunset()
{
    return solarnoon() + halfday;
}


double CSolar::decl()
{
	double gamma = 2*M_PI*(jday-1)/365; // day angle
	double delta = 0.0;
    //   Decl = -23.45*PI/180*Cos(2*PI*(JDay+10)/365); // Goudriaan 1977
    //   Decl = 23.5*PI/180*Cos(2*PI*(JDay-172)/365);       // Resenberg, blad, verma 1982
	delta = 0.006918-0.399912*cos(gamma)+0.070257*sin(gamma)-0.006758*cos(2*gamma)+0.000907*sin(2*gamma)
		-0.002697*cos(3*gamma)+0.00148*sin(3*gamma); // Spencer equation, Iqbal (1983) Pg 7 Eqn 1.3.1. Most accurate among all
    // delta = 23.45*M_PI/180 * sin(2*M_PI / 365 * (284 + jday)); // Iqbal (1983) Pg 10 Eqn 1.3.3, and sundesign.com
	return delta;
    //   Decl = arcsin(0.39785*sin(DegToRad(278.97+0.9856*Jday+1.9165*sin(DegToRad(356.6+0.9856*Jday)))));  // Campbell and Norman, p168
}  

double CSolar::elev()
{
    //   if Sin_elev < 0 then Elev = 0
    //   else
    return asin(sin_elev());
}


double CSolar::sin_elev() // When time gets the same as solarnoon, this function fails. 3/11/01 ??
{
	double value = 0.0;
	try
	{
		value = sin(latitude)*sin(decl())+cos(latitude)*cos(decl())*cos((time-solarnoon())*2*M_PI/24);
	}
	catch(...)
	{
	}
	return value;
}



double CSolar::cos_elev()
{
	double value = 0.0;
	try
	{
		//value = sqrt(1.0-sin_elev*sin_elev);
	}
	catch(...)
	{
	}
	return value;
}
 

double CSolar::azimuth()
{
    /*	
    The solar azimuth angle is the angular distance between due South and the
    projection of the line of sight to the sun on the ground.
    View point from south, morning: +, afternoon: -
	See An introduction to solar radiation by Iqbal (1983) p 15-16
	Also see http://www.susdesign.com/sunangle/index.html
    */
	double value = 0.0;
	try
	{
		value = acos((sin_elev()*sin(latitude)-sin(decl()))/(cos(latitude)*cos_elev()));
	}
	catch(...)
	{
	}
    if (time < solarnoon())
	{
        return  value;
	}
    else if (time > solarnoon())
	{
		return  -value;
	}
    else return 0;
}



double CSolar::Sg() // Gourdriaan and van Laar's global solar radiation
{
	return sc*tau*sin_elev()*(1+0.033*cos(2*M_PI*(jday-10)/365));
     // tau: atmospheric transmissivity, Gourdriaan and van Laar (1994) p 30
}


double CSolar::Stot() // Cambell and Norman's global solar radiation, this approach is used here
{
    return Sdr() + Sdf();
}

double CSolar::Sdr()
{
    if (sin_elev() >= 0) //  result = Sg-Sdf
	{
        return sc*pow(tau,m())*sin_elev()*(1+0.033*cos(2*M_PI*(jday-10)/365));
	}
    else return 0;
}


double CSolar::pressure() // atmospheric pressure in kPa
{
	double value = 100.0;
	try
	{
		value = 101.3*exp(-altitude/8200); // campbell and Norman (1998), p 41
	}
//	catch(ArithmeticException * pe)
	catch (int code)
	{
//		Console::WriteLine(pe->ToString());
		value = 100.0;
	}
	return value;
}

double CSolar::m() // the optical air mass number
{
	double xx = __max(0.0001, sin_elev());
	return pressure()/(101.3*xx); // campbell and Norman (1998), p 173
}

double CSolar::Sdf()
{
    if (sin_elev() >= 0)
	{
		// result = sg*Fdif
        return 0.3*(1- pow(tau,m()))*sc*sin_elev()*(1+0.033*cos(2*M_PI*(jday-10)/365));   // campbell and Norman (1998), p 173
        // result = sc*sin_elev*(1-power(tau,m))/(1-1.4*ln(tau))/2 //Takakura (1993), p 5.11
	}
    else return 0;
}

double CSolar::PARfr() // PAR fraction
{

    // if elev <= 0.0 then result = 0 else
    if (tau >= 0.7) return 0.45;         // clear sky (tau >= 0.7): 45% is PAR      Goudriaan and van Laar (1994)
    else if (tau <= 0.3) return  0.55;    // cloudy sky (<= 0.3): 55% is PAR
    else return 0.625 - tau*0.25;
}

double CSolar::PARtot() // total PAR (umol m-2 s-1) on horizontal surface
{
	if (isInputPFD) return PAR;
    else return Stot()*PARfr()*PAR2PFD;
}


double CSolar::PARdif()
{
    if (isInputPFD) return PARtot()*Fdif();
    else return Sdf()*PARfr()*PAR2PFD;
}

double CSolar::PARdir()
{
	double tmp = 0.0;
    tmp = PARtot()*(1-Fdif());
    if (isInputPFD)
	{
		return tmp;
	}
    // PAR is measured using a flat sensor. Canopy rad trans model requires the flux density on the surface normal to the beam
    else return Sdr()*PARfr()*PAR2PFD;
}

double CSolar::Fdif() //Fraction of diffused light
{
	double df, dr;
	/*
    if (tau >= 0.7) return 0.2;  // clear sky : 20% diffuse
    else if tau <= 0.3 return 1;  // cloudy sky: 100% diffuse
    else return (-2*tau+1.6);         // inbetween
	*/
    // df = (1-power(tau,m))/(1-1.4*ln(tau))/2;
    df = (1-pow(tau,m()))*0.3;
    dr = pow(tau,m());
    return df/(dr+df);
}


