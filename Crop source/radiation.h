#pragma once
#define PAR2PFD 4.6  // conversion factor from W/m2 to PFD (umol m-2 s) for PAR waveband (median 550 nm of 400-700 nm) of solar radiation, see Campbell and Norman (1994) p 149
	// 4.55 is a conversion factor from W to photons for solar radiation, Goudriaan and van Laar (1994)
	// some use 4.6 i.e., Amthor 1994, McCree 1981, Challa 1995.
#define sc 1370.0   // solar constant
#include "stdafx.h"
//#using <mscorlib.dll>
//using namespace System;

class CSolar
{
public:
	CSolar(void);
	~CSolar(void);
    void SetVal(int day, double tm, double lat, double longi, double alti, double t, double I0);
//    void SetVal(DateTime day, DateTime Tm, double lat, double longi, double t, double I0);
    void SetVal(int day, double tm, double lat, double longi, double t, double I0);
    void SetVal(int day, double tm, double lat, double longi, double t);
    void SetVal(int day, double lat, double longi); // only daily global solar radiation and location data available, typical format for available weather data
	double daylength();

private:
    int jday;
	double time, latitude, longitude, altitude, tau, solarRad, PAR, NIR, halfday, stdMeridian; // To do: incorporate the computation of global irradiance including NIR
    bool isInputDaily, isInputPFD, isInputObs, isEast; // true if observed solar radiations are input
    double Sg();
	double PARfr();
    double Stot();
    double Sdr();
    double Sdf();
    double m();
    double Fdif();
protected:
	double solarnoon();
	double azimuth();
	double elev();
	double decl();
	double sin_elev();
	double cos_elev();
	double sunrise();
	double sunset();
	double PARtot();
	double PARdif();
	double PARdir();
	double pressure();
};
