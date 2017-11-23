// Class CSolar
//PFD Photon Flux Density
//tau atmospheric transmissivity, HGJones 1991 considered a constant, it can be calculated from
// measured radiation/potential radiation at atmospheric surface
//SolRad Solar Radiation
// PARfr
//UseObs UseTau
//PAR Photosynthetically Active Radiation (what is PFD)
// Altitude
//time - time of day
//CosTheta Zenith angle see \eqn 11.1 Campbell and Norman (1998)
// consider using equation from GLYCIM elevation angle = 90 - zenith angle since the are measured
// differently
#pragma once
#define cPFD 4.6  // conversion factor from PAR (W/m2) to PFD (umol m-2 s) for solar radiation, see Campbell and Norman (1994) p 149
// 4.6 is a conversion factor from W to photons for visible solar radiation.
// Amthor 1994, McCree 1981, Challa 1995, Campbell and Norman 1998.
// Some use 4.55, see Goudriaan and van Laar (1994)
#define     SolarConst 1367.0   // solar constant, Iqbal (1983) Watts m-2
#define     FDIV_GUARD 1.0e-8 // floating point divison guard (tolerance)
#define     PI 3.1415

class CSolar
{
private:
	int JDay;
// Environmental Components
	double Time, Latitude, Longitude, altitude, tau;
	double DayLength, Declination, SolarNoon, SinElevation, Azimuth;
	double Sunrise, Sunset,HalfDay, Elevation, CosElevation;
	double CosTheta;
// solar components
	double PAR, SolarRadiation, PotentialSolarTotal,PotentialSolarDirect, PotentialSolarDiffuse;
	double NIRTotal,NIRDiffuse,NIRDirect,PARTotal,PARDiffuse,PARDirect;
	double PFD, NIR, PARFraction, NIRFraction;
	double PotentialPARDiffuse,PotentialPARDirect, PotentialPARTotal;
	double PotentialNIRDiffuse,PotentialNIRDirect, PotentialNIRTotal;
	double SolarDirect, SolarDiffuse, FracDiffuse,  PFDDirect, PFDDiffuse;
	double FracPARDirect, FracPARDiffuse, FracNIRDirect,FracNIRDiffuse,FracSolarDiffuse,FracSolarDirect;
	double FracPFDDirect, FracPFDDiffuse;
	double RowAzimuth; //Angle of row orientation measured from 180?

	bool useObs, useTau;
	static double SDERP[9]; //used for declination calcs
	// To do: incorporate the computation of global irradiance including NIR
	// PAR is in Wm-1, PFD is in umol m-2 s-1 for visible wavebands

	void SetDayLength();
	void SetDeclination();// Solar declination as a f(Jday)
	void SetSolarNoon();
	void SetSolarElevation() ;    // Solar height (=elevation)



	void SetAzimuth() ; // Solar azimuth measured from ?
	double press();
	double m()   ;

    void SetPotentialSolar();
	void SetPotentialPAR();
	void SetPotentialNIR();

	void SetNIRTotal();
    void SetNIRDirect();
	void SetNIRDiffuse();

    void SetPARTotal();
    void SetPARDirect();
	void SetPARDiffuse();


	void SetFracPARDirect();
	void SetFracNIRDirect();
	void SetFracDiffuse();  //Set fraction of diffuse in total light
	void SetPARFraction(double fr); // when measured PAR and Radiation are both available
	void SetPARFraction();  //When only radiation or PAR are available




public:
	CSolar(void);
	~CSolar(void);
	void SetVal(int Day, double Time, double Lat,  double Longi, double Alti, double SolRad0);

	double GetDayLength()           {return DayLength;}
	double GetDeclination()         {return Declination;}
	double GetSolarNoon()           {return SolarNoon;}
	double GetSunrise()             {return Sunrise;}
	double GetSunset()              {return Sunset;}
	double GetSinElevation()        {return SinElevation;}
	double GetSolarElevation()      {return Elevation;}
	double GetAzimuth()             {return Azimuth;}

	double GetPotentialSolarDiffuse()        {return PotentialSolarDiffuse;}
	double GetPotentialSolarTotal()          {return PotentialSolarTotal;}
	double GetPotentialSolarDirect()         {return PotentialSolarDirect;}
	double GetSolarRadiation()               {return SolarRadiation;}

	double GetPAR()                 {return PAR;}
	double GetPARFraction()         {return PARFraction;}
	double GetFracPARDirect ()      {return FracPARDirect;}
	double GetFracPARDiffuse()      {return 1.0-FracPARDirect;}

	double GetFracNIRDirect()       {return FracNIRDirect;}
	double GetNIRFraction()         {return 1.0-PARFraction;}
	double GetPFD()                 {return PAR*cPFD;}
	double GetPFDDiffuse()          {return PAR*cPFD*(1.0-FracPARDirect);}
	double GetPFDDirect()           {return PAR*cPFD*FracPARDirect;}
	double GetPFDTotal()            {return PAR*cPFD;} // sane as GetPFD
	double GetNIR()                 {return NIR;}
	double GetNIRTotal()            {return NIRTotal;}
	double GetNIRDiffuse()          {return NIRDiffuse;}
	double GetNIRDirect()           {return NIRDirect;}


};
