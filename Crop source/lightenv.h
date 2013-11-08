#pragma once
#ifndef _DELPHI_LIGHTENV_H_
#define _DELPHI_LIGHTENV_H_
#ifndef DLLIMPORT
#define DLLIMPORT __declspec(dllimport)  __stdcall
#endif
//import delphi dll for canopy light environment calculations
//TODO: to be translated to C++
//#using <mscorlib.dll>

extern "C"
{
	void DLLIMPORT  setSolarEnv(int Day, double Lati, double Longi);
	double DLLIMPORT  getDaylength();
	double DLLIMPORT  getSunrise();
	double DLLIMPORT  getSunset();
	void DLLIMPORT  radTrans(int, double, double, double, double, double, double);
	void DLLIMPORT  radTrans2(int, double, double, double, double, double, double, double);
	void DLLIMPORT  setRad(int, double, double, double, double);
	void DLLIMPORT  setRad2(int jday, double tm, double Lati, double Longi, double I0, double PFD0);
	void DLLIMPORT  setPFD(int Jday, double Time, double Lati, double Longi, double Alti, double tau, double I0, double LAI, double LAF); // LAF: leaf angle factor = 1.37 for corn leaves
	double DLLIMPORT  sunlitPFD();
	double DLLIMPORT  shadedPFD();
	double DLLIMPORT  sunlitLAI();
	double DLLIMPORT  shadedLAI();
	double DLLIMPORT  K_b();
	double DLLIMPORT  getPFDtot();
	double DLLIMPORT  getPFDdir();
	double DLLIMPORT  getPFDdif();
	double DLLIMPORT  getNIRtot();
	double DLLIMPORT  getNIRdir();
	double DLLIMPORT  getNIRdif();
	double DLLIMPORT  getSolRad();
}
#endif