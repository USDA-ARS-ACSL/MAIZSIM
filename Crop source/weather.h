#pragma once
#ifndef _WEATHER_H_
#define _WEATHER_H_
/* I use this structure to hold variables transfered to and from 2DSOIL */
struct TWeather
{
public:
	TWeather()
	{
		year=2004,jday=1,time=0.0, daytime = jday+time, CO2=370.0, airT=20.0, PFD=0.0, solRad=0.0,
		RH=50.0, wind=1.5, rain=0.0, dayLength=12.0, soilT = airT;
	}
	int year;
	int jday;
	double time, daytime;
	//dt added pcrl and pcrq to check balances
	double CO2, airT, PFD, solRad, RH, wind, rain, dayLength, soilT,ET_supply,LeafWP,pcrs,pcrl,
		pcrq, TotalRootWeight;
	float MaxRootDepth,ThetaAvail;
	//double CO2, airT, PFD, solRad, RH, wind, rain, dayLength, soilT;
	// DT added output flags from 2DSOIL
	int DailyOutput, HourlyOutput;
};
#endif