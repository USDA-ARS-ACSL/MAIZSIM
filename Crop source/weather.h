#pragma once
#ifndef _WEATHER_H_
#define _WEATHER_H_
/* I use this structure to hold variables transfered to and from 2DSOIL */
struct TWeather
{
public:
	TWeather()
	{
		year = 2004,
			jday = 1,
			time = 0.0,
			daytime = jday + time,
			CO2 = 370.0,
			airT = 20.0,
			PFD = 0.0,
			solRad = 0.0,
			RH = 50.0,
			wind = 1.5,
			rain = 0.0,
			dayLength = 12.0,
			soilT = airT,
			canopyT = airT,
			SoilMP_med = -100,
			PredawnLWP = -0.05,
			pcrs = 0.0,
			pcrl = 0.0,
			pcrq = 0.0,
			TotalRootWeight = 0.0,
			SoilMP_med = -.05,
			MaxRootDepth = 0.0,
			ThetaAvailRZ = 10.0,
			DailyOutput=0, 
			HourlyOutput=1;
	}
	int year;
	int jday;
	double time, daytime;
	//dt added pcrl and pcrq to check balances
	double CO2, airT, PFD, solRad, RH, wind, rain, dayLength, soilT,ET_supply,
		LeafWP, PredawnLWP, canopyT,
		pcrs,pcrl,
		pcrq, TotalRootWeight, SoilMP_med;
	float MaxRootDepth,ThetaAvailRZ;
	//double CO2, airT, PFD, solRad, RH, wind, rain, dayLength, soilT;
	// DT added output flags from 2DSOIL
	int DailyOutput, HourlyOutput;
};
#endif