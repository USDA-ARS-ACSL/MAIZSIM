#pragma once
#ifndef _INITINFO_H_
#define _INITINFO_H_
#define MINUTESPERDAY (24*60);
#ifndef FLOAT_EQ
#define EPSILON 0.001   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#define MINUTESPERDAY (24*60);
#endif

struct TInitInfo
{
public:
	TInitInfo()
	{
		char description[255] = "\0";
		char cultivar[255]="\0";
		GDD_rating = 1331;
		genericLeafNo=15;
		latitude = 38.0; longitude = 0.0; altitude = 50.0;
		sowingDay = 150;
		beginDay = 1; endDay = 365;
		year = 2004;
		timeStep=5.0;
		plantDensity = 8.0;
		CO2 = 370.0;
		Rmax_LIR=.0978;
		Rmax_LTAR = 0.53;
		DayLengthSensitive=true;
		PhyllochronsToSilk=8;
	}
	char description[255];
	char cultivar[255];
	int GDD_rating; // GDD or GTI rating of the cv, see Stewart 1999 for conversion between MRMR and other ratings
	int genericLeafNo; // leaf number at the end of juvenile phase independent of environmental ques of leaf initiation
	double plantDensity;
	double latitude, longitude, altitude;
	int sowingDay, beginDay, endDay;
	double CO2;
	int year;
	double timeStep;
	bool DayLengthSensitive; //1 if daylength sensitive
	double Rmax_LIR, Rmax_LTAR; //  Maximum Leaf tip initiation and appearance rates
	double PhyllochronsToSilk; //number of phyllochrons from tassel initiation for 50% silking. 
	// routines for date manipulation from fortran
	
};
#endif