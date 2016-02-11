//Class Controller
//
// Controller.h
//
//Based on CPM simulation_controller
#pragma once
#ifndef _CONTROLLER_H_
#define _CONTROLLER_H_
#include "timer.h"
#include "development.h"
#include "plant.h"
#include "weather.h"
#include "initinfo.h"
#include "radiation.h"
#ifndef FLOAT_EQ
#define EPSILON 0.001   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#define MINUTESPERDAY (24*60);
#endif


class CController
{
private:
	enum InputDataFormat {DDSOIL, ICASA, SPAR, FACE}; // input data format
	bool cropEmerged, cropHarvested;
	int	firstDayOfSim,	lastDayOfSim, SowingDay;
	// cdt 12/06/2010 - I had to increase the size of these from 132 to 133. We couldn't use
	// strcpy in vs 2008 so I changed to strcpy_s but apparently it requires enough space for the null
	// terminatino character in the destination string
	char varietyFile[133], outputFile[11], cropFile[133], logFile[132],LeafFile[133];
	char DebugFile[133];
	int iCur, // current record number
		errorFlag;


	TInitInfo			initInfo;
	CPlant*             plant;
	Timer*				time;
	TWeather*			weather;
	CDevelopment*       develop;
	InputDataFormat     weatherFormat;
public:
	CController(const char*, const char*, const char*, TInitInfo);
    ~CController();
	void setErrStatus(int ier) {errorFlag = ier;}
	int getErrStatus() {return errorFlag;}
//	char* getWeatherFile() {return weatherFile;}
	char* getInitFile() {return varietyFile;}
//	char* getOutputFile() {return outputFile;}
	char* getLogFile() {return logFile;}
	CPlant * getPlant() {return plant;}

	void addOutputMessage(char*);
    void readWeatherFrom2DSOIL(const TWeather &);
	void outputToCropFile();
	void outputToLeafFile();
	void outputToDebug();

	Timer* getTime() {return time;}

	double getFirstDayOfSim() {return firstDayOfSim;}
	double getLastDayOfSim() {return lastDayOfSim;}
	double getSowingDay()    {return SowingDay; }
	double RootWeightFrom2DSOIL;
	float MaxRootDepth, AvailableWater;

	TInitInfo getInitInfo() {return initInfo;}

	void initialize();
	int run(const TWeather &, double PredawnLWP);
	double ET_supply;   //actual supply of water to plant mol m-2 (leaf) s-1


};
#endif
