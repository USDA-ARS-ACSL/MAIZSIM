#include "stdafx.h"
#include "controller.h"
#include "initinfo.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#ifndef FLOAT_EQ
#define EPSILON 0.001   // floating point comparison tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#endif
#define comma ","
#define MINUTESPERDAY (24*60);

// const a = 17.27; b = 237.7; //constant in deg C
inline double E_sat(double T){return 0.6105*exp(17.27*T/(237.7+T));}
//#using <mscorlib.dll>
using namespace std;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CController::CController(const char* varFile, const char* outfile, const char* LFile, TInitInfo iniInfo)
{
	time			 = NULL;
	weather 		 = NULL;
	plant			 = NULL;
	
	strcpy_s(varietyFile, varFile);
    strcpy_s(cropFile, outfile);
	char *pch=(char*)calloc(133,sizeof(char));
	char *next=(char*)calloc(133,sizeof(char));
	char *ext=".dbg";
	// TODO fix so it does not need a buffer
	char temp[133]; 
	//get root of debug file from crop file, store in temp at first
	strcpy_s(temp,cropFile);
	// get char pointer to "." in the variable temp
	// and put root of filename into pch
	pch=strtok_s(temp,".",&next);
	//attach extension to root filename for debug file
	strcat(pch,ext);  // or strcat_s(pch,133,ext) or strncat(pch,ext,133)
	//add extension to debug file
	strcpy_s(DebugFile,pch);
	strcpy_s(LeafFile, LFile);
	initInfo = iniInfo;
	iCur = 0;
    initialize();
	errorFlag = 0;
	
}

CController::~CController()
{
	if ( time != NULL )
		delete time;
	if ( weather != NULL )
		delete [] weather ;
	if ( plant != NULL )
		delete plant;
//	if ( output != NULL )
//		delete output;
}
//**********************************************************************
void CController::initialize()
{
	#ifndef _DEBUG_FILE
	ofstream DebugOut(DebugFile, ios::out);
	DebugOut
    		 << setw(10) << "date"
			<< setw(6) << "jday"
			<< setw(7) << "time"
			<< setw(3) << "sunlit"
            << setw(4) << "shaded"
			<< setw(5) << "Pn"
			<<setw(7)  << "Pg"
			<<setw(7)  << "Lvs Init"
			<<setw(9)  << "Lvs apprd"
			<<setw(9)  << "Lvs grwg"
			<<setw(9)  << "C to Lvs"
			<<setw(9)  <<"C2_effect"
			<<setw(9)  <<"sunlitRatio"
			<<endl
        ;
	DebugOut.close();
#endif

	Timer dConvert; //object to convert dates 
	int mm, dd,yy;  // for calendar dates
	char* Buffer=(char*)calloc(256,sizeof(char)); //to hold duumy strings from variety file
	cout << "Initializing Controller object...." << endl <<endl <<endl ; 
	cout <<setiosflags(ios::left) << endl
		<< " ***********************************************************" << endl
		<< " *          MAIZSIM: A Simulation Model for Corn           *" << endl
		<< " *                     VERSION  1.1.42 2016                *" << endl
		<< " *   USDA-ARS, CROP SYSTEMS AND GLOBAL CHANGE LABORATORY   *" << endl
		<< " *   U of Washington, Environmental and Forest Sciences    *" << endl
		<< " ***********************************************************" << endl
		<< endl << endl;
// output headings of files
// This is the leaf file for output (see function "output to leaffile"
	{
		ofstream LeafOut(LeafFile, ios::out);
		LeafOut << setiosflags(ios::left)
			<< setiosflags(ios::fixed)
            << setw(10) << "date,"
			<< setw(6) << "jday,"
			<< setw(7) << "time,"
			<< setw(9) << "Lvs_Init,"
			<< setw(9) << "Lvs_Apr,"
			<< setw(9) << "Leaf_#,"
			<< setw(7) << "area,"
			<< setw(10) << "mass,"
			<< setw(10) << "Sen_Area,"
			<< setw(10) << "Pntl_Area,"
			<< setw(9) << "Elong_age,"
			<< setw(9) << "CarbRat,"
			<< setw(9) << "SLA,"
			<< setw(9) << "dropped,"
			<< setw(9) << "state,"
			<< setw(9) << "GDDSum"
			<< endl;
	}

//This is the plant file for output (see function "output to crop file"
	{
		ofstream cropOut(cropFile, ios::out); 
		cropOut << setiosflags(ios::left) 
			<< setiosflags(ios::fixed)
 			<< setw(9) << "date," 
 			<< setw(6) << "jday," 
			<< setw(8) << "time,"
			<< setw(8) << "Leaves,"
			<< setw(11)<< "MaturLvs,"
			<< setw(8) << "Dropped,"
			<< setw(9) << "LA/pl,"
			<< setw(9) << "LA_dead,"
			<< setw(8) << "LAI,"
			<< setw(8) << "RH,"
			<< setw(8) << "LeafWP,"
			<< setw(8) << "PFD,"
			<< setw(8) << "SolRad,"
			<< setw(8) << "SoilT,"
			<< setw(8) << "Tair,"
			<< setw(8) << "Tcan,"
            << setw(11) << "ETdmd,"
			<< setw(11) << "ETsply,"
			<< setw(5) << "Pn,"
			<< setw(8) << "Pg,"
			<< setw(10) << "Respir,"
			<< setw(8) << "av_gs,"
			<< setw(12) << "sunlit_LAI,"
			<< setw(12) << "shaded_LAI,"
			<< setw(12) << "sunlit_PFD,"
			<< setw(12) << "shaded_PFD,"
			<< setw(12) << "sunlit_An,"
			<< setw(12) << "shaded_An,"
			<< setw(12) << "sunlit_Ag,"
			<< setw(12) << "shaded_Ag,"
			<< setw(12) << "sunlit_gs,"
			<< setw(12) << "shaded_gs"
			<< setw(9) << "VPD,"
			<< setw(10) << "Nitr,"
			<< setw(10) << "N_Dem,"
			<< setw(10) << "NUpt,"
			<< setw(10) << "LeafN,"
			//<< setw(10) << "N_Effect,"
			<< setw(10) << "PCRL,"
		    << setw(8) << "totalDM,"
			<< setw(8) << "shootDM,"
			<< setw(8) << "earDM,"
			<< setw(10) << "TotleafDM,"
			<< setw(10) << "DrpLfDM,"
			<< setw(8) << "stemDM,"
			<< setw(8) << "rootDM,"
			<< setw(8) << "SoilRt,"
			<< setw(8) << "MxRtDep,"
            << setw(8) << "AvailW,"
			<< setw(9) << "solubleC,"
			<< setw(9) << "Note"
		    << endl;

	}

// Read variety file here
	try
	{
		ifstream cfs(varietyFile, ios::in);
		if (!cfs)
		{
			throw "Variety File not found.";
		}
		cfs.getline(initInfo.description, sizeof(initInfo.description),'\n');
		cfs.getline(initInfo.cultivar, sizeof(initInfo.cultivar),'\n') ;
		
		cfs.getline(Buffer, 255,'\n');
		cfs.getline(Buffer, 255,'\n');
        cfs >> initInfo.GDD_rating >> initInfo.genericLeafNo >> initInfo.DayLengthSensitive 
			>>initInfo.stayGreen >>initInfo.LM_min
			 >>initInfo.Rmax_LTAR >> initInfo.Rmax_LIR >> initInfo.PhyllochronsToSilk;

// end reading variety file


	    dConvert.caldat(initInfo.sowingDay,mm,dd,yy);
		if (cfs.eof()) cfs.close();
		cout << "Read variety file: " << varietyFile << endl <<endl;
		cout << "Simulation information for:" << endl;
		cout << setiosflags(ios::left)
			<< setw(20) << "Description: " << initInfo.description << endl <<endl
			<< setw(10)	<< "Cultivar: " << initInfo.cultivar << endl
		    << setw(6) << "Generic Leaf Number: " << initInfo.genericLeafNo << endl
			<< setw(6) << "Day Length Sensitive: " << initInfo.DayLengthSensitive << endl
			<< setw(6) << "Stay Green Parameter: " <<initInfo.stayGreen << endl
			<< setw(6) << "Maximum length of largest leaf: " << initInfo.LM_min <<endl
			<< setw(6) << "Rmax Leaf initiation rate: " << initInfo.Rmax_LIR << "  " << "Rmax Leaf tip appearance rate: " << initInfo.Rmax_LTAR << endl
			<< setw(6) << "Phyllochrons to Silk: " << initInfo.PhyllochronsToSilk << endl <<endl
			<< setw(6) << "Year: " << initInfo.year << endl
			<< setw(6) << "Sowing day: " <<  mm << "/" <<dd <<"/" << yy << endl 
			<< setw(6) << "TimeStep (min): " << initInfo.timeStep << endl
			<< setw(6) << "average [CO2]: " << initInfo.CO2 << endl << endl
			;

	}
	catch(const char* message)
	{
		cerr << message << "\n";
		exit(1);
	}

	firstDayOfSim = initInfo.beginDay;
	lastDayOfSim = initInfo.endDay;
	SowingDay    =initInfo.sowingDay;
	
	cropEmerged = false;
	cropHarvested = false;
    dConvert.caldat(firstDayOfSim,mm,dd,yy);
	//not sure if we need this class
    time = new Timer(dd, mm, yy,initInfo.timeStep/60.0); // Timer class gets stepsize in hours
	//dt modified this after modifying code to use julian day 
	int dim = (int)(((lastDayOfSim+1)-firstDayOfSim)*(24*60/initInfo.timeStep)); // counting total records of weather data
    weather = new TWeather[dim];

	plant	= new CPlant(initInfo);
}



void CController::readWeatherFrom2DSOIL(const TWeather & wthr)
{
	weatherFormat = DDSOIL;
	weather[iCur] = wthr;
	ET_supply = wthr.ET_supply;
	weather[iCur].daytime = wthr.jday + wthr.time;
}



int CController::run(const TWeather & wthr)
{
	readWeatherFrom2DSOIL(wthr);
    if (weather[iCur].jday >= initInfo.sowingDay && weather[iCur].jday <= lastDayOfSim)
	{
		plant->update(weather[iCur]);
		RootWeightFrom2DSOIL=wthr.TotalRootWeight;
		MaxRootDepth=        wthr.MaxRootDepth;
		AvailableWater=      wthr.ThetaAvail;
	//	plant->update(weather[iCur]);
// dt added ability to output daily based on 2dsoil output
// Always hourly for now - have to add code to average values
//		if ((weather[iCur].DailyOutput==1)&&(int(weather[iCur].time*24.0)==6))
//			{
			outputToCropFile();
//			if (plant->get_develop()->Germinated())
				outputToLeafFile();
#ifndef _DEBUG_FILE
				if (plant->get_develop()->Germinated()) outputToDebug();
#endif
//			}
//		if (weather[iCur].HourlyOutput==1)
//			{
//			outputToCropFile();
//			if (plant->get_develop()->Germinated())
//				outputToLeafFile();
//			}

		iCur++;
		time->step();
	} 
    return 0;
}
void CController::outputToCropFile()
{
	{
			int mm,id,iyyy;
			string DateForOutput;
			double av_gs = plant->get_conductance();
			if (!plant->get_develop()->Emerged())
			{
				av_gs=0;
			}
			double vpd = plant->get_VPD();
			if(vpd<0)
			{
				vpd=0;
			}
			time->caldat(weather[iCur].jday, mm, id,iyyy);
#if 0
			DateForOutput.Format("%.2d/%.2d/%4i",mm,id,iyyy);
#else
			char DateForOutputBuff[16];
			sprintf(DateForOutputBuff, "%.2d/%.2d/%4i", mm, id, iyyy);
			DateForOutput = DateForOutputBuff;
#endif
			string s = "";
			if (plant->get_develop()->Matured()) {s="Matured";}
			else if (plant->get_develop()->GrainFillBegan())   {s="grainFill";}
			else if (plant->get_develop()->Silked())           {s="Silked";}
			else if (plant->get_develop()->Flowered())         {s="Flowered";}
			else if (plant->get_develop()->Tasseled())         {s = "Tasseled"; }
			else if (plant->get_develop()->TasselInitiated())  {s="Tasselinit";}
			else if (plant->get_develop()->Emerged())          {s="Emerged";}
			else if (plant->get_develop()->Germinated())       {s="Germinated";}
    		else if (plant->get_develop()->Dead())             {s="Inactive";}
			else {s="none";}
		//	if (FLOAT_EQ(plant->get_develop()->emergence.daytime,weather[iCur].daytime)){s = "Emergence";}
		////	if (FLOAT_EQ(plant->get_develop()->tasselInitiation.daytime,weather[iCur].daytime)){s = "Tassel Initiation";}
		//	if (FLOAT_EQ(plant->get_develop()->anthesis.daytime, weather[iCur].daytime)){s = "Anthesis";}
		//	if (FLOAT_EQ(plant->get_develop()->silking.daytime, weather[iCur].daytime)){s = "Silking";}
		//	if (FLOAT_EQ(plant->get_develop()->beginGrainFill.daytime, weather[iCur].daytime)){s = "Begin grain filling";}
		//	if (FLOAT_EQ(plant->get_develop()->maturity.daytime,weather[iCur].daytime)){s = "Begin grain filling";}

		    ofstream ostr(cropFile, ios::app);
			ostr << setiosflags(ios::right) 
				<< setiosflags(ios::fixed)
				<< setw(9) << DateForOutput << comma
 				<< setw(6) << weather[iCur].jday <<comma
				<< setw(8) << setprecision(0) << weather[iCur].time*24.0 << comma 
				<< setw(8) << setprecision(2) << plant->get_develop()->get_LvsAppeared() << comma
				<< setw(8) << setprecision(2) << plant->get_nodalUnit()->get_leaf()->get_TotalMatureLeaves() << comma
				<< setw(8)  << setprecision(2) << plant->get_nodalUnit()->get_leaf()->get_TotalDroppedLeaves() << comma
				<< setw(9) << setprecision(2) << plant->calcGreenLeafArea() << comma
				<< setw(9) << setprecision(2) << plant->calcSenescentLeafArea() << comma
				<< setw(8) << setprecision(2) << plant->calcGreenLeafArea()*initInfo.plantDensity/(100*100) << comma
				<< setw(8) << setprecision(2) << weather[iCur].RH << comma
				<< setw(8) << setprecision(4) << weather[iCur].LeafWP << comma  //print out leaf water potential Yang 8/22/06
				<< setw(8) << setprecision(2) << weather[iCur].PFD << comma
				<< setw(8) << setprecision(2) << weather[iCur].solRad << comma
				<< setw(8) << setprecision(2) << weather[iCur].soilT << comma
				<< setw(8) << setprecision(2) << weather[iCur].airT << comma
				<< setw(8) << setprecision(2) << plant->get_tmpr() << comma
                << setw(10) << setprecision(3) << plant->get_ET_Old() << comma
			    << setw(10) << setprecision(3) << ET_supply << comma //in both cases transpiration is grams per plant per hour 
			
				<< setw(8) << setprecision(4) << plant->get_Pn() << comma   //g Carbo per plant per hour
				<< setw(8) << setprecision(4) << plant->get_Pg() << comma
				<< setw(8) << setprecision(4) << plant->get_MaintenanceRespiration() << comma //dt 03/2011 added to better calc mass balance g carbon per plant per hour
				<< setw(8) << setprecision(4) << av_gs << comma  //return average stomatal conductance Yang 10/31/06
				<< setw(12) << setprecision(3) << plant->get_sunlit_LAI() << comma
				<< setw(12) << setprecision(3) << plant->get_shaded_LAI() << comma
				<< setw(12) << setprecision(2) << plant->get_sunlit_PFD() << comma
				<< setw(12) << setprecision(2) << plant->get_shaded_PFD() << comma
				<< setw(12) << setprecision(4) << plant->get_sunlit_A_net() << comma
				<< setw(12) << setprecision(4) << plant->get_shaded_A_net() << comma
				<< setw(12) << setprecision(4) << plant->get_sunlit_A_gross() << comma
				<< setw(12) << setprecision(4) << plant->get_shaded_A_gross() << comma
				<< setw(12) << setprecision(4) << plant->get_sunlit_gs() << comma
				<< setw(12) << setprecision(4) << plant->get_shaded_gs() << comma
			    << setw(9) << setprecision(3) << vpd << comma
				<< setw(10) << setprecision(4) << plant->get_N() << comma
				<< setw(10) << setprecision(4) << plant->get_CumulativeNitrogenDemand() << comma
				<< setw(10) << setprecision(4) << plant->get_CumulativeNitrogenSoilUptake() << comma
				<< setw(10) << setprecision(4) << plant->get_LeafN() << comma //return mass of N in leaves YY
				//<< setw(10) << setprecision(4) << plant->
				<< setw(10)<< setprecision(4)<< plant->get_roots()->get_ActualCarboIncrement() << comma
				<< setw(8) << setprecision(3) << plant->get_mass() << comma
				<< setw(8) << setprecision(3) << plant->get_shootMass() << comma  //masses are grams per plant
				<< setw(8) << setprecision(2) << plant->get_earMass() << comma
				<< setw(8) << setprecision(2) << plant->get_leafMass() << comma
				<< setw(8) << setprecision(2) << plant->get_DroppedLeafMass() << comma
				<< setw(8) << setprecision(2) << plant->get_stemMass() << comma
				<< setw(8) << setprecision(3) << plant->get_rootMass() << comma
				<< setw(8) << setprecision(3) << RootWeightFrom2DSOIL << comma
				<< setw(8) << setprecision(1) << MaxRootDepth << comma
				<< setw(12) << setprecision(3) << AvailableWater << comma
				<< setw(8) << setprecision(2) << plant->get_C_reserve() << comma
				<< setw(20)<< setiosflags(ios::skipws) << "\"" + s + "\"" 
			//	<< setw(9) << setprecision(2)  << weather[iCur].dayLength
            //    << setw(9) << setprecision(2) << plant->get_develop()->get_GDDsum()
				<< endl; 
		ostr.close();
	}
}
void CController::outputToLeafFile()
{

	int mm,id,iyyy;
	CNodalUnit*  nU;
	CDevelopment* myDevelop=plant->get_develop();

	string DateForOutput;
	time->caldat(weather[iCur].jday, mm, id,iyyy);
#if 0
	DateForOutput.Format("%.2d/%.2d/%4i",mm,id,iyyy);
#else
	char DateForOutputBuff[16];
	sprintf(DateForOutputBuff, "%.2d/%.2d/%4i", mm, id, iyyy);
	DateForOutput = DateForOutputBuff;
#endif
	ofstream ostr(LeafFile, ios::app);
	ostr << setiosflags(ios::right) 
		<< setiosflags(ios::fixed);
	for (int i = 1; i <= myDevelop->get_LvsInitiated(); i++)
	{ 
		nU=&plant->get_nodalUnit()[i]; // note the use of "&" I needed to call a function
		                              // outside the class as opposed to calling a function in  a class
		ostr << setw(11) << DateForOutput << comma
			<< setw(7)   << weather[iCur].jday << comma
			<< setw(3)   << setprecision(0) << weather[iCur].time*24.0 << comma
			<< setw(9)   << setprecision(2) << plant->get_develop()->get_LvsInitiated() << comma
			<< setw(9)   << setprecision(2) << plant->get_develop()->get_LvsAppeared() << comma
			<< setw(9)   << setprecision (0)<< nU->get_leaf()->get_Rank() << comma
			<< setw(9)   << setprecision(3) << nU->get_leaf()->get_greenArea() << comma
//			<< scientific
			<< setw(10)   << setprecision(4) << nU->get_leaf()->get_mass() << comma
			<< fixed
			<< setw(9)   << setprecision(3) << nU->get_leaf()->get_senescentArea() << comma
            << setw(9)   << setprecision(3) << nU->get_leaf()->get_potentialArea() << comma
			<< setw(9)   << setprecision(3) << nU->get_leaf()->get_Elongation_Age() << comma
			<< setw(9)   << setprecision(3) << nU->get_leaf()->get_N_content() << comma
			<< setw(9)   << setprecision(1) << nU->get_leaf()->get_SLA() << comma
			<< setw(9)   << setprecision(3) << nU->get_leaf()->isDropped() << comma
			<< setw(9)   << setprecision(3) << nU->get_leaf()->isGrowing() << comma
			<< setw(9) << setprecision(2) << plant->get_develop()->get_GDDsum()
			<< endl; 
	}
	ostr.close();
	nU=NULL;

}

void CController::outputToDebug()
{
	// needed for saving information on carbon allocation and assimilation
	// Can be modified for other variables. 
	// only called if _DEBUG_FILE is defined.
	CNodalUnit*  nU;
	CDevelopment* myDevelop=plant->get_develop();
	int mm,id,iyyy,i;
	string DateForOutput;
	time->caldat(weather[iCur].jday, mm, id,iyyy);
#if 0
	DateForOutput.Format("%.2d/%.2d/%4i",mm,id,iyyy);
#else
	char DateForOutputBuff[16];
	sprintf(DateForOutputBuff, "%.2d/%.2d/%4i", mm, id, iyyy);
	DateForOutput = DateForOutputBuff;
#endif
	


	ofstream DebugOut(DebugFile, ios::app);
	 DebugOut << setiosflags(ios::right) 
		<< setiosflags(ios::fixed);
     DebugOut << setw(11) << DateForOutput
			<< setw(7)   << weather[iCur].jday 
			<< setw(3)   << setprecision(0) << weather[iCur].time*24.0
			<< setw(7)   << setprecision(2) << plant->get_sunlit_LAI()
	        << setw(7)   << setprecision(2) << plant->get_shaded_LAI()
			<< setw(8)   << setprecision(2) << plant->get_Pn()*1000.0   //g Carbo per plant per hour
		    << setw(8)   << setprecision(2) << plant->get_Pg()*1000.0 
			<< setw(8)   << setprecision(0) << myDevelop->get_LvsInitiated()
		    << setw(8)   << setprecision(0) << myDevelop->get_LvsAppeared()
			<< setw(8)   << setprecision(0) << plant->get_nodalUnit()->get_leaf()->get_TotalGrowingLeaves()
			<< setw(8)   << setprecision(2) << plant->get_leafPart()*1000.0 
			<< setw(6)   << setprecision(2) << plant->getC2_effect()
			<< setw(6)   << setprecision(2) <<plant->getSunlitRatio()
			<< endl;
	 myDevelop=NULL;
	DebugOut.close();
    
}
