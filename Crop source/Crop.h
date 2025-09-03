// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the PLANT_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// PLANT_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.

// See creating a pocket PC dll using C++ article for information
// on how to do this
#ifdef _WIN32
#ifdef MYPLANT_EXPORTS
#define PLANT_API __declspec(dllexport)
#else
#define PLANT_API __declspec(dllexport)
#endif
#else
#define PLANT_API
#endif


#include "initinfo.h"

// Other global definitions go here

    double p_VMAX,
		  PopSlab,
		  Trel,            // relative time since emergence
		  Period,          // one hour of time
		  Emergence;       // Day of year plant emerged
   //double PredawnLWP=-0.05; moved to tweather
   
	TInitInfo			initInfo; //Initial input data stored here
  double old_shoot_weightPerM2 = 0; //Declare and initilize a variable to save 
                               //the above-ground biomass in each time step YY  

	bool  bEmergence=false; //boolean for emergence occurring
	int   ModNum=0, //Module Number
		  errPlant=1    // error number to return to 2DSOIL
          ;
	const double CO2=370, Press=98;
	const double cm2perm2=10000;
	double  TotalTillerLeafArea, TotalLeafArea, TillerLeafArea, 
		    MainStemLeafArea;
	double netPhotosyn,respiration, GPhotosyn;
    int nNodes, nLeaves, ntLeaves, nTiller;
	const double maxAge=47;
	double WaterUptake=0; // hourly water uptake from 2dsoil
	double NitrogenUptake=0; // nitrogen uptake value from 2dsoil accumulated between time steps mg/plant
	double NitrogenUptakeOld=0;
	double CumulativeNitrogenDemand=0.0; //grams plant-1
	double CumulativeActualNFromSoil=0.0; //grams plant-1
    double U_N, U_M, U_P, U_D; //U_N  maximum observed N uptake rate (g N m-2 ground d-1) (Lindquist et al, 2007) YY
	                           //U_M maximum uptake rate as limited by maximum N fraction per unit (Equation 2 in Lindquist et al., 2007)
                               // U_P potential rate of N accumulation (g N m-2 ground d-1) (Lindquist et al. 2007)
                               //U_D U uptake rate (g N m-2 d-1) as limited by the difference between potential and actual amount of N 
				               //in existing biomass, equation 3 in Lindquist et al. 2007)
	

	double d = 0.075; //d: shape coefficient in the logistic function to simulate cumulative N uptake (Equation 9 in Lindquist et al. 2007)
	double q_n = 0.032; //q_n the maximum ratio of daily N uptake to measured daily growth rate (g N g-1) (Lindquist et al., 2007)
	double N_min = 4.10; //maximum nitrogen concentration for C4 species: 4.1% (Lindquist et al. 2007)
	double N_shape = 0.5; //shape coefficient in calculation of nitrogen concentration in relation to up-ground biomass (equation 4 Lindquist et al, 2007) YY
//Common Structures defined here
	double massIncrease=0; //increase in shoot biomass in each time step YY
	double shoot_weightPerM2=0;
	double SLNmin = 0.5; //SLNmin: Base specific leaf nitrogen content 
	double CurrentNUptakeError=0;
	double CumulativeNUptakeError=0;
	double ET_diff=0; //for debugging

    float LAMDAS, LAMDAC;
	const int NumNPD=4000, NumElD=3500, NumBPD=600, NSeepD = 2,
              NumSPD= 30, NumSD =10, NDrainD=2, NumDR=30, 
			  NumGD = 3, NumPlD=100, 
              NMatD=15, NumModD=20, MBandD=15,
	          NumSurfDatD=3+NumGD+NumSD;
	int compare(const void *arg1, const void *arg2); //Function for comparing two numbers in a sort routine
#pragma pack(2)
 struct ShootCommon{
    double PCRL,PCRQ,PCRS,HourlyCarboUsed,ET_demand,LCAI,Cover,Convr;
	float MaxRootDepth,Shade,Height,LAI,AWUPS,nitroDemand;
	float xBStem,yBStem,SGT,PSIM,
      LAREAT,PopRow,RowSp,RowAng,PopArea,CEC,
      EORSCS,AWUPSS,SolRad,Total_Eor,
      Total_Pcrs,SIncrSink,Psild,
      OsmFac, EOMult,LeafWP, NDemandError, CumulativeNDemandError, 
	  TotalRootWeight, InitialRootCarbo,
	  ConstI[2],constK[2], Cmin0[2];
	int isGerminated, isEmerged;
  };
// DT Made IR float (from int) 
 
 //Weather
 struct WeathCommon{
	 int   MSW1,MSW2,MSW3,MSW4,MSW5,MSW6,MSW7;
	 float BSOLAR,ETCORR,BTEMP,ATEMP,ERAIN,BWIND,BIR,WINDA, IRAV;
	 int   JDAY, NCD,JDLAST;
     float CLDFAC,DEL[24],RINT[24],RNS,RNC,RAIN,IR;
	 float WIND,CO2,TDUSK,TDUSKY,CPREC[NumSD],TAIR[24],VPD[24],ROUGH,
           RADINT[24],WATTSM[24],DIFINT[24],ROWINC[24];
	 float CLOUD, SHADOW[24],DIFWAT[24],DIRINT[24];
	 float WATACT, WATRAT, WATPOT, RNLU;
	 int   NumF[40],NumFP;
	 float hFur[40],QF;
	 int   IFUR;
	 float GAIR[NumGD],PG,LATUDE,Longitude, Altitude,
		         RI,par[24],parint[24],daylng;
	 float AutoIrrigAmt;
	 int   AutoIrrigateF;

 };
 //grid
 struct GridCommon{
	     int     NumNP, NumEl, IJ, KAT, MBand,Nmat, KX[4][NumElD];
		 float   x[NumNPD], y[NumNPD], Area[NumElD], nodeArea[NumNPD];
 };

//nodal
 struct NodeCommon{
	        int NumSol,NumG,ListN[NumNPD],ListNE[NumNPD],MatNumN[NumNPD];
			float hNew[NumNPD], ThNew[NumNPD], Vx[NumNPD], Vz[NumNPD], 
                  Q[NumNPD], Conc[NumSD][NumNPD], g[NumGD][NumNPD], 
                  Tmpr[NumNPD],Con[NumNPD],TcsXX[NumNPD],RO[NumNPD],
                  hNew_org[NumNPD], QAct[NumNPD],ThetaAvailRZ, ThetaFullRZ;
			float ThAvail[NumNPD],ThFull[NMatD];
			bool  lOrt;
			float QGas[NumGD][NumNPD], ThetaAir[NumNPD];
 };

//elements
 struct ElementCommon{
	      
	       float Sink[NumNPD], cSink[NumSD][NumNPD],
                 gSink[NumGD][NumNPD],tSink[NumNPD], 
                 RTWT[NumNPD], 
			     RMassM[NumNPD],RDenM[NumNPD],
			     RMassY[NumNPD], RDenY[NumNPD] ;
		    int    MatNumE[NumElD];
			float gSink_OM[NumGD][NumNPD], cSink_OM[NumGD][NumNPD],
				gSink_root[NumGD][NumNPD], gSink_rootY[NumGD][NumNPD],
				gSink_rootM[NumGD][NumNPD], gSink_N2O[NumGD][NumNPD];
 };

 //boundary
 struct  BoundaryCommon{
	       int  NumBP, NSurf, NVarBW,NVarBS,NVarBT,NVarBG,
                NumSurfDat, NSeep, NSP[NSeepD], NP[NumSPD][NSeepD],
				NDrain,NDR[NDrainD],ND[NumDR][NDrainD],
                KXB[NumBPD];
		   int CodeW[NumNPD],CodeS[NumNPD],CodeT[NumNPD],
			     CodeG[NumNPD],PCodeW[NumNPD];
		   float Width[NumBPD], VarBW[3][NumBPD],
			     VarBS[NumSD][NumBPD],
				 VarBT[4][NumBPD], VarBG[3][NumGD][NumBPD],EO,Tpot, gridWidth;
		   int PondingbyHead, PondingbyFlux;
 };
// Time

 struct TimeCommon{
	         double tNext[NumModD],dtMx[4],Time, Step, dtOpt,dtMin, 
				     dMul1, dMul2,tTDB[4],Tfin,tAtm;
			 float  Tinit;
			 int    lInput,Iter;
			 int   DailyOutput, HourlyOutput,RunFlag, DailyWeather, HourlyWeather;
			 int   beginDay, sowingDay, endDay, OutputSoilNo, OutPutSoilYes, Year;
		     int iTime,iDawn,iDusk;
			 double TimeStep;
 };


 //modules
 struct   ModuleCommon{
	              int NumMod,Movers[4], NShoot;
 };
 struct   ErrorCommon{
	              int errPlant;
 };

 struct FileCommon{
	           double starter;
	           char WeatherFile[256], TimeFile[256], BiologyFile[256],
               ClimateFile[256], NitrogenFile[256], SoluteFile[256],
               ParamGasFile[256],SoilFile[256], 
               ManagementFile[256],IrrigationFile[256],DripFile[256],
               WaterFile[256], WaterBoundaryFile[256], 
               GraphicsFile[256], InitialsFile[256],VarietyFile[256],
			   NodeGraphics[256],ElemGraphics[256],NodeGeomFile[256],
			   GeometryFile[256], SurfaceGraphics[256],
			   FluxGraphics[256], MasssBalanceFile[256],
			   MassBalanceFileOut[256],LeafFileIn[256],
			   OrganicMatterGraphics[256],
			   RunFile[256], MassBalanceRunoffFileOut[256],
			   MulchFile[256],MassBalanceMulchFileOut[256];
 };






#pragma pack()


 
#ifdef __cplusplus
extern "C" {
#endif

// Your exported function headers go here 
#ifdef _WIN32
PLANT_API void _stdcall CROP(struct ShootCommon    *, WeathCommon    *,
#else
PLANT_API void crop(struct ShootCommon    *, WeathCommon    *,
#endif
							         GridCommon     *, NodeCommon     *,
									 ElementCommon  *, BoundaryCommon *,
									 TimeCommon     *, ModuleCommon   *,
									 FileCommon     *);
#ifdef __cplusplus
}
#endif


 

