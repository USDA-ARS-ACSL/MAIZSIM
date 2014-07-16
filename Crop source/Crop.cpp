// myplant.cpp : Defines the entry point for the DLL application.
//#define MYPLANT_EXPORTS

#include "stdafx.h"
#include "crop.h"
#include "controller.h"

#include "time.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;
#include <cmath>
#define endl "\n"

			 // note that we have to dereference the variable in order to
			 // assign a value that can be passed back to 2DSOIL. This is 
			 // because the FORTRAN program expects a pointer rather than
			 // a value. I don't think this applies to structures as it does 
            // variables that may be in the arguments list.
          // note use of lower case names. Upper and lower case conversions between
           // fortran and C++ don't matter here because these are arguments and not
           // a function name. CROP must be upper case because it is a function
           // name
void _stdcall CROP(struct 
				   ShootCommon        *SHOOTR, 
				   WeathCommon        *Weather,
				   GridCommon         *grid_public,
				   NodeCommon         *node_public,
				   ElementCommon      *ele_public,
				   BoundaryCommon     *bound_public,
				   TimeCommon         *time_public,
				   ModuleCommon       *module_public,
				   FileCommon         *file_public
				   )

{

	// I think pLeaf can be local since the leaf array holds the list
	// and the pointers are not lost between invocations of the procedure
	//First read input data if start of simulation
	char* Buffer=(char*)calloc(256,sizeof(char));
	
	// DEL LATER ::NR dConvert;
	ofstream ostr;
	static CController* pSC; //SK, declare as static to ensure only one copy is instantiated during 2DSOIL execution
	static char varFile[132], GraphicFile[132], LeafFile[132]; //need initials file name too?
	// varFile contains variety information, GraphicFile holds output,LeafFile holds individual leaves

	if (time_public->lInput==1) //SK, initialiing crop module
	{
// Parse the file names from the FORTRAN strings passed from 2dsoil
		int size=sizeof(file_public->VarietyFile);
		int i;
		for (i=0 ;i<size;i++)
		{
			if (file_public->VarietyFile[i]!='\n')
			{
				varFile[i]=file_public->VarietyFile[i];
			}
		}
		varFile[i]='\0'; // append null char to end of string

		size=sizeof(file_public->GraphicsFile);
		for (int i=0 ;i<size;i++)
		{
			if (file_public->GraphicsFile[i]!='\n')
			{
				GraphicFile[i]=file_public->GraphicsFile[i];
			}
		}
		GraphicFile[i]='\0';
   
	    size=sizeof(file_public->LeafFileIn);
		for (i=0 ;i<size;i++)
		{
			if (file_public->LeafFileIn[i]!='\n')
			{
				LeafFile[i]=file_public->LeafFileIn[i];
			}
		}
		LeafFile[i]='\0'; // append null char to end of string


		std::cout << "Crop object reading initials file..." << endl;
		size=sizeof(file_public->InitialsFile);
		char iniFile[132]; //only read here so we can create a temporary variable
		for (i=0; i<size; i++)
		{
			if (file_public->InitialsFile[i]!='\n')
			{
				iniFile[i]=file_public->InitialsFile[i];
			}
		}
		iniFile[i]='\0'; 
		//  Read initial file.

		ifstream in(iniFile);
		if (!in) 
		{
			cout << "error opening file" << endl ; 
			// acctually need to call the 2DSOIL error message
			// routine
			//			 error_public->errPlant=1;
			return;
		}
		CString Junk;
		char ttmp[10];
//DT for some reason need to read the buffer twice when reading text after a line of numbers
// TODO - Plant density can be calculated from poprow and row spacing.
// TODO need to add the initinfo components into SHOOTR and read from there. 
		in.getline(Buffer,100) ;
		in.getline(Buffer,100) ;
		in >> SHOOTR->PopRow >> SHOOTR->RowSp >>initInfo.plantDensity >> SHOOTR->RowAng ;
		in >> SHOOTR->xBStem >> SHOOTR->yBStem >> SHOOTR->CEC >>SHOOTR->EOMult >>initInfo.CO2;
		in.getline(Buffer,10);
	    in.getline(Buffer,100);
		in >> initInfo.latitude >> initInfo.longitude >> initInfo.altitude;
		in.getline(Buffer,10);
		in.getline(Buffer,100);
		in >> Weather->AutoIrrigate;
		in.getline(Buffer,10);
		in.getline(Buffer,100);

        
		int mm,dd,yy;
	
		Timer  dConvert;

		in >> ttmp;
		sscanf(ttmp, "%d/%d/%d", &mm,&dd,&yy);
		initInfo.beginDay = dConvert.julday(mm,dd,yy);
		initInfo.year=yy;


		in >> ttmp;
		sscanf(ttmp, "%d/%d/%d", &mm,&dd,&yy);
		initInfo.sowingDay = dConvert.julday(mm,dd,yy);

		in >> ttmp;
		sscanf(ttmp, "%d/%d/%d", &mm,&dd,&yy);
		Emergence = dConvert.julday(mm,dd,yy);
		
		in >> ttmp;
		sscanf(ttmp, "%d/%d/%d", &mm,&dd,&yy);
		initInfo.endDay = dConvert.julday(mm,dd,yy);
        in >> initInfo.timeStep;


		if (in.gcount() <10)
		{
			//			 errPlant=1;
			cout << "not enough data" << endl; // rudimentary error checking
		}
		in.close();  
		


		SHOOTR->iTime=1;
		SHOOTR->LCAI=0.0;
		SHOOTR->LAREAT=0.0;
		SHOOTR->Height=0.0;
		//dt change for debugging purposes
		SHOOTR->Convr=1.0; // was 0.1 or 0.38 should be 1.0 as dry matter is used in all cases
		SHOOTR->AWUPS = 0.0;  //initialize AWUPS, AWUPS_old and psil_ in 2DSOIL Yang 8/15/06
		SHOOTR->psil_ = -0.5; 
		SHOOTR->PCRS = 0.0;
		SHOOTR->ET_demand = 0.0;
		SHOOTR->HourlyCarboUsed=0;  //it is also zero'd upon initialization in 2dsoil
		Period=1.0/24.0;
		PopSlab=SHOOTR->PopRow/100.0*SHOOTR->EOMult;
            //PopSlab=SHOOTR->PopRow/100*SHOOTR->RowSp*SHOOTR->EOMult;
		Popare=SHOOTR->PopRow*100.0/SHOOTR->RowSp;

// A new plant model object is created and initialized (calls initialize function) here
//  ***************************************************************************
 		pSC = new CController(varFile, GraphicFile, LeafFile, initInfo);
//  ***************************************************************************

		//can pass root parameters back to 2DSOIL here
		/*
		This block is not used now, all 2dsoil variables are read by 2dsoil 
		This should be deleted after we are sure of the structure
        SHOOTR->RRRM=initInfo.RRRM;
		SHOOTR->RRRY=initInfo.RRRY;
		SHOOTR->RVRL=initInfo.RVRL;
		SHOOTR->ALPM=initInfo.ALPM;
		SHOOTR->ALPY=initInfo.ALPY;
        SHOOTR->RTWL=initInfo.RTWL;
		SHOOTR->RtMinWtPerUnitArea=initInfo.RtMinWtPerUnitArea;
        SHOOTR->Wl=initInfo.Wl;
		SHOOTR->Wa=initInfo.Wa;
	    SHOOTR->Wr=initInfo.Wr;
		SHOOTR->Wb=initInfo.Wb;
		*/
// -----------
		NitrogenUptake=pSC->getPlant()->get_N()*PopSlab; //initialize nitrogen uptake with what is already in the plant
		//SK 8/20/10: this is curious but OK

		SHOOTR->NDemandError=0;
		SHOOTR->CumulativeNDemandError=0;
		time_public->RunFlag=1;

		{
		module_public->NumMod=module_public->NumMod+1 ;
		ModNum=module_public->NumMod;
		
		time_public->tNext[ModNum-1]=pSC->getSowingDay();
	} // end if
	} //end initialization


//	if((module_public->NShoot == 0) && (fabs(time_public->Time-pSC->getSowingDay()))<0.001)
		//	{
//		module_public->NumMod=module_public->NumMod+1 ;
//		ModNum=module_public->NumMod;
//		module_public->NShoot=1;
//		time_public->tNext[ModNum-1]=time_public->Time+Period;
//	} // end if
	//SK:  Running the crop module step by step
	if (module_public->NShoot>0)
	{
		WaterUptake=WaterUptake+SHOOTR->AWUPS*time_public->Step;  //g water per slab taken up in an hour
		NitrogenUptake=NitrogenUptake+SHOOTR->SIncrSink/1.0e6; //Cumulative N (mass, g (plant slab)-1) in this time step
	}
		
// Note that SIncrSink has been multiplied by time step in the solute uptake routing
// the 1000 scales from ug to mg.
		if(fabs(time_public->Time-time_public->tNext[ModNum-1])< fabs(0.001*time_public->Step))
		{
		//If the sowing date has come and there is not plant, let the program know so other calculations are not done
			if((module_public->NShoot == 0) && (fabs(time_public->Time-pSC->getSowingDay()))<0.001)
				{
				   module_public->NShoot=1;
				}

			//DT added soil temperature from around seed area (right now it is taken as default - 20)
			// find soil temperature of surface 5 cm
			double Es;
            //WaterUptake=WaterUptake*24;   //Water uptake rate per day per slab
			// calculate error for demand and actual uptake, if negative, demand is greater then uptake
			CurrentNUptakeError=NitrogenUptake/PopSlab-pSC->getPlant()->get_CumulativeNitrogenDemand();
			CumulativeNUptakeError+=CurrentNUptakeError;

			TWeather wthr;
			{
			    wthr.HourlyOutput= time_public->HourlyOutput;
				wthr.DailyOutput=time_public->DailyOutput;
				wthr.jday = Weather->JDAY;
				wthr.time = time_public->Time-Weather->JDAY;
				wthr.CO2 = initInfo.CO2;
				wthr.airT = Weather->TAIR[SHOOTR->iTime-1];
				wthr.PFD = Weather->par[SHOOTR->iTime-1]*4.6; // conversion from PAR in Wm-2 to umol s-1 m-2
				wthr.solRad = Weather->WATTSM[SHOOTR->iTime-1]; //conversion from Wm-2 to J m-2 in one hour
				Es = (0.611*exp(17.502*wthr.airT/(240.97+wthr.airT))); // saturated vapor pressure at airT
				wthr.RH = (1-(Weather->VPD[SHOOTR->iTime-1]/Es))*100.0; // relative humidity in percent
				wthr.rain = Weather->RINT[SHOOTR->iTime-1];
				wthr.wind = Weather->WIND*(1000.0/3600.0); // conversion from km hr-1 to m s-1
				wthr.dayLength = Weather->daylng;
				wthr.psil_ = SHOOTR->psil_/10;  //and leaf water potential information into MAIZESIM Yang 8/15/06
				//since psil_ in 2dsoil is in bar but in maizesim is in MPa, so, have to
				//divide it by 10 to convert it into MPa before passing the value to Maizesim 1 bar=10kPa

				if (abs(wthr.time-0.2083)<0.0001) //If time is 5 am, then pass the leaf water potential (the predawn leaf water potential)
					//from SHOOTR to the wthr object. YY

				{
					lwpd=SHOOTR->psil_; //Here psil_ is in bar. Since the LWPeffect in leaf.cpp uses leaf water potential
					//in bar, so here lwpd is in bar, instead of being scaled to MPa. YY
				}

// pass actual carbohydrate amount used in 2dsoil back to the plant
				//ToDo - make pcrs a new variable (ActualRootCarboUsed) and make it a member of plant.
				//dt here I changed this temporarily for debugging
				//don't need to divide by 24 since the value has been integrated over an hour
				 wthr.pcrs = SHOOTR->HourlyCarboUsed/PopSlab;   //original
				 SHOOTR->HourlyCarboUsed=0.0; 
		 
				//dividing it by PopSlab converts it to g/day/plant;
				//ToDo: need to document this better, what is pcrs being used for.

				//SHOOTR->PCRS in 2dsoil is the actual rate of carbon supplied to roots in a soil slab, it is in g/day;
				//dividing it by (ShOOTR->Rowsp*1)/10000 converts it to g/day/m^2;
				//further dividing it by weather->daylng converts it to g/hour/m^2;
				//then dividing it by Popare, which is the plant density, converts it to g/hour/plant, 
				//which is the unit of the wthr.pcrs in maizesim. Yang. 10/27/06

				//Pass through nitrogen uptake (total mg per slab in the one hour) from 2DSOIL. 
				wthr.TotalRootWeight=SHOOTR->TotalRootWeight/PopSlab;
				wthr.MaxRootDepth=SHOOTR->MaxRootDepth;
				// Available water is cm per profile - should be divided by PopSlab
				wthr.ThetaAvail=node_public->ThetaAvail/PopSlab;
				
				
				if (NitrogenUptake >0 ) 
				{
					double nuptake = NitrogenUptake/PopSlab;
					double leafloss = pSC->getPlant()->get_droppedLfArea();
					double nloss = leafloss*SLNmin;
					
					//SK 8/20/10: Here seems to be the only place where totalN of the plant is set. NitrogenUptake is initiated from get_N at the begining of the timestep so OK. 


                    pSC->getPlant()->set_N(NitrogenUptake/PopSlab);  // Units are converted from g slab-1 to g plant -1 YY
// need to look at loss of N in the dropped leaf (plant N goes negative?)				                                                         
					//pSC->getPlant()->set_N(NitrogenUptake/PopSlab-pSC->getPlant()->get_droppedLfArea()*SLNmin);
				   
                   
				}
				 //SLNmin base Specific leaf nitrogen content; for now assume it's 0.5 YY
// debuging dt 1-30-2012
				
				if (SHOOTR->LAI==0) 
				{
					wthr.ET_supply=0;
				}
				else
				{
				    //Note water uptake has been summed over the past hour so it is an hourly amount
					  //into MAIZESIM Yang 8/15/06, dt 4/24/2011
				      wthr.ET_supply = WaterUptake/(SHOOTR->EOMult*SHOOTR->PopRow)*100; //units area gr per plant per hour
					//dt 4-24-2011 I replaced SHOOTR->AWUPS with WaterUptake. AWUPS is an instantaneous value.
					/* 
					ET_Supply is the actual amount of water that can be taken from the soil slab ( derived from AWUPS, g day-1 slab-1)
					To compare this variable with the et rate in maizesim it has to be converted into grams water per plant
					To do this multiply by EOMULT to double slab width if plant is at the edge. Then multiply by 100/PopRow
					to get area inhabited by the plant. This provides a per plant estimate from area.
					*/
					  //debug
					  ET_diff=wthr.ET_supply*24-SHOOTR->ET_demand;

				}
			}
			int end=grid_public->NumNP;
			double soilT=0, maxY=0;
			int count=0;
			double LowerBoundary=5;
			// first find top of grid
			for (int i=0; i<end; i++)
			{ 
				if (grid_public->y[i] > maxY) maxY=grid_public->y[i];
			}
			LowerBoundary=maxY-LowerBoundary;
			// now find average temperature in layer between surface and lower boundary

			for (int i=0; i<end; i++)
			{
				if ( grid_public->y[i]>=LowerBoundary)
				{
					soilT=soilT+node_public->Tmpr[i];
					count++;
				}
			}

			wthr.soilT=soilT/count;
// The model code to simulate growth ect begins here when the plant object is called :
//TODO add some error catching code here
			int ier = pSC->getErrStatus();
			if ( ier == 0 ) 
			{
				ier = pSC->run(wthr, lwpd); //Pass both weather and leaf water potential into the "run" function
				//of the controller pSC YY
				// need to get rid of other run module (with only wthr) done?
				// need to add lwpd to wthr structure
			}



			{
// Assumes that germination takes place about halfway through the sowing date
				if (wthr.time >= 0.49 && wthr.time <= 0.51)
				{
					if (!pSC->getPlant()->get_develop()->Germinated())
					{
						cout << wthr.jday <<endl;
					} 
				}
			}

			if (pSC->getPlant()->get_develop()->Emerged())
				// pass appropriate data to 2DSOIL file structures 
			{
				//ActualCarboIncrement is calculated from "assimilate", which in turn is calculated from photosynthsis_net in
				//plant; the unit of assimilate then is in g/plant/hour, thus, at this point, pcrl has unit g/plant/hour
				// Multiply by 24 to get g plant-1 day-1; multiply by popslab to get g Carbo slab-1 day-1
				// dt 03/14/2011- I added a pool of carbo to hold leftover carbon from root growth, here it is implemented - see also plant
                double pool=pSC->getPlant()->get_C_pool_root(); //This holds any carbon not used for root growth in the previous time step
				if ((pSC->getPlant()->get_C_pool_root()>0) && (pSC->getPlant()->get_rootPart()<0.00001)) // this assures the pool is only used at night
					                                                                                    // minimizes complexity when pcrq has a value
																										// since we have to add leftover carbo from pcrq to the shoot
				{ 
					
					SHOOTR->PCRL=(pSC->getPlant()->get_rootPart()+pool)*24*PopSlab;
		    		pSC->getPlant()->set_C_pool_root(0.0);
				}

				else
				{
					SHOOTR->PCRL=(pSC->getPlant()->get_rootPart())*24*PopSlab;
					
				}
                SHOOTR->PCRQ=(pSC->getPlant()->get_rootPart()+ pSC->getPlant()->get_shootPart())*24*PopSlab;
				//SHOOTR->PCRQ=SHOOTR->PCRL; //for debugging now remove later
			//dt 03/2011 added these two for debugging now - need to calculate mass balcance of carbo sent to root
				//can drop them later
				wthr.pcrl=  SHOOTR->PCRL/PopSlab/24;
				wthr.pcrq=  SHOOTR->PCRQ/PopSlab/24;

				
				SHOOTR->LCAI =pSC->getPlant()-> calcGreenLeafArea()* pSC->getInitInfo().plantDensity/(100*100);
				SHOOTR->Cover= 1.0 - exp (-0.79*SHOOTR->LCAI);
				SHOOTR->Shade=(float)SHOOTR->Cover*SHOOTR->RowSp;
				SHOOTR->Height=min(SHOOTR->Shade,SHOOTR->RowSp);
				SHOOTR->ET_demand = (pSC->getPlant()->get_ET()*24);//pass ET demand from shoot to root. Yang
				/*In GasExchange, the unit of ET is mmol m-2(leaf) sec-1
				
				need to convert to grams plant-1
				Here, multiplying ET by 0.018 and 3600*24 converts it to g m-2(ground) day-1
				dividing it by plantdensity converts it to g plant-1 day-1 */
				SHOOTR->LAI=(float)pSC->getPlant()->calcGreenLeafArea()*(float)pSC->getInitInfo().plantDensity/(100*100);
				//Pass LAI from maizesim into 2dsoil
                
				shoot_weightPerM2 = pSC->getPlant()->get_shootMass()*pSC->getInitInfo().plantDensity; //Calculate total shoot mass per meter aquared YY
                massIncrease = (shoot_weightPerM2 - old_shoot_weightPerM2); //Calculated increase in above-ground biomass per m2 YY
				
				//The biomass  returned by getPlant()->get_shootMass() is the weight of each single plant (g/plant), 
				//to convert it into (g/m-2), it has to be 
				//multiplied by pSC->getIniInfo().plantDensity

				/*
				SK 8/20/10: Perhaps it's good idea to merge all plant N business into plant.cpp where C allocation is taken care of.
				For now I'm not modifying any codes here. TODO - dt 03/15/2011 Much of this nitrogen stuff can be cleaned up as some of
				Yang's original code is no longer used (i.e., U_N etc).

				*/

             	double NitrogenRatio;     //optimal N ratio according to N Dilution ratio
				
				if (shoot_weightPerM2<100)
				{
					NitrogenRatio = a/100; //when shoot weight is lower than 100 g m-2, then nitrogen concentration is assumed to by .0410 g/g
				}
				else
				{
					
                    NitrogenRatio = a/100.0 *pow(shoot_weightPerM2,(1.0-b)); //sqrt(shoot_weightPerM2); //Calcualte above ground potential N concentration 
				//concentration as a function of aboveground biomass (Greenwood et al., 1990; Lindquist et al., 2007) YY
				}
                pSC->getPlant()->set_NitrogenRatio(NitrogenRatio/10.0);
			  //	double d=075;  //d: shape coefficient in the logistic function to simulate cumulative N uptake (Equation 9 in Lindquist et al. 2007)
				U_N = 0.359*d/4; //U_N maximum observed N uptake rate (g N m-2 ground d-1) (Lindquist et al, 2007) YY
			    //The unit of U_N is g N m-2 ground d-1

				U_M = q_n*massIncrease*24; //U_M maximum uptake rate as limited by maximum N fraction per unit (Equation 2 in Lindquist et al., 2007)
				//	double q_n = 0.032; //q_n the maximum ratio of daily N uptake to measured daily growth rate (g N g-1) (Lindquist et al., 2007)
                //unit of U_M is also g N m-2 ground d-1; however, the unit of massIncrease is g m-2/step (here one hour). 
				//Here in the simulation, the default length of one step is an hour; so, we have to scale it up to one
				//day by multiplying it by 24

				if (shoot_weightPerM2<100) //if shoot weight<100 (g m-2) then U_P is calculated this way 
				{
			
					U_P = (a/100.0)*massIncrease*24; // U_P potential rate of N accumulation (g N m-2 ground d-1) (Lindquist et al. 2007)
				}
				else //otherwise, it is calculated like this (Equation 6, Lindquist et al., 2007) YY
				{
					U_P = ((1.0-b)*a*10.0/100.0)*pow(shoot_weightPerM2,-b)*massIncrease*24;
				}
                //unit of U_P is also g N m-2 ground d-1; however, the unit of massIncrease is g/step. Here
				//in the simulation, the default length of one step is an hour; so, we have to scale it up to one
				//day by multiplying it by 24

				U_D = a*10/100*pow(shoot_weightPerM2,-b)-pSC->getPlant()->get_N()*(pSC->getInitInfo().plantDensity/(100*100));
				//U_D U uptake rate (g N m-2 d-1) as limited by the difference between potential and actual amount of N 
				//in existing biomass, equation 3 in Lindquist et al. 2007)
				//the returned value from get_N() is in g N/plant. It has to be converted to g/m-2 ground
				//that's why the actual n content is mulitpled by pSC->getIniInfo().plantDensity/(100*100) YY

				// set up accounting of N here.
				// first hourly//Actual and needed N uptake in the last hour per plant per day
                double HourlyActualNFromSoil =(NitrogenUptake-NitrogenUptakeOld)/PopSlab; //houly rate per day
				double HourlyNitrogenDemand= max(U_P, 0)/pSC->getInitInfo().plantDensity/24.0; //Determine the nitrogen demand (equation 1 Lindquist et al. 2007) in grams plant-1
                pSC->getPlant()->set_HourlyNitrogenSoilUptake(HourlyActualNFromSoil);
				pSC->getPlant()->set_HourlyNitrogenDemand(HourlyNitrogenDemand);

                // now do cumulative amounts
				CumulativeNitrogenDemand+=HourlyNitrogenDemand; // grams plant-1 day-1
                pSC->getPlant()->set_CumulativeNitrogenDemand(CumulativeNitrogenDemand); //units are g plant day-1  
				pSC->getPlant()->set_CumulativeNitrogenSoilUptake(NitrogenUptake/PopSlab);


				double OldNDemand;
				OldNDemand=SHOOTR->nitroDemand/PopSlab/1e6/24;
				SHOOTR->nitroDemand = (float)HourlyNitrogenDemand*(float)PopSlab*1e6*24; //Pass the nitrogen demand into 2dsoil YY
				                                                               //Units are ug slab-1
				old_shoot_weightPerM2 = shoot_weightPerM2; //Save the value of the above_ground biomass of this time-step
				NitrogenUptakeOld=NitrogenUptake; // save the cumulative N uptake from this time step;
				SHOOTR->NDemandError=(float)CurrentNUptakeError;
				SHOOTR->CumulativeNDemandError=(float)CumulativeNUptakeError;

			}

			if (pSC->getPlant()->get_develop()->Matured()) 
			{
				cout << "Completing crop simulation..." <<endl;
				module_public->NShoot=0; //tell 2dsoil that crops harvested
				time_public->tNext[ModNum-1]=1e12; // set the next time so the model
				pSC = NULL; // if matured points to nothing
				delete pSC;
				time_public->RunFlag=0;
				ostr.close();
			}
			else
			{
				time_public->tNext[ModNum-1]=time_public->Time+Period;
				WaterUptake=0;
				//NitrogenUptake=0;
			}
		}// end hourly calculations code
	 // end if NShhot>0 section 

	return;
}