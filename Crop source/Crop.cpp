// myplant.cpp : Defines the entry point for the DLL application.
//#define MYPLANT_EXPORTS

#include "stdafx.h"
#include "crop.h"
#include "controller.h"

#include "time.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;
#include <cmath>
#define endl "\n"
#define comma ","

			 // note that we have to dereference the variable in order to
			 // assign a value that can be passed back to 2DSOIL. This is 
			 // because the FORTRAN program expects a pointer rather than
			 // a value. I don't think this applies to structures as it does 
            // variables that may be in the arguments list.
          // note use of lower case names. Upper and lower case conversions between
           // fortran and C++ don't matter here because these are arguments and not
           // a function name. CROP must be upper case because it is a function
           // name
int compare(const void *arg1, const void *arg2)
{
	/* Compare all of both strings: */
	if (*(double*)arg1 > *(double*)arg2) return 1;
	else if (*(double*)arg1 < *(double*)arg2) return -1;
	return 0;
};
#ifdef _WIN32
void _stdcall CROP(struct 
#else
void crop(struct 
#endif
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
	

	
	static CController* pSC; //SK, declare as static to ensure only one copy is instantiated during 2DSOIL execution
	// varFile contains variety information, GraphicFile holds output,LeafFile holds individual leaves

	if (time_public->lInput==1) //SK, initialiing crop module
	{
// Parse the file names from the FORTRAN strings passed from 2dsoil
//KY looks like GNU Fortran handle linebreak differently, making filename detection unusable
//KY this new macro based on std::string should work on both platforms with smaller code
#define SETSTR(s, n) std::string s(n, sizeof(n)); s.erase(s.find_last_not_of(" \n\r\t")+1);
		SETSTR(varFile, file_public->VarietyFile);
		SETSTR(GraphicFile, file_public->GraphicsFile);
		SETSTR(LeafFile, file_public->LeafFileIn);

// DT 10/07/2014 moved reading of initials file to soil model (Init.for)
		
		initInfo.plantDensity=SHOOTR->PopArea; 
        initInfo.latitude=Weather->LATUDE;
		initInfo.longitude=Weather->Longitude;
		initInfo.altitude=Weather->Altitude;
	    initInfo.year=time_public->Year;
		initInfo.sowingDay=time_public->sowingDay;
		initInfo.beginDay=time_public->beginDay;
		initInfo.endDay=time_public->endDay;
		initInfo.timeStep=time_public->TimeStep;
		initInfo.CO2 = Weather->CO2;
    	time_public->iTime=1;
		SHOOTR->LCAI=0.0;
		SHOOTR->LAREAT=0.0;
		SHOOTR->Height=0.0;
		//dt change for debugging purposes
		SHOOTR->Convr=1.0; // was 0.1 or 0.38 should be 1.0 as dry matter is used in all cases
		SHOOTR->AWUPS = 0.0;  //initialize AWUPS, AWUPS_old and LeafWP in 2DSOIL Yang 8/15/06
		SHOOTR->LeafWP = -0.5; 
		SHOOTR->PCRS = 0.0;
		SHOOTR->ET_demand = 0.0;
		SHOOTR->HourlyCarboUsed=0;  //it is also zero'd upon initialization in 2dsoil
		Period=time_public->TimeStep/60.0/24.0; // period should be in days, input in minutes
		PopSlab=SHOOTR->PopRow/100.0*SHOOTR->EOMult;
		SHOOTR->isEmerged=SHOOTR->isEmerged=0;
		/*These two lines show the relationships among some of the space variables.
              --PopSlab=SHOOTR->PopRow/100*SHOOTR->RowSp*SHOOTR->EOMult;
		      --PlantDensity=SHOOTR->PopRow*100.0/SHOOTR->RowSp;
	     */

// A new plant model object is created and initialized (calls initialize function) here
//  ***************************************************************************
 		pSC = new CController(varFile.c_str(), GraphicFile.c_str(), LeafFile.c_str(), initInfo); //Consider putting file names in initInfo
//  ***************************************************************************

		NitrogenUptake=pSC->getPlant()->get_TotalN()*PopSlab; //initialize nitrogen uptake accumulator with what is already in the plant
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

	//SK:  Running the crop module step by step
	if (module_public->NShoot>0)
	{
		WaterUptake=WaterUptake+SHOOTR->AWUPS*time_public->Step;  //g water per slab taken up to this time
		NitrogenUptake=NitrogenUptake+SHOOTR->SIncrSink/1.0e6; //Cumulative N (mass, g per slab in this time step 
		                                              //- 1.0e6 converts from ug to g) time step is accounted for in
		                                              //the solute uptake module

	}
// Note that SIncrSink has been multiplied by time step in the solute uptake routing
// the 1.0e6 scales from ug to g.


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
				wthr.daytime = time_public->Time;
				wthr.CO2 = Weather->CO2; 
				//if (Weather->CO2<=0) 
				//{
				//	wthr.CO2=initInfo.CO2;         //Can set CO2 in initials for specific simulations where CO2 is constant
				//}
				wthr.airT = Weather->TAIR[time_public-> iTime-1];
				wthr.canopyT = pSC->getPlant()->get_tmpr();
				wthr.PFD = Weather->par[time_public->iTime-1]*4.6; // conversion from PAR in Wm-2 to umol s-1 m-2
				wthr.solRad = Weather->WATTSM[time_public->iTime-1]; //one hour Total Radiation incident at soil surface (Wm-2)
				Es = (0.611*exp(17.502*wthr.airT/(240.97+wthr.airT))); // saturated vapor pressure at airT
				wthr.RH = (1-(Weather->VPD[time_public->iTime-1]/Es))*100.0; // relative humidity in percent
				wthr.rain = Weather->RINT[time_public->iTime-1];
				wthr.wind = Weather->WIND*(1000.0/3600.0); // conversion from km hr-1 to m s-1
				wthr.dayLength = Weather->daylng;
				wthr.LeafWP = SHOOTR->LeafWP/10;  //and leaf water potential information into MAIZESIM Yang 8/15/06 MPa
				wthr.pcrl=SHOOTR->PCRL/PopSlab/24.;
				wthr.pcrq=SHOOTR->PCRQ/PopSlab/24.;
				//since LeafWP in 2dsoil is in bar but in maizesim is in MPa, so, have to
				//divide it by 10 to convert it into MPa before passing the value to Maizesim 1 bar=10kPa

				if (abs(wthr.time-0.2083)<0.0001) //If time is 5 am, then pass the leaf water potential (the predawn leaf water potential)
					//from SHOOTR to the wthr object. YY

				{
					wthr.PredawnLWP=SHOOTR->LeafWP; //Here LeafWP is in bar. Since the LWPeffect in leaf.cpp uses leaf water potential
					//in bar, so here PredawnLWP is in bar, instead of being scaled to MPa. YY
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
				//then dividing it by plant density, converts it to g/hour/plant, 
				//which is the unit of the wthr.pcrs in maizesim. Yang. 10/27/06

				//Pass through nitrogen uptake (total mg per slab in the one hour) from 2DSOIL. 
				wthr.TotalRootWeight=SHOOTR->TotalRootWeight/PopSlab;
				wthr.MaxRootDepth=SHOOTR->MaxRootDepth;
				// Available water is cm per profile - should be divided by PopSlab
				wthr.ThetaAvailRZ=node_public->ThetaAvailRZ/PopSlab;
				
				
				if (NitrogenUptake >0 ) 
				{
					double nuptake = NitrogenUptake/PopSlab; //uptake is now grams per plant
					double leafloss = pSC->getPlant()->get_droppedLfArea();
					double nloss = leafloss*SLNmin;
					
					//SK 8/20/10: Here seems to be the only place where totalN of the plant is set. NitrogenUptake is initiated from get_N at the begining of the timestep so OK. 


                    pSC->getPlant()->set_TotalN(NitrogenUptake/PopSlab);  // Units are converted from g slab-1 to g plant -1 YY
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
					  ET_diff=wthr.ET_supply*24.0-SHOOTR->ET_demand;

				}
			}
			int end=grid_public->NumNP;
			double soilT=0, maxY=0;
			double soilMP[400];
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
					soilMP[i] = node_public->hNew[i];
					count++;
				}
			}

			wthr.soilT=soilT/count;
			std::qsort((void *)soilMP, count, sizeof(double),compare);
			wthr.SoilMP_med = soilMP[count / 2];

			
// The model code to simulate growth ect begins here when the plant object is called :
//TODO add some error catching code here
			int ier = pSC->getErrStatus();
			if ( ier == 0 ) 
			{
				ier = pSC->run(wthr); //Pass weather  into the "run" function
				                         //of the controller pSC YY
				
			}



			{
// Assumes that germination takes place about halfway through the sowing date
				if (wthr.time >= 0.49 && wthr.time <= 0.51)
				{
					if (!pSC->getPlant()->get_develop()->Germinated())
					{
						cout << "Germinating:" <<wthr.jday <<endl;
					} 
				}
			}
// The remaining groups of code handle carbon and nitrogen exchange between 2dsoil and maizsim
            if ((pSC->getPlant()->get_develop()->Germinated()) && (!pSC->getPlant()->get_develop()->Emerged()))
				// begin root growth at germination
			{
				if (!SHOOTR->isGerminated)
				{
					SHOOTR->isGerminated=1;
					SHOOTR->InitialRootCarbo=pSC->getPlant()->get_rootMass()*PopSlab; // get initial root mass to distribute over initial nodes
				}

				
			}
			if (pSC->getPlant()->get_develop()->Emerged())
				// pass appropriate data to 2DSOIL file structures 
			{
				if (!SHOOTR->isEmerged) SHOOTR->isEmerged=1;
				//ActualCarboIncrement is calculated from "assimilate", which in turn is calculated from photosynthsis_net in
				//plant; c_pool is 'assimilate' that has been converted from mass CO2 to mass CH2O
				//  the unit of c_pool then is in g Carbo/plant/hour, thus, at this point, pcrl has unit g/plant/hour
				// Multiply by 24 to get g plant-1 day-1; multiply by popslab to get g Carbo slab-1 day-1
				// dt 03/14/2011- I added a pool of carbo to hold leftover carbon from root growth, here it is implemented - see also plant
                double pool=pSC->getPlant()->get_C_pool_root(); //This holds any carbon not used for root growth in the previous time step
				if ((pSC->getPlant()->get_C_pool_root()>0) && (pSC->getPlant()->get_rootPart()<0.00001)) // this assures the pool is only used at night
					                                                                                    // minimizes complexity when pcrq has a value
																										// since we have to add leftover carbo from pcrq to the shoot
				{ 
					
					SHOOTR->PCRL=(pSC->getPlant()->get_rootPart()+pool)*24.0*PopSlab;
		    		pSC->getPlant()->set_C_pool_root(0.0);
				}

				else
				{
					SHOOTR->PCRL=(pSC->getPlant()->get_rootPart())*24.0*PopSlab;
					
				}
				bool gf = pSC->getPlant()->get_develop()->GrainFillBegan();
							
				SHOOTR->PCRQ=(pSC->getPlant()->get_rootPart()+ (pSC->getPlant()->get_shootPart()))*24.0*PopSlab;
                  if (gf)
				  {
                      SHOOTR->PCRQ=(pSC->getPlant()->get_rootPart())*24.0*PopSlab +
						  (0.75*pSC->getPlant()->get_shootPart())*24.0*PopSlab;
				  }
                //DT 09/19/14 under strong water stress mid season too much carbon is allocated to the roots, we
				// try to limit it here.
				//SHOOTR->PCRQ=SHOOTR->PCRL; //for debugging now remove later
			//dt 03/2011 added these two for debugging now - need to calculate mass balcance of carbo sent to root
				//can drop them later
				wthr.pcrl=  SHOOTR->PCRL/PopSlab/24.0;
				wthr.pcrq=  SHOOTR->PCRQ/PopSlab/24.0;

				
				SHOOTR->LCAI =pSC->getPlant()-> calcGreenLeafArea()* pSC->getInitInfo().plantDensity/(100.0*100.0);
				SHOOTR->Cover= 1.0 - exp (-0.79*SHOOTR->LCAI);
				SHOOTR->Shade=(float)SHOOTR->Cover*SHOOTR->RowSp*SHOOTR->EOMult;
				SHOOTR->Height=min(SHOOTR->Shade,SHOOTR->RowSp);
				SHOOTR->ET_demand = (pSC->getPlant()->get_ET()*24.0);//pass ET demand from shoot to root. Yang
				/*In GasExchange, the unit of ET is mmol m-2(leaf) sec-1
				
				need to convert to grams plant-1
				Here, multiplying ET by 0.018 and 3600*24 converts it to g m-2(ground) day-1
				dividing it by plantdensity converts it to g plant-1 day-1 */
				SHOOTR->LAI=(float)pSC->getPlant()->calcGreenLeafArea()*(float)pSC->getInitInfo().plantDensity/(100.0*100.0);
				//Pass LAI from maizesim into 2dsoil
                
				shoot_weightPerM2 = pSC->getPlant()->get_shootMass()*pSC->getInitInfo().plantDensity; //Calculate total shoot mass per meter aquared YY
                massIncrease = (shoot_weightPerM2 - old_shoot_weightPerM2); //Calculated increase in above-ground biomass per m2 YY
				massIncrease = max(0.0,massIncrease); // was going zero 2/6/2016 DT
				//The biomass  returned by getPlant()->get_shootMass() is the weight of each single plant (g/plant), 
				//to convert it into (g/m-2), it has to be 
				//multiplied by pSC->getIniInfo().plantDensity

				/*
				SK 8/20/10: Perhaps it's good idea to merge all plant N business into plant.cpp where C allocation is taken care of.
				For now I'm not modifying any codes here. TODO - dt 03/15/2011 Much of this nitrogen stuff can be cleaned up as some of
				Yang's original code is no longer used (i.e., U_N etc).

				*/

             	double NitrogenRatio;     //optimal N ratio according to N Dilution ratio as g N/ g dry matter
				double p_func;
				if (shoot_weightPerM2<100)
				{
					NitrogenRatio = N_min/100.0; //when shoot weight is lower than 100 g m-2, then nitrogen concentration is assumed to by .0410 g/g
				}
				else
				{
					p_func = pow(shoot_weightPerM2, (- N_shape));
                    NitrogenRatio = N_min *pow(shoot_weightPerM2,(-N_shape))/10.0; //sqrt(shoot_weightPerM2); //Calcualte above ground potential N concentration 
				//concentration as a function of aboveground biomass (Greenwood et al., 1990; Lindquist et al., 2007) YY
				}
                pSC->getPlant()->set_NitrogenRatio(NitrogenRatio);
			  //	double d=075;  //d: shape coefficient in the logistic function to simulate cumulative N uptake (Equation 9 in Lindquist et al. 2007)
				U_N = 0.359*d/4.0; //U_N maximum observed N uptake rate (g N m-2 ground d-1) (Lindquist et al, 2007) YY
			    //The unit of U_N is g N m-2 ground d-1

				U_M = q_n*massIncrease; //U_M maximum uptake rate as limited by maximum N fraction per unit (Equation 2 in Lindquist et al., 2007)
				//	double q_n = 0.032; //q_n the maximum ratio of daily N uptake to measured daily growth rate (g N g-1) (Lindquist et al., 2007)
                //unit of U_M is also g N m-2 ground d-1; however, the unit of massIncrease is g m-2/step (here one hour). 
				//Here in the simulation, the default length of one step is an hour; so, we have to scale it up to one
				//day by multiplying it by 24

				if (shoot_weightPerM2<100.0) //if shoot weight<100 (g m-2) then U_P is calculated this way 
				{
			
					U_P = (N_min/100.0)*massIncrease; // U_P potential rate of N accumulation (g N m-2 ground d-1) (Lindquist et al. 2007)
				}
				else //otherwise, it is calculated like this (Equation 6, Lindquist et al., 2007) YY
				{
					U_P = ((1.0-N_shape)*N_min/10.0)*pow(shoot_weightPerM2,-N_shape)*massIncrease;
				}
                //unit of U_P is also g N m-2 ground d-1; however, the unit of massIncrease is g/step. Here
				//in the simulation, the default length of one step is an hour; so, we have to scale it up to one
				//day by multiplying it by 24
				double t1 = N_min * 10.0/100.0 * pow(shoot_weightPerM2, -N_shape) * shoot_weightPerM2;
				double t2 = pSC->getPlant()->get_TotalN() * pSC->getInitInfo().plantDensity;
				//U_D = t1 - t2;
				double t3 = pSC->getPlant()->get_TotalN();
				U_D = N_min*10.0/100.0*pow(shoot_weightPerM2,-N_shape)*shoot_weightPerM2-pSC->getPlant()->get_TotalN()*pSC->getInitInfo().plantDensity;
				U_D=max(U_D,0.0)/24.0;
				//Since U_D is the difference between potential and actual amount of N in existing biomass
				// potentially the plant must take this up in a day if it could to meet the differenc. Thus it has implicit units of 
				// per day
				//U_D U uptake rate (g N m-2 d-1) as limited by the difference between potential and actual amount of N 
				//in existing biomass, equation 3 in Lindquist et al. 2007)
				//the returned value from get_TotalN() is in g N/plant. It has to be converted to g/m-2 ground
				//that's why the actual n content is mulitpled by pSC->getIniInfo().plantDensity YY

				// set up accounting of N here.
				// first hourly//Actual and needed N uptake in the last hour per plant per day
				//U_P is g N m-2, 
				// NitrogenUptake is g per slab, HourlyActualNFromSoil is grams per plant need to convert to g N m-2 to 
				//be consistent with other U_# variables.
				double HourlyActualNFromSoil = (NitrogenUptake - NitrogenUptakeOld)/PopSlab * pSC->getInitInfo().plantDensity*24.0;
				double HourlyNitrogenDemand= max(min(U_P, U_M,U_N, U_D),0.0); //Determine the nitrogen demand (equation 1 Lindquist et al. 2007) in grams plant-1
                pSC->getPlant()->set_HourlyNitrogenSoilUptake(HourlyActualNFromSoil);
				pSC->getPlant()->set_HourlyNitrogenDemand(HourlyNitrogenDemand);

                // now do cumulative amounts
				CumulativeNitrogenDemand+=HourlyNitrogenDemand; // grams m-2 day-1
				// send these to plant model as grams plant-1 day-1
                pSC->getPlant()->set_CumulativeNitrogenDemand(CumulativeNitrogenDemand/pSC->getInitInfo().plantDensity); //units are g plant day-1  
				pSC->getPlant()->set_CumulativeNitrogenSoilUptake(NitrogenUptake/PopSlab);


				double OldNDemand;
				OldNDemand = SHOOTR->nitroDemand / PopSlab /1.0e6;
				SHOOTR->nitroDemand = (float)HourlyNitrogenDemand*(float)PopSlab*1.0e6 * 24.0; //Pass the nitrogen demand into 2dsoil YY
				                                                               //Units are ug slab-1
				old_shoot_weightPerM2 = shoot_weightPerM2; //Save the value of the above_ground biomass of this time-step
				NitrogenUptakeOld=NitrogenUptake; // save the cumulative N uptake from this time step;
				SHOOTR->NDemandError=(float)CurrentNUptakeError;
				SHOOTR->CumulativeNDemandError=(float)CumulativeNUptakeError;

			}

			if (pSC->getPlant()->get_develop()->Dead()) 
//			if (pSC->getLastDayOfSim() <= pSC->getTime()->get_day_of_year()) 
			{
				cout << "Completing crop simulation..." <<endl;
				module_public->NShoot=0; //tell 2dsoil that crops harvested
				time_public->tNext[ModNum-1]=1.0e12; // set the next time so the model
				pSC = NULL; // if matured points to nothing
				delete pSC;
				time_public->RunFlag=0;
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