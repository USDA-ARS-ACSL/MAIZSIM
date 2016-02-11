//
#include "stdafx.h"
#include "plant.h"
#include "gas_exchange.h"
#include "lightenv.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#define PRIMORDIA 5   //primordia initiated in the embryo, conservative species-wide, (Poethig, 1994)
#define CO2_MW 44.0098
#define C_MW 12.011
#define CH2O_MW 30.03

using namespace std;


CPlant::CPlant(const TInitInfo& info )
{
	nodalUnit = NULL;
	ear = NULL;
	roots = NULL;
	develop = NULL;
	seedMass = mass =  CH2O = 0.275; // seed weight g/seed
	C_content = 0.40; // 40% C, See Kim et al. (2007) EEB
	C_pool = 0.0;  //mass*C_content;
	C_reserve = 0.0;
	C_ReserveLeaf=0.0;
	C_demand = C_supply = 0.0;
	C_pool_root=0.0;
	maintRespiration = 0.0;
	C2_effect=1.0;
	SunlitRatio=0.0;

	// initialize plant part sizes //
	shootPart = 0.60;
	rootPart = 0.40;
	earMass=droppedLeafmass=rootMass=0.0;
	//rootMass = rootPart*seedMass; This will be implemented later when roots grow from germination
	;
	shootMass=seedMass*(1.0-rootPart);
	rootMass=seedMass-shootMass;  // need to distribute this in a few seed elements. at about 1.059e-4 wt/L =780 cm
	leafMass = activeLeafMass=0.90*shootMass;
	stemMass=0.10*shootMass;
	shootPart_old = shootPart; rootPart_old = rootPart;
    TotalNitrogen = 0.275*3.4/100.0; //assume nitrogen concentration at the beginning is 3.4% of the total weight of the seed
	// need to check with Yang. This doesn't look correct
    leaf_N = 0;
	leaf_NFraction = 0;
	N_pool = 0.0;
	CumulativeNitrogenDemand=0;
	HourlyNitrogenDemand=0;
	CumulativeNitrogenSoilUptake=0;



    SunlitLAI=0.0;
	ShadedLAI=0.0;
	sowingDay = 1.0;
	age = 0.0;
	initInfo = info;
	roots = new CRoots();
	ear = new CEar();
	develop = new CDevelopment(initInfo);
	nodalUnit = new CNodalUnit[initInfo.genericLeafNo+10]; // create enough leaf nodes for now, to be replaced by dynamic collection
	for (int i=0; i <= PRIMORDIA; i++) // leaf[0] is a coleoptile
	{
		nodalUnit[i].initialize(i, develop);
		nodalUnit[i].get_leaf()->set_mass(leafMass/PRIMORDIA); //assing initial mass for each primordium
       	nodeNumber = i;
	}
	finalNodeNumber = info.genericLeafNo;
	leafArea =greenLeafArea = actualGreenLeafArea = senescentLeafArea = potentialLeafArea = droppedLfArea=0.0;
	previousDroppedlfArea = currentDroppedLfArea = droppedLfArea = 0;
	
	temperature = develop->get_Tcur();
	
	photosynthesis_net =photosynthesis_gross = transpiration = transpirationOld=assimilate =  0.0;
		VPD=0;
    conductance=0;
	emerge_gdd = 0;
	SunlitRatio=0.0;
}

CPlant::~CPlant() 
{
	if (nodalUnit != NULL) delete [] nodalUnit;
	if (roots != NULL) delete roots;
	if (develop != NULL) delete develop;
	if (ear != NULL) delete ear;
}



void CPlant::update(const TWeather & weather, double PredawnLWP) 
{
	int TotalGrowingLeaves = 0;
	int TotalDroppedLeaves = 0;
	double PlantNitrogenContent;
	double percentN;
	if (develop->Emerged())
	{
		calcRed_FRedRatio(weather);
	}

	develop->update(weather);
	
	
	finalNodeNumber = develop->get_youngestLeaf();
	//SK: get N fraction allocated to leaves, this code is just moved from the end of the procedure, this may be taken out to become a separate fn
	{
	    //Calculate faction of nitrogen in leaves (leaf NFraction) as a function of thermal time from emergence
	    //Equation from Lindquist et al. 2007 YY
		//SK 08/20/10: TotalNitrogen doesn't seem to be updated at all anywhere else since initialized from the seed content
		//SK: I see this is set in crop.cpp ln 253 from NUptake from 2dsoil
		//but this appears to be the amount gained from the soil for the time step; so how does it represent totalNitrogen of a plant?

       double thermal_time = develop->get_GDDsum()-emerge_gdd;//record thermal time from emergency YY
	    //Calculate faction of nitrogen in leaves (leaf NFraction) as a function of thermal time from emergence
	    //Equation from Lindquist et al. 2007 YY

        leaf_NFraction = 0.79688-0.00023747*thermal_time-0.000000086145*thermal_time*thermal_time;
		if (leaf_NFraction <=0)
		{
			leaf_NFraction = 0;//fraction of leaf n in total shoot n can't be smaller than zero. YY
		}
		leaf_N = leaf_NFraction*this->get_N(); //calculate total nitrogen amount in the leaves YY units are grams N in all the leaves
		leaf_N_content = leaf_N/(this->greenLeafArea/10000); //SK 8/22/10: set avg greenleaf N content before update in g/m2; 
		PlantNitrogenContent=this->get_N();
		percentN=this->get_N()/shootMass;
		double Ratio=0;
		if (NitrogenRatio>0)
			{
		      Ratio= percentN/NitrogenRatio;
			}
	}

	
	if (!develop->Germinated())
	{
		temperature = develop->get_Tcur();
		return;
	}
	else if(!develop->Emerged()) //after germination but emergence
	{
// here we calculate the mass of roots that were initialized in the soil model (read in with the element data)
// This is so there is no discontinuity in root mass (relative to the carbon supplied by the plant later)
// Jan 9, 2015, added root growth before emergence 
				// these are roots taht grew from the seed
		if (!this->get_roots()->GetInitialized())
			{
				// now gets root mass from seed partitioning
				this->get_roots()->SetInitialized();
                this->get_roots()->import_CH2O(rootMass);
			    //rootMass=weather.TotalRootWeight;
			}

		 temperature = develop->get_Tcur();
	     C_reserve = seedMass*C_content;
		//				C_pool += C_reserve*(1/20)*(1/24)*(initInfo.timeStep/60); // assume it takes 20 days to exhaust seed C reserve 
		 C_pool = C_reserve;
 	for (int i = 1; i <= develop->get_LvsInitiated()  ; i++)
		{
			
 			if(!nodalUnit[i].isInitiated())
			{
				nodalUnit[i].initialize(i,develop);
                nodeNumber = i;
			}
			else
			{
				nodalUnit[i].update(develop,PredawnLWP);
				// from germination to emergence, C supply is from the seed
				if (nodalUnit[i].get_leaf()->isGrowing()) 
				 {
					 TotalGrowingLeaves++;
				}
				if (nodalUnit[i].get_leaf()->isDropped()) 
				{
					TotalDroppedLeaves++;
				}
					
			}
		}	
// calculate relative area increases for leaves now that they are updated

		//calcMaintRespiration(weather); // commented for now, have to test this for seedling
		nodalUnit[0].get_leaf()->set_TotalGrowingLeaves(TotalGrowingLeaves);
		nodalUnit[0].get_leaf()->set_TotalDroppedLeaves(TotalDroppedLeaves);

		if(nodalUnit[1].get_leaf()->isInitiated() && !nodalUnit[1].get_leaf()->isAppeared())
		{
			calcMaintRespiration(weather);
			C_allocation(weather);
			seedMass = __max(0.0,seedMass-C_supply); // need setmass
			setMass();
		}
		else if (nodalUnit[1].get_leaf()->isAppeared())
		{
	        calcPerLeafRelativeAreaIncrease();
			calcLeafArea();
			calcGasExchange(weather);
			double c_pool1 = C_pool;
			double c_pool2 = assimilate*CH2O_MW/CO2_MW; // convert from grams CO2 to grams carbohydrate (per hour per plant)
            C_pool += c_pool2;
			calcMaintRespiration(weather);
			C_allocation(weather);
			seedMass = __max(0.0,seedMass-C_supply*c_pool1/C_pool); // get seedmass reduction only for those C from the seed
			setMass();
		}
		return;
	}
	else if(!develop->Dead())
	{
		emerge_gdd = develop->get_EmergeGdd(); //if the plant emergies, then pass the emerge_gdd value from the develop object into the plant object

		for (int i = 1; i <= develop->get_LvsInitiated(); i++)
		{
			if(!nodalUnit[i].isInitiated())
			{
				nodalUnit[i].initialize(i,develop);
				nodeNumber = i;
			}
			else
			{
         		nodalUnit[i].get_leaf()->set_N_content(leaf_N_content); //SK 8/22/10: set leaf N content before each nodal unit is updated
                nodalUnit[i].update(develop, PredawnLWP); //Pass the predawn leaf water potential into a nodel
				//to enable the model to simulate leaf expansion with the
				//effect of predawn leaf water potential YY
			}
			if (nodalUnit[i].get_leaf()->isGrowing()) 
			{
				TotalGrowingLeaves++;
				
			}
			if (nodalUnit[i].get_leaf()->isDropped()) 
			{
				TotalDroppedLeaves++;
			}

		}
 //SK 8/22/10: Why is this infor being stored in node0? this is confusing because it is supposed to be the coleoptile

		nodalUnit[0].get_leaf()->set_TotalGrowingLeaves(TotalGrowingLeaves); //TODO make this a class variable and an CPlant:update method
		nodalUnit[0].get_leaf()->set_TotalDroppedLeaves(TotalDroppedLeaves);

		if (TotalDroppedLeaves >= develop->get_totalLeaves()) 
		{
			develop->death.done = true;
			develop->death.daytime = weather.daytime;
		}


//SK: Below {} lumps N related codes
	{
		double shootmass = this->get_shootMass();

		if (shootmass <=(100/initInfo.plantDensity)) //100 g m-2
		{
			double temp = 0.063*this->get_shootMass(); //when shoot biomass is lower than 100 g/m2, 
			                                                //the maximum [N] allowed is 6.3%; 
			                                               
			                                                //shoot biomass and Nitrogen are in g
			if (TotalNitrogen>temp)  //need to adjust demand or else there will be mass balance problems
			{
				this->set_N(temp);
			}
		}

		
        calcPerLeafRelativeAreaIncrease(); // calculate relative area increases for leaves 

    	calcPotentialLeafArea();
		calcGreenLeafArea();
		calcActualGreenLeafArea();
		leaf_N_content = leaf_N/(this->greenLeafArea/10000); //Calculate leaf nitrogen content of per unit area; 
		calcLeafArea();
		//leaf_N_content = leaf_N/this->actualGreenLeafArea; //Calculate leaf nitrogen content of per unit area; 
		//defining leaf nitrogen content this way, we did not consider the difference in leaf nitrogen content
		//of sunlit and shaded leaf yet YY
        
//		calcSenescentLeafArea();
		droppedLfArea = (1-greenLeafArea/potentialLeafArea)*potentialLeafArea; //calculated dropped leaf area YY
		//SK 8/20/10: Changed it to get all non-green leaf area
        currentDroppedLfArea = droppedLfArea - previousDroppedlfArea; //leaf dropped at this time step 
	//	this->set_N((this->get_N()-(leaf_N/leafArea)*currentDroppedLfArea)); //calculated the total amount of nitrogen after the dropped leaves take some nitrogen out 
		                                                    //no nitrogen remobilization from old leaf to young leaf is considered for now YY

    //SK 8/22/10: N remobilization is implicitly done by adjusting greenleaf area after determining senesced leaf area. Currently it is assumed that all N is moved from the senesced to the active
		previousDroppedlfArea = droppedLfArea;
		//when less than 5% of total leaf area is green, physiological maturity is reached. SK
		//see http://www.agry.purdue.edu/ext/corn/news/timeless/TopLeafDeath.html
		if ((greenLeafArea <= 0.05*leafArea) && !develop->maturity.done) 
		{
			develop->maturity.done = true;
			develop->maturity.daytime = weather.daytime;
			cout << "* Physiological maturity " << develop->get_GDDsum() << " T growth: " << develop->get_Tgrow() 
				<< " green leaf %: " << greenLeafArea/leafArea*100 << endl;
		}

	}


		calcGasExchange(weather);
		calcMaintRespiration(weather);
		C_allocation(weather);
		//	C_allocation(weather, potentialCarbonDemand); //Calculating carbon allocation with carbon demand known YY
		if (abs(weather.time)<0.0001)
		{
			C_reserve += __max(0, C_pool);
			C_pool = 0.0; //reset shorterm C_pool to zero at midnight, needs to be more mechanistic
		}
		else
		{
            C_pool += assimilate*CH2O_MW/CO2_MW; // convert from grams CO2 to grams carbohydrate (per hour per plant)
		}
		setMass();
 	}
}

void CPlant::setMass()
//TODO: allocate biomass into individual organs, currently it is allocated as a bulk to leaf, stem, and so on
//so individual leaf doesn't have mass but the first or the last one has it all
{
	double m = 0;
	
/*
double agefn=1;
	if (potentialLeafArea>0)
		{
			agefn= (leafArea/potentialLeafArea); // as more leaves drop leaf biomass should go down 
		}
    //agefn=1;
	
	
	//double agefn = ((LeafArea-senescentLeafArea)/potentialLeafArea);  // dt 8/26/2011 trying this again - will relate it to 
		                                                             // dropped leaf area? Yang, 6/22/2003
	{
		double lf = 0.0;
		// apply leaf drop
        lf = this->get_nodalUnit()->get_leaf()->get_mass()/greenLeafArea;
		this->get_nodalUnit()->get_leaf()->setMass(lf*agefn); 
		droppedLfmass = lf*(1-agefn); //droppedLfmass is the mass of dropped leaves. In this model, in the Cleaf class
		                          //the potential area of each leaf (and the sum of potential area of each leaf) is determined, so is the acutal area of 
		                          //of each leaf. If a leaf is determined to be a "dropped" leaf (when the "terminated" booling variable is set to be "true",
		                          //then the area of this leaf is set to 0. Therefore the total acutal leaf area may be smaller than
		                          //potential leaf area because of the senescense.But the mass of the dropped leaf is not taken out of the total leaf mass
		                          //until here in this setMass() function. The variable agefn actually records the ratio between 
		                          //actual leaf area and potential leaf area. In this sense, the droppedLfmass should be (1-agefn)*lf  YY
	   
	}
*/
	// dt the addition of C_reserve here only serves to maintain a total for the mass. It could have just as easily been added to total mass.
	// C_reserve is added to stem here to represent soluble TNC, SK
    stemMass = this->get_nodalUnit()->get_stem()->get_mass() + C_reserve;
    earMass = this->get_ear()->get_mass();
   // need to iterate here to set leaf mass
	leafMass=calcTotalLeafMass();
	leafMass = this->get_leafMass();
	activeLeafMass = calcActiveLeafMass();
	droppedLeafmass=calcDroppedLeafMass();
	rootMass = this->get_roots()->get_mass();
	shootMass = seedMass + stemMass + leafMass + earMass;
	mass = shootMass + rootMass;
   
	return;
}



double CPlant::calcLeafArea()
{
	double area = 0.0; 
	double dL = 0.0;
	for (int i = 1; i <= develop->get_LvsInitiated(); i++)
	{
 		dL = nodalUnit[i].get_leaf()->get_area();
		area += dL ;
	}
	leafArea = area;
	return area;
}

double CPlant::calcGreenLeafArea()
{
	double area = 0.0;
	for (int i = 1; i <= develop->get_LvsInitiated(); i++)
	{
		area += nodalUnit[i].get_leaf()->get_greenArea();
	}
	greenLeafArea = area;
	return area;
}

double CPlant::calcActualGreenLeafArea()
{
	double area = 0.0;
	for (int i= 1; i<=develop->get_LvsInitiated(); i++)
	{
		area +=nodalUnit[i].get_leaf()->get_actualgreenArea();
	}
	actualGreenLeafArea=area;
	return area;
}

double CPlant::calcSenescentLeafArea()
{
	double area = 0.0;
	for (int i = 1; i <= develop->get_LvsInitiated(); i++)
	{
		area += nodalUnit[i].get_leaf()->get_senescentArea();
	}
	senescentLeafArea = area;
	return area;
}

double CPlant::calcPotentialLeafArea()
{
	double area = 0.0;
	for (int i = 1; i <= develop->get_LvsInitiated(); i++)
	{
		area += nodalUnit[i].get_leaf()->get_potentialArea();
	}
	potentialLeafArea = area;
	return area;
}

double CPlant::calcPotentialLeafAreaIncrease()
{
	double areaIncrease = 0.0;
	for (int i=1; i<=develop->get_LvsInitiated(); i++)
	{
		areaIncrease +=nodalUnit[i].get_leaf()->get_potentialAreaIncrease();
	}
	potentialLeafAreaIncrease = areaIncrease;
	return areaIncrease;
}

// calculate relative area increases for leaves now that they are updated

void  CPlant::calcPerLeafRelativeAreaIncrease()
	{
	  double PotentialLeafAreaIncrease=calcPotentialLeafAreaIncrease();
	  double RelativeLeafAreaIncrease=0;

      for (int i = 1; i <= develop->get_LvsInitiated() ; i++)
		{
			RelativeLeafAreaIncrease=0.0; //Relative area increase is used for carbon partitioning later.
			if (PotentialLeafAreaIncrease>0)
			    {
				    RelativeLeafAreaIncrease=nodalUnit[i].get_leaf()->get_potentialAreaIncrease()/PotentialLeafAreaIncrease;
					nodalUnit[i].get_leaf()->set_RelativeAreaIncrease(RelativeLeafAreaIncrease);
    			}
    	}
	}

double CPlant::calcActiveLeafMass() //this is the total mass of active leaves that are not entirely dead (e.g., dropped).
//It would be slightly greather than the green leaf mass because some senesced leaf area is included until they are complely aged (dead), SK
{
	double Mass = 0.0;
	for (int i=1; i<=develop->get_LvsInitiated(); i++)
	{
		if (!nodalUnit[i].get_leaf()->isDropped())
			{
				Mass +=nodalUnit[i].get_leaf()->get_mass();
			}
	}
	activeLeafMass=Mass;
	return Mass;
}
double CPlant::calcDroppedLeafMass() //It has been pointed that corn leaves don't really drop and is still likely to be attached to the node. 
// this is true so a better term would be "dead" leaves (entirely senesced) but I am fine with continuing to use dropped leaf as long as we are clear about the meaning, SK
{
	double Mass = 0.0;
	for (int i=1; i<=develop->get_LvsInitiated(); i++)
	{
		if (nodalUnit[i].get_leaf()->isDropped())
			{
				Mass+=nodalUnit[i].get_leaf()->get_mass();
			}
	}
	droppedLeafmass=Mass;
	return Mass;
}
double CPlant::calcTotalLeafMass()
{
	double TotalMass = 0.0;
	double Mass = 0.0;
	double Area = 0.0;
	for (int i=1; i<=develop->get_LvsInitiated(); i++)
	{
		Mass = nodalUnit[i].get_leaf()->get_mass();
		Area = nodalUnit[i].get_leaf()->get_area();
		nodalUnit[i].get_leaf()->set_SLA(Area/Mass); // Set SLA from current leaf area and mass, SK
		TotalMass += Mass;
	}
	leafMass=TotalMass; //this should equal to activeLeafMass + droppedLeafMass;
	return TotalMass;
}

double CPlant::calcPotentialCarbondemand()
{ //this will only be used for total leaf area adjustment. If the individual leaf thing works 
  // out this will be deleted. 
	double SLA=200; //Just a mocking value for now. Need to find a more mechanistic way to simulate change in SLA YY
	// SK 8/20/10: changed it to 200 cm2/g based on data from Kim et al. (2007) EEB
	 

	double LeafMassDemand = potentialLeafAreaIncrease/SLA; //units are biomass not carbon
	//PotentialCarbonDemand=carbondemand; //for now only carbon demand for leaf is calculated. 
	return LeafMassDemand;
}
void CPlant::calcGasExchange(const TWeather & weather)
{
	const double tau = 0.50; // atmospheric transmittance, to be implemented as a variable => done
	const double LAF = 1.37; // leaf angle factor for corn leaves, Campbell and Norman (1998)
	const double leafwidth = 5.0; //to be calculated when implemented for individal leaves
	const double atmPressure= 100.0; //kPa, to be predicted using altitude
	double activeLeafRatio = greenLeafArea/leafArea;
	double LAI = greenLeafArea*initInfo.plantDensity/(100.0*100.0);

	//CGas_exchange * sunlit = new CGas_exchange("Sunlit");
	//CGas_exchange * shaded = new CGas_exchange("Shaded");

	  CGas_exchange * sunlit = new CGas_exchange("Sunlit", this->leaf_N_content);
	  CGas_exchange * shaded = new CGas_exchange("Shaded", this->leaf_N_content);

//TODO: lightenv.dll needs to be translated to C++. It slows down the execution, 3/16/05, SK
	radTrans2(weather.jday, weather.time, initInfo.latitude, initInfo.longitude, weather.solRad, weather.PFD, LAI, LAF);
    double temp5, temp6, temp7,temp8, temp9;
	temp5=sunlitLAI();
	temp6=shadedLAI();
	temp7=getNIRtot();
    temp8=sunlitPFD();
	temp9=shadedPFD();
	 	
	  //Calculating transpiration and photosynthesis without stomatal control Y
	   // sunlit->SetVal_NC(sunlitPFD(), weather.airT, weather.CO2, weather.RH, 
	       //           weather.wind, atmPressure, leafwidth, weather.LeafWP, weather.ET_supply); 
	    // shaded->SetVal_NC(shadedPFD(), weather.airT, weather.CO2, weather.RH, 
	          //     weather.wind, atmPressure, leafwidth, weather.LeafWP, weather.ET_supply);

	//Calculating transpiration and photosynthesis with stomatal controlled by leaf water potential LeafWP Y
	    sunlit->SetVal_psil(sunlitPFD(), weather.airT, weather.CO2, weather.RH, 
	                weather.wind, atmPressure, leafwidth, weather.LeafWP, weather.ET_supply*initInfo.plantDensity/3600/18.01/LAI); 
	    shaded->SetVal_psil(shadedPFD(), weather.airT, weather.CO2, weather.RH, 
	          weather.wind, atmPressure, leafwidth, weather.LeafWP, weather.ET_supply*initInfo.plantDensity/3600/18.01/LAI);

	
	photosynthesis_gross = (sunlit->A_gross*sunlitLAI() + shaded->A_gross*shadedLAI());//plantsPerMeterSquare units are umol CO2 m-2 ground s-1 ;
	photosynthesis_net = (sunlit->A_net*sunlitLAI() + shaded->A_net*shadedLAI());//
	SunlitLAI = sunlitLAI();
	ShadedLAI = shadedLAI();
	transpirationOld=transpiration;  // when outputting the previous step transpiration is compared to the current step's water uptake
	transpiration=0;
	if (sunlitLAI()>0) transpiration=sunlit->ET*sunlitLAI();
	if (shadedLAI()>0) transpiration+=shaded->ET*shadedLAI();
	transpiration=transpiration/(initInfo.plantDensity)*3600.0*18.01;//plantsPerMeterSquare units are grams per plant per hour ;
	// Units of Transpiration from sunlit->ET are mol m-2 (leaf area) s-1
	// Calculation of transpiration from ET involves the conversion to gr per plant per hour 
	temperature = (sunlit->Tleaf*sunlitLAI() + shaded->Tleaf*shadedLAI())/LAI;
	//psi_l = (sunlit->get_psi()*sunlitLAI() + shaded->get_psi()*shadedLAI())/LAI;
	this->VPD = sunlit->get_VPD();
	// photosynthesis_gross is umol CO2 m-2 leaf s-1
	// in the following we convert to g C plant-1 per hour
	assimilate = (photosynthesis_gross*CO2_MW/1.0e6)*(60.0*initInfo.timeStep)/initInfo.plantDensity; // grams CO2 per plant per hour
	photosynthesis_gross=photosynthesis_gross*CH2O_MW/1.0e6*(60.0*initInfo.timeStep)/initInfo.plantDensity; //grams carbo per plant per hour
	photosynthesis_net=  photosynthesis_net*CH2O_MW/1.0e6*(60.0*initInfo.timeStep)/initInfo.plantDensity; //grams carbo per plant per hour
	double temp1=sunlitLAI();
	double temp2=shadedLAI();
	double temp3=sunlit->get_gs();
	double temp4=shaded->get_gs();
	if (sunlitLAI()!=0 && shadedLAI() !=0 && LAI !=0)
	{
		this->conductance=(sunlit->get_gs()*sunlitLAI()+shaded->get_gs()*shadedLAI())/LAI; //average stomatal conductance Yang
		if (this->conductance<0)
		{
			conductance=0;
		}
	
	}
	else 
	{
		this->conductance =0;
	}
	 delete sunlit;
	 delete shaded;
}

void CPlant::C_allocation(const TWeather & w)
// this needs to be f of temperature, source/sink relations, nitrogen, and probably water
// a valve function is necessary because assimilates from CPool cannot be dumped instantanesly to parts
// this may be used for implementing feedback inhibition due to high sugar content in the leaves
// The following is based on Grant (1989) AJ 81:563-571
{
   double b1=2.325152587; // Normalized (0 to 1) temperature response fn parameters, Pasian and Lieth (1990)
                        // Lieth and Pasian Scientifica Hortuculturae 46:109-128 1991
                        // parameters were fit using rose data -
   double b2=0.185418876; // I'm using this because it can have broad optimal region unlike beta fn or Arrhenius eqn
   double b3=0.203535650;
   const double Td = 48.6; //High temperature compensation point
  // double shootPart = 0.0;
  // double rootPart = 0.0;
   double shootPart_real = 0.0;
   double rootPart_real = 0.0;
   //double leafPart = 0.0; made class variable for now to transfer, may have to do this for all
   double sheathPart = 0.0;
   double stalkPart = 0.0;
   double reservePart = 0.0;
   double huskPart = 0.0;
   double cobPart = 0.0;
   double grainPart = 0.0;
   double flag = 0.0;

   double g1=1+exp(b1-b2*w.airT);
   double g2=0.0;
   if (w.airT<Td) g2=1-exp(-b3*(Td-w.airT));

   double tmprEffect = g2/g1;   
   double grofac = 1/(5*60/initInfo.timeStep); // translocation limitation and lag, assume it takes 1 hours to complete, 
                                              // 0.2=5 hrs
   // this is where source/sink (supply/demand) valve can come in to play
   // 0.2 is value for hourly interval, Grant (1989)
	double scale = 0.0; // see Grant (1989), #of phy elapsed since TI/# of phy between TI and silking
	scale = develop->get_phyllochronsFromTI()/(develop->get_youngestLeaf() - develop->get_LvsAtTI());
	//		if (w.time == 0.0) std::cout << scale << endl;

    C_supply = 0.0;  // daily mobilization of carbon
	const double C_min = 1.0;


    if (C_pool > C_demand) //C_demand does not enter into equations until grain fill
	{
        C_supply = __max(C_pool*tmprEffect*grofac, 0); //CADD from Grant
        C_pool -= C_supply; // C_Pool is what is left over from the carbon available for growth
	}
	else
    if (abs(C_pool)<0.0001)       //C_pool is zero 
	{
		if (C_reserve >0)
		{
			C_supply = __max(C_reserve*tmprEffect*grofac, 0); //All the reserve is not available
			C_reserve-=C_supply; //reduce reserve pool for used carbon
		}
	}
	else
	if ((C_reserve > C_demand) && (C_demand>0))
	{
		// conversion and translocation from long term reserve should be less efficient, apply nother glucose conversion factor
		// translocation from the soluble sugar reserve
		if (C_pool <= 0.0)
		{
			C_reserve += C_pool; // C_pool negative here, add instead of subtract
			C_pool = 0;
		}
		// deplete C_pool first* tmprEffect
		C_supply = __max(C_demand* tmprEffect*grofac, 0);
		C_reserve -= C_supply; //reserve C is used to increase grain mass
		C_reserve += C_pool; // send remaining C (not enough to meet the demand) in shorterm pool to reserve
		C_pool = 0; // empty the C_pool
	}
	else
	if (C_pool > (maintRespiration))
	{
        C_supply = maintRespiration;
		C_pool -= __max(C_supply,0);

	}
	else
	if (C_reserve > (maintRespiration))
	{
        C_supply = maintRespiration;
		C_reserve -= __max(C_supply,0);
		C_reserve += C_pool; // send remaining C (not enough to meet the demand) in shorterm pool to reserve
		C_pool = 0; // empty the C_pool
		//In this way, the c_reserve is used to satisfy the maintainance respiration demand. Then C_pool is 
		//used to replenish C_reserve and then C_pool is set to 0. In other words, C_pool is depleted first
		//and then C_reserve is used to supplement whatever demand that's left
	}
	else
	{
		C_reserve += C_pool;
		C_pool = 0.0;
		C_supply = __min(C_reserve, maintRespiration);
	}
//	std::cout  << setw(15)<< setprecision(8) <<w.daytime  
//		<< " " 
//		<< setw(10)<< std::setprecision(5)<< C_pool << " " 
//		<< setw(10)<< std::setprecision(5)<< C_demand <<" " 
//		<< setw(10)<< std::setprecision(5)<< C_supply << " " 
//		<< setw(10)<< std::setprecision(5)<<C_reserve <<endl ;
	double Fraction = __min(0.925, 0.50 + 0.50*scale); // eq 3 in Grant
//	const double convFactor = 1/1.43; // equivalent Yg, Goudriaan and van Laar (1994)
	double Yg = 0.750; // synthesis efficiency, ranges between 0.7 to 0.76 for corn, see Loomis and Amthor (1999), Grant (1989), McCree (1988)
//  double Yg = 0.74;
	
	// this is the same as (PhyllochronsSinceTI - lvsAtTI/(totalLeaves - lvsAtTI)
	shootPart = __max(0,Yg*(Fraction*(C_supply-maintRespiration))); // gCH2O partitioned to shoot
	rootPart = __max(0,Yg*((1-Fraction)*(C_supply-maintRespiration))); // gCH2O partitioned to roots

    if (!develop->Germinated())
   {
	   return;
   }
    else if (!develop->TasselInitiated())
   {
	
   // shootPart was reduced to 0.37; rootPart was 0.43 in sourcesafe file yy
// SK, commenting it out. Yg needs to be multiplied here because it represents growth respiration.
//		   shootPart = __max(0, 0.67*(C_supply-maintRespiration)); //these are the amount of carbons allocated with no drought stress
//         rootPart = __max(0, 0.33*(C_supply-maintRespiration));  //Yang, 6/22/2003


	   if (w.pcrs>rootPart_old) // if in time step t-1, the value of pcrs is higher than that of pcrl
	   { 
// give a half of carbon from shoot needed to meet root demand? SK
	      shootPart_real = __max(0, shootPart-(w.pcrs-rootPart_old)); //than take the difference between the two out of the carbon allocation to shoot at time step t
		  rootPart_real = rootPart+ (w.pcrs-rootPart_old);  //and put that amount of carbon into carbon allocation to root at time step t.          
	   }
	   else 
	   {
		   C_pool_root +=__max(0,rootPart_old-w.pcrs);
		   shootPart_real = shootPart;
		   rootPart_real = rootPart-__max(0,rootPart_old-w.pcrs); //subtract out carbon sent to the root pool
		   
	   }

	   rootPart_old = rootPart;
 
	   leafPart = shootPart_real*0.725;
	   sheathPart = shootPart_real*0.275;
	   stalkPart = 0.0;
	   reservePart = 0.0;
	   huskPart = 0.0;
	   cobPart = 0.0;
	   grainPart = 0.0;
   }
	
   else if (!develop->GrainFillBegan())
   {
		//C partitioning
	shootPart = __max(0,Yg*((Fraction)*(C_supply-maintRespiration))); // gCH2O partitioned to shoot
	//Everything is partitioned to the shoot after grain filling begins. I removed Fraction from the Equation above
	rootPart = __max(0,Yg*((1-Fraction)*(C_supply-maintRespiration))); // gCH2O partitioned to roots
//	 rootPart=0.0;

		if (w.pcrs>rootPart_old)
	   { 
		  shootPart_real = __max(0, shootPart-(w.pcrs-rootPart_old));
		  rootPart_real = rootPart+(w.pcrs-rootPart_old);
	   }
	   else 
	   {
		   C_pool_root +=__max(0,rootPart_old-w.pcrs);
		   shootPart_real = shootPart;
		   rootPart_real = rootPart-__max(0,rootPart_old-w.pcrs); //subtract out carbon sent to the root pool		   
	   }

	    rootPart_old = rootPart;
		leafPart = shootPart_real*__max(0.725 - 0.775*scale,0);
		sheathPart = shootPart_real*__max(0.275- 0.225*scale,0);
// allocate shootPart into components
		if (scale <=0.85)
		{
		  stalkPart = shootPart_real*1.10*scale;
		}
		else
		{
		 reservePart = shootPart_real*__max(2.33-0.6*exp(scale), 0);
		}

		if (scale <= 1.0)
		{
		  huskPart = shootPart_real*__max(exp(-7.75+6.60*scale),0);
		}
		else 
		{
		 huskPart = shootPart_real*__max(1.0 - 0.675*scale,0);
		}

		if (scale <= 1.125) 
		{
			cobPart = shootPart_real*exp(-8.4+7.0*scale);
		}
		else
		{
            cobPart = shootPart_real*0.625;
		}
		//give reserve part what is left over, right now it is too high
		if (reservePart>0)
		{
			double sum=cobPart+huskPart+leafPart+stalkPart+sheathPart;
			reservePart=__max(0,shootPart-sum);
		}
   }
   else if (!develop->Dead())//no acutally kernel No. is calculated here ? Yang, 6/22/2003
   {
       // here only grain and root dry matter increases root should be zero but it is small now. 
	   const int maxKernelNo = 800; // assumed maximum kerner number per ear
	   double maxKernelFillRate = 0.012*(initInfo.timeStep/(24*60)); // 
	                                                                  //max kernel filling rate = 0.012g Kernel-1 day-1, Grant (1989)
	   C_demand = maxKernelNo*maxKernelFillRate*tmprEffect*C_content; //dt added c_content
//	   shootPart = __max(0,Yg*(Fraction*(C_supply-maintRespiration))); // gCH2O partitioned to shoot
//	   rootPart = __max(0,Yg*((1-Fraction)*(C_supply-maintRespiration))); // gCH2O partitioned to roots
	   shootPart = __max(0,Yg*(C_supply-maintRespiration)); // gCH2O partitioned to shoot
	   rootPart=0.0; // no more partitioning to root during grain fill

      
	   if (w.pcrs>rootPart_old)
	   { 
	      shootPart_real = __max(0, shootPart-(w.pcrs-rootPart_old));
		  rootPart_real = rootPart+ (w.pcrs-rootPart_old);
	   }
	   else 
	   {
		   C_pool_root +=__max(0,rootPart_old-w.pcrs);
		   shootPart_real = shootPart;
		   rootPart_real = rootPart-__max(0,rootPart_old-w.pcrs); //subtract out carbon sent to the root pool
		   flag = 0;
	   }

	   rootPart_old = rootPart;
	   
	   grainPart = shootPart_real*1.0;
   }

   else
   {

   }
   double stemPart = sheathPart + stalkPart; //TODO: sheath and stalk haven't been separated in this model
                                             //reservePart needs to be added later
   double earPart = grainPart + cobPart + huskPart;
   double sum= stemPart+earPart + leafPart;
// here we can allocate leaf part among leaves
   CLeaf* leaf;
 

// Partition carbon to leaves relative to growth rate
// now need to find increment of Carbo to add to each leaf
// update leaf's mass
// first find the number of growing leaves. This is saved as a variable in the [0] nodal unit
   double LeafPartSum=leafPart;
   double PCarboDemandPerLeaf;
   
   for (int i = develop->get_LvsInitiated(); i >=1  ; i--)
   {
	   //need to find demand first
	   leaf= nodalUnit[i].get_leaf();  
	 //Update SLA based on current shoot growth rate. Do this for every leaf
	   // no sla adjustment until after emergence
//	   double SLA_est;
//	   if (shootPart_real > 0 && develop->Emered()) 
//	   {
//		   SLA_est=1/(25.0+150.0*shootPart_real)*10000; // see Grant, 1989, eq 9. 10000 converts from
//		                                                // m2 to cm2
//	       if (leaf->isGrowing()) leaf->set_SLA(min(400.0,SLA_est));
//	   }
	    double totLA;

	    if (nodalUnit[i].isInitiated()&& LeafPartSum >= 0.0)
		{
			
			

			if (!leaf->isDead())
			{
			
                totLA = calcPotentialLeafArea();
				PCarboDemandPerLeaf=__max(0.0, leaf->get_potentialArea()/totLA*leafPart); //Adjusting C allocation based on leaf size if not aging. doing it based on current growth rate is more mechanistic but seems to have issue now. To revisit. SK
//			   PCarboDemandPerLeaf=__max(0.0, leaf->get_RelativeAreaIncrease()*leafPart); // carbon allocated according to growth rate

				if (LeafPartSum >= PCarboDemandPerLeaf)
			   {
					leaf->import_CH2O(PCarboDemandPerLeaf);
					LeafPartSum-=PCarboDemandPerLeaf;
			   }
			   else
			   {
					leaf->import_CH2O(__min(LeafPartSum, PCarboDemandPerLeaf));
					LeafPartSum=0.0;
			   }

			}
		}
   }
 
   //cout <<"leafPartSum: " <<LeafPartSum <<"LeafPart: " <<leafPart<<endl;

   this->get_roots()->set_ActualCarboIncrement(rootPart_real);
   this->get_nodalUnit()->get_stem()->import_CH2O(stemPart);
//   this->get_nodalUnit()->get_leaf()->import_CH2O(leafPart);
   this->C_reserve += reservePart; // everything is carbohydrate now
   // before emergence root weight has been initialized. Just dump this carbon for now.
   if (develop->Emerged()) this->get_roots()->import_CH2O(rootPart_real);
   this->get_ear()->import_CH2O(earPart);

   double partSum = stemPart + earPart + leafPart; // checking the balance if sums up to shootPart
//   if (w.time == 0.5) std::cout << "adding CH2O: " << grainPart << " to grain " << endl;
//   if (w.time == 0.5) std::cout << "Sum of part coeff is " << partSum << " and shoot CH2O is " << shootPart <<  endl;
    //cout <<C_pool << " " << C_ReserveLeaf <<endl;
}


void CPlant::calcMaintRespiration(const TWeather & w)
// based on McCree's paradigm, See McCree(1988), Amthor (2000), Goudriaan and van Laar (1994)
// units very important here, be explicit whether dealing with gC, gCH2O, or gCO2
  {
	const double Q10 = 2.0; // typical Q10 value for respiration, Loomis and Amthor (1999) Crop Sci 39:1584-1596
	//const double Q10 = 2.1; // where does this value come from?
	double dt = initInfo.timeStep/(24*60);
//	const double maintCoeff = 0.015; // gCH2O g-1DM day-1 at 20C for young plants, Goudriaan and van Laar (1994) Wageningen textbook p 54, 60-61
	const double maintCoeff = 0.018;
	double agefn = (greenLeafArea+1.0)/(leafArea+1.0); // as more leaves senesce maint cost should go down, added 1 to both denom and numer to avoid division by zero. 
	//no maint cost for dead materials but needs to be more mechanistic, SK
	agefn=1.0;
	double q10fn = pow(Q10,(w.airT - 20.0)/10); // should be soil temperature
	maintRespiration = q10fn*maintCoeff*mass*dt;// gCH2O dt-1, agefn effect removed. 11/17/14. SK.
}
void CPlant::calcRed_FRedRatio(const TWeather &weather)
// this function calculates an estimate of the Red to Far red light ratio from sunlit and shaded ET. This 
// ration is used to estimate the effects of plant density on leaf expansion and LAI. 
// A daily mean ratio is calculated. We use a 3 parameter sigmoid function to model the effect
{
	//double Xo=0.43, B=0.05, A=1.2;
	double Xo=0.6, B=0.13, A=2.0;
	double dt=initInfo.timeStep/(24.0*60);
	double C2_effectTemp;
      //First set counter to 0 if it is the beginning of the day. 
	if (abs(weather.time)<0.0001)
		{// have to rename C2_effect to Light_effect
			//Journal of Experimental Botany, Vol. 65, No. 2, pp. 641–653, 2014
			C2_effectTemp=exp(-(SunlitRatio-Xo)/B);
			C2_effect=__min(1.0,A/(1.0+C2_effectTemp));
			develop->set_shadeEffect(C2_effect);
			SunlitRatio=0.0;
			
		}
	else
	{
		if (SunlitLAI/(SunlitLAI+ShadedLAI)>0.05) // calculate from emergence
		{

			SunlitRatio+=SunlitLAI/(SunlitLAI+ShadedLAI)*dt;
	
		}
	
		else SunlitRatio+=1.0*dt;

	}
     
	
}


void CPlant::writeNote(const TWeather & w)
{
	ostringstream oStr;
	string s = "";
	if (FLOAT_EQ(develop->germination.daytime,w.daytime)){s = "Germinated";}
	if (FLOAT_EQ(develop->emergence.daytime,w.daytime)){s = "Emergence";}
	if (FLOAT_EQ(develop->tasselInitiation.daytime,w.daytime)){s = "Tassel Initiation";}
	if (FLOAT_EQ(develop->anthesis.daytime, w.daytime)){s = "Anthesis";}
	if (FLOAT_EQ(develop->silking.daytime, w.daytime)){s = "Silking";}
	if (FLOAT_EQ(develop->beginGrainFill.daytime, w.daytime)){s = "Begin grain filling";}
	if (FLOAT_EQ(develop->maturity.daytime,w.daytime)){s = "Begin grain filling";}

	if (s != "")
	{
        oStr << s << " "<< w.jday;
	}
	//note.swap(oStr.str());
	note = oStr.str();
}