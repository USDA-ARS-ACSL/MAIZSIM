
#include "gas_ex_species_param.h"

/*! \class CGasExchange 
* \brief Class for gas exchange calculations\n
* This class simulates gas exchange in plant leaves for C3 and C4 plants. \n
* \par Usage
	- Use <b>SetParams</b> to initialize the GasExchange object with parameters for a specific variety of plant
	- Use <b>SetVal</b> to pass environmental variables for a simulation and return a structure with output. \n
* See \ref Interface for details on input and output	
    

	*/

 
class CGasExchange  
	{


	public:
		CGasExchange(std::string sType_in, double n_content, const TGasExSpeciesParam& photoparam);

		~CGasExchange(void);



		/*!    \struct tparms gas_exchange.h
		       \brief Structure to hold plant parameters for the photosynthesis model. \n
			   Some parameters are specific for C3 or C4 type Plants \n
		*/   
		
		/*!
		//                parameter        description
		@param 		Vcm25	4		Photosynthetic Rubisco Capacity at 25C (umol m-2 s-1)
		@param 		Jm25	5		Potential Rate of electron transport at 25C  (umol m-2 s-1)
		@param 		Vpm25	6		C4 Carboxylation rate at 25C (C4, umol m-2 s-1)
		@param 		TPU25	7		Rate if Triose Phosphate Utilization at 25C (C3, umol m-2 s-1)
		@param 		RD25	8		Mitochondrial respiration in the light at 25C (umol m-2 s-1)
		@param 		Theta   9		Initial slope of CO2 response (umol m2 s-1) - de Pury (1997)
		@param 		EaVc   10		Activation energy for Arrhenius function used to calculate temperature dependence for Vcmax (kJ mol-1)	
		@param 		Eaj    11		Activation energy for Arrhenius function used to calculate temperature dependence for J (kJ mol-1)
		@param 		Hj     12		Curvature parameter of the temperature dpendence of Jmax (kJ mol-1)
		@param		Sj	   13		Electron transport temperature response parameter for Jmax (J mole-1 K-1)
		@param      Hv     14       Curvature parameter of the temperature dependence of Vcmax (J mole-1)
		@param 	    EaVp   15		Activation energy for Arrhenius function used to calculate temperature dependence for Vpmax (kJ mol-1)
		@param 	    Sv	   16		Electron transport temperature response parameter for Vcmax (J mole-1 K-1)
		@param 		EAP	   17		Activation energy for Arrhenius function used to calculate temperature dependence for TPU (kJ mol-1)
		@param 		EAR	   18	 	Activation energy for Arrhenius function used to calculate temperature dependence for respiration (kJ mol-1) 
		@param 		g0	   19		Minimum stomatal conductance to water vapor at the light compensation point in the BWB model (mol m-2 s-1)	
		@param 	    g1	   20   	Empirical coefficient for the sensitivity of StomatalConductance to A, Cs and hs in BWB model (no units)
		@param 		StomRatio	21	Stomatal Ratio (fraction)
		@param  	LfWidth     22	Leaf Width (m)
		@param 		LfAngFact	23	Leaf Angle Factor		
		@param 		Remark		24	Text 
		*/	
	
		//These functions return results of calculations 

		double get_ANet() {return AssimilationNet;}      //!< return net photosynthesis (umol CO2 m-2 s-1) 
		double get_AGross() {return AssimilationGross;}  //!< return gross photosynthesis  (umol CO2 m-2 s-1)
		double get_Transpiration()     {return Transpiration;}   //!< return transpiration rate (umol H2O m-2 s-1)
		double get_LeafTemperature() {return Tleaf;} //!< return leaf temperature (C)
		double get_Ci()    {return Ci;}        //!< return internal CO2 concentration (umol mol-1)
		double get_StomatalConductance()    {return StomatalConductance;}    //!< return stomatal conductance to water vapor (mol m-2 s-1)
		double get_BoundaryLayerConductance()    {return BoundaryLayerConductance;}    //!< return boundary layer conductance (mol m-2 s-1)
		double get_Respiration() {return DarkRespiration;}   //!< return respiration rate (umol CO2 m-2 s-1)
		double get_VPD(){ return VPD;}         //!< return vapor pressure deficit (kpa)
		double get_leafpEffect() { return this->leafpEffect; }
		double set_leafpEffect(double pressure);

		void SetVal(double PhotoFluxDensity, double Tair, double CO2, double RH, 
			double wind, double Press, 
			double width, double leafp, double ET_supply);
		//!< sets input values for calculations for a particular set of environmental variables



	private:
		void GasEx(double leafp, double ET_supply);  //!<Main module to calculate gas exchange rates
		void Photosynthesis(double Ci);  //!<Function to calculate C4 photosynthesis
		void EnergyBalance();     //!<calculates leaf temperature and transpiration
		double SearchCi(double CO2i);          //!< called iterively to find optimal internal CO2 concentration returns optimal internal CO2 concentration (CO2i) 
		double EvalCi(double Ci);   //!<Calls photosynthesis modules to evaluate Ci dependent calculations returns difference between old Ci and new Ci
		double CalcStomatalConductance(double pressure);             //!< Stomatal conductance (mol m-2 s-1)
		double CalcTurbulentVaporConductance();             //!<  Conductance for turbulant vapor transfer in air - forced convection (mol m-2 s-1)
		double Es(double Temperature);      //!< Saturated vapor pressure at given temperature. kPa
		double Slope(double Temperature);  //!<  Slope of the vapor pressure curve
		double QuadSolnUpper (double a, double b, double c ); //!< Upper part of quadratic equation solution
		double QuadSolnLower (double a, double b, double c ); //!< Lower part of quadratic equation solution
		double minh(double fn1,double fn2,double theta2); //!< Hyperbolic minimum
		double get_CiCaRatio()  {return Ci_Ca;} //!< Ratio between internal and atmospheric CO2
		double leafp, leafpEffect;  //leafp is leaf water pressure, and leafpEffect is the effect of leafp on stomatal conductance. Yang 8/20/06
		double lfNContent; //lfNContent is leaf nitrogen content in unit of g m-2(leaf) YY

		// variables passed as arguments to the constructor
		std::string sType;
		TGasExSpeciesParam sParms;
		
		double  PhotoFluxDensity,  //!< Photosynthetic Flux Density (umol photons m-2 s-1 
			R_abs, //!< Absorbed incident radiation (watts m-2)        
			Tair,  //!< Air temperature at 2m, (C) 
			CO2,   //!< CO2 concentration (umol mol-1 air) 
			RH,   //!<  Relative Humidity (%, i.e., 80) 
			wind, //!<  Windspeed at 2 meters (km s-1) 
			width, //!< Leaf width (m) 
			Press,  //!<  Air pressure (kPa) 
			Theta; //Initial slope of CO2 response(umol m2 s - 1) - de Pury(1997)

		// These variables hold the output of calculations 
		double AssimilationNet,    //!< Net photosynthesis (umol CO2 m-2 s-1)
			AssimilationGross, //!< Gross photosynthesis (umol CO2 m-2 s-1) (Adjusted for respiration)
			Transpiration,     //!< Transpiration mol H2O m-2 s-1 
			Tleaf,  //!< Leaf temperature C 
			Ci,     //!< Internal CO2 concentration umol mol-1 
			StomatalConductance,     //!< Stomatal conductance umol m-2 s-1 
			BoundaryLayerConductance,    //!< Boundary layer conductance umol m-2 s-1 
			DarkRespiration,    //!< Plant respiration    umol m-2 s-1 
			VPD,    //!< Vapor Pressure Density, kPa */
			Ci_Ca;  //!< Ratio of internal to external CO2, unitless
	 
		double errTolerance; /*!< Error tolerance for iterations */
		double eqlTolerance; /*!< Equality tolerance */

		int iter_total;      //!< Total number of iterations */
		///  @param PlantType string that holds the type of plant, C3 or C4

		bool ConstantLeafTemperature; //!< if true, uses constant temperature - if true, does not solve for leaf temperature */
		int iter1, iter2;         //!< Iteration counters
		int  iter_Ci;   /*!< Iteration value for Ci umol mol-1, internal CO2 concentration */ 
		bool isCiConverged; /*!< True if Ci (internal CO2 concentration) iterations have converged */
		inline double Square(double a) { return a * a; } /*!< Squares number */ 
		inline double Min(double a, double b, double c) {return (__min(__min(a,b),c));} /*!< Finds minimum of three numbers */
		
	};


 
