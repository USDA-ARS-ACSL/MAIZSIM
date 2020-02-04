#pragma once
#ifndef _GAS_EX_SPECIES_PARAM_H_
#define  _GAS_EX_SPECIES_PARAM_H_


struct TGasExSpeciesParam
{
	// TODO these need to be documented, 1 variable per line
public:
	TGasExSpeciesParam()
	{

		Vpm25,
		Vcm25,
		Jm25,
		Rd25,
		EaVp,
		EaVc,
		Eaj,
		Hj,
		Sj,
		Ear,
		g0,
		g1,
		SC_param, 
		BLC_param,
		gamma1,
		Gamma_gsw,
		sf,
		phyf, 
		stomaRatio, 
		widthPara,
		internalCO2Ratio, //holds factor to multiply CO2 by to get first estimate of Ci
		f,  //spectral correction
		scatt,
		Kc25,
		Ko25,
			Kp25,
			gbs,
			gi,
			phyf,
			LfWidth,
			LfAngFact;

	}
public:
	double Vpm25, Vcm25, Jm25, Rd25, EaVp, EaVc, Eaj, Sj, Hj, Ear, g0, g1;
	double SC_param, BLC_param;
	double   gamma1, Gamma_gsw,sf, phyf, stomaRatio, widthPara;
	double   internalCO2Ratio;
	double f,  
		scatt;//spectral correction
	int
		Kc25,
		Ko25,
		Kp25;
	double
		gbs,
		gi,
		LfWidth,
		LfAngFact;
	

};
#endif