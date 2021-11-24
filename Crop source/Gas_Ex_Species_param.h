#pragma once
#ifndef _GAS_EX_SPECIES_PARAM_H_
#define  _GAS_EX_SPECIES_PARAM_H_


struct TGasExSpeciesParam
{
	// TODO these need to be documented, 1 variable per line
public:
	TGasExSpeciesParam()
	{

		Vpm25=0,
		Vcm25 = 0,
		Jm25 = 0,
		Rd25 = 0,
		EaVp = 0,
		EaVc = 0,
		Eaj = 0,
		Hj = 0,
		Sj = 0,
		Ear = 0,
		g0 = 0,
		g1 = 0,
		SC_param = 0,
		BLC_param = 0,
		gamma1 = 0,
		Gamma_gsw = 0,
		sf = 0,
		phyf = 0,
		stomaRatio = 0,
		widthPara = 0,
		internalCO2Ratio = 0, //holds factor to multiply CO2 by to get first estimate of Ci
		f = 0,  //spectral correction
		scatt = 0,
		Kc25 = 0,
		Ko25 = 0,
			Kp25 = 0,
			gbs = 0,
			gi = 0,
			phyf = 0,
			LfWidth = 0,
			LfAngFact = 0;

	}
public:
	double Vpm25, Vcm25, Jm25, Rd25, EaVp, EaVc, Eaj, Sj, Hj, Ear, g0, g1;
	double SC_param, BLC_param;
	double   gamma1, Gamma_gsw,sf, phyf, stomaRatio, widthPara;
	double   internalCO2Ratio;
	double f,  
		scatt;//spectral correction
	int
		Kc25 = 0,
		Ko25 = 0,
		Kp25 = 0;
	double
		gbs = 0,
		gi = 0,
		LfWidth = 0,
		LfAngFact = 0;
	

};
#endif