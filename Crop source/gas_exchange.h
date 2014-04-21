#pragma once
double QuadSolnUpper (double a, double b, double c );
double QuadSolnLower (double a, double b, double c );

class CGas_exchange
{
public:
	CGas_exchange(string sType_in, double n_content);
	~CGas_exchange(void);

private:
  double PFD, R_abs, Tair, CO2, RH, wind, age, SLA, width, Press, N; //leafp is leaf water pressure passed from 2dsoil. Yang 8/20/06
  double leafp, leafpEffect;  //leafp is leaf water pressure, and leafpEffect is the effect of leafp on stomatal conductance. Yang 8/20/06
  double lfNContent; //lfNContent is leaf nitrogen content in unit of g m-2(leaf) YY
  void GasEx_psil(double leafp, double ET_supply);

   void Photosynthesis(double pressure); 
   void EnergyBalance(double Jw);
   void getParms();
   double gbw();
   double gsw(double pressure);
   double Es(double T);
   double Slope(double T);
   string sType;
public:
   void SetVal_psil(double PFD, double Tair, double CO2, double RH, double wind, double Press, double width, double leafp, double ET_supply);
   double get_gs(){return gs;}
  double get_leafpEffect(){return this->leafpEffect;}
  double set_leafpEffect(double pressure);
  double get_VPD(){ return VPD;}
  struct tparms
  {
	  double Vpm25, Vcm25, Jm25, Rd25, EaVp, EaVc, Eaj, Sj, Hj, Ear, g0, g1, g2;
	  double a_ABA, beta_ABA, delta, lamda_l, lamda_r, K_max; 
  } Parms;
  double A_gross, A_net, ET, Tleaf, Ci, gs, gb, Rd, iter1, iter2, VPD; 
  
  //double g_s;
};
