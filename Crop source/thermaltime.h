#pragma once
#define MINUTESPERDAY (24.0*60.0);

class CThermalTime
{
public:
	CThermalTime(void);
	virtual ~CThermalTime(void);
	double get_Tcur() {return Tcur;}
	double get_Tbase() {return Tbase;}
	double get_Topt() {return Topt;}
	double get_Tmax() {return Tmax;}
	double get_sum() {return sum;}
	double get_dTmpr() {return dTmpr;}
	double get_timeStep() {return timeStep;}
	double get_actualAge() {return actualAge;}
	void set_Tcur(double x) {Tcur = x;}
	void set_Tbase(double x) {Tbase = x;}
	void set_Topt(double x) {Topt=x;}
	void set_Tmax(double x) {Tmax=x;}
	void set_timeStep(double x) {timeStep=x;}
	void set_temperatures(double Tb, double To, double Tm) {Tbase = Tb; Topt = To; Tmax = Tm;}
	void add(double x);
    void update(double Tmpr, double step);
	void initialize(double step);

private:
	double Tcur;
	double Tbase;
	double Topt;
	double Tmax;
	double sum;
	double dTmpr;
	double timeStep;
	double actualAge;

};
