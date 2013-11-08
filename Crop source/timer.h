//Timer.h
//Class Timer vs. 1.0
/*****************************************************
//		Author:	R.L. Olson 03/12/97
//		Programmer: R.L. Olson
//      Modified by: S.-H. Kim 10/15/04
/*****************************************************
//		Description
/*****************************************************
	A class that implements a simulation timer
/*****************************************************
//		Data Dictionary
/*****************************************************
	int
		day_of_year			// day number
	double
		hour				// hour of day
		step_size			// step-size in hours
	CPMDate *date			// date object
/*****************************************************
		Interface
/*****************************************************
	Timer()					// constructor
	Timer(double ss)			// constructor (ss = step size)
	Timer(int day, int year, double ss)
							//day = day number, ss = step size

	char* get_date()			//returns date as string
	char* date_from_day_of_year(int d)
			// returns date as string from day number
	int get_day_of_year()		// accessing
	double get_hour()			// accessing
	double get_step_size()		// accessing
	void set_hour(double d)		// accessing
	void step()					//move forward one timestep

/*****************************************************
/*****************************************************/
#ifndef __TIMER_H__
#define __TIMER_H__
#include "cpmdate.h"

class Timer
{
private:
	CPMDate *date;
	double hour, step_size, steps_per_day;
	int day_of_year;

public:
	char* date_from_day_of_year(int d);
	CPMDate* get_date();
	int get_day_of_year() {return date->get_day_number();}
	double get_hour() {return hour;}
	double get_step_size() {return step_size;}
	double get_steps_per_day() {return steps_per_day;}
//BASept99 deleted next method - not used
//	void set_hour(double d) {hour = d;}
	bool step();
	 ~Timer();
	 Timer();
	 Timer(double ss);
	 Timer(int day, int month, int year, double ss);

 	int julday(const int mm, const int id, const int iyyy);
	void caldat(const int julian, int &mm, int &id, int &iyyy);

};
#endif