//Timer.cpp
//Member functions for Class Timer
#include "stdafx.h"
#include "timer.h"
#include <cmath>
#include <iostream>

Timer::Timer()
{
	date=NULL;
}

Timer::Timer(int dy, int my, int yr, double ss)
{
	day_of_year=julday(my,dy,yr)-julday(1,1,yr)+1;// day_number
	hour = 0.0;	// hr
	
	date = new CPMDate(day_of_year, yr);
	step_size = ss;	// hr
	steps_per_day = 24.0/step_size;
}

Timer::~Timer()
{
	//dt Dec 11, 2007 I'm not sure if this code is necessary. Time is not used to manage dates in this
	// way and we will never have a need to declared a date object. We have to declare it as null in
	// any constructor for Timer except if we want a date.
   if (date != NULL)  delete date;
}

Timer::Timer(double ss)
{
	step_size = ss;
	hour = 0.0;
	steps_per_day = 24.0/step_size;
	date=new CPMDate();
}

CPMDate* Timer::get_date()
{
	return date;
}

bool Timer::step()
{
//	hour = 0.0;
	hour += step_size;
//	cout << day_of_year << " " << hour << "\n";
	if (hour >= 24.0)
	{
		hour = 0.0;
		day_of_year += 1;
		date->increment();
		return true;
	}
	return false;
}

char* Timer::date_from_day_of_year(int d)
{
	CPMDate *dt = new CPMDate(d, (date->get_year()));
//	return dt->get_date();
	return dt->get_digit_string(); // set date format = mm/dd/yy, SK
}
int Timer::julday(const int mm, const int id, const int iyyy)
{
	const int IGREG=15+31*(10+12*1582);

	int ja,jul,jy=iyyy,jm;

	if (jy == 0) throw("julday: there is no year zero.");
	if (jy < 0) ++jy;
	if (mm > 2) {
		jm=mm+1;
	} else {
		--jy;
		jm=mm+13;
	}
	jul = int(floor(365.25*jy)+floor(30.6001*jm)+id+1720995);
	if (id+31*(mm+12*iyyy) >= IGREG) {
		ja=int(0.01*jy);
		jul += 2-ja+int(0.25*ja);
	}
	//jul=jul-2415019;    // 1/1/1900 JD zero
	jul=jul-2415079;   // on 3/1/1900 jday = 1
	return jul;
}
void Timer::caldat(int julian, int &mm, int &id, int &iyyy)
{
	const int IGREG=2299161;
	int ja,jalpha,jb,jc,jd,je;

	julian=julian+2415019 ;      //provide a reference of 1/1/1900 was 79 for 3/1/1900
	if (julian >= IGREG) {
		//jalpha=int((DP(julian-1867216)-0.25)/36524.25);
		jalpha=int(((double)(julian-1867216)-0.25)/36524.25);
		ja=julian+1+jalpha-int(0.25*jalpha);
	} else if (julian < 0) {
		ja=julian+36525*(1-julian/36525);
	} else
		ja=julian;
	jb=ja+1524;
	jc=int(6680.0+((double)(jb-2439870)-122.1)/365.25);
	jd=int(365*jc+(0.25*jc));
	je=int((jb-jd)/30.6001);
	id=jb-jd-int(30.6001*je);
	mm=je-1;
	if (mm > 12) mm -= 12;
	iyyy=jc-4715;
	if (mm > 2) --iyyy;
	if (iyyy <= 0) --iyyy;
	if (julian < 0) iyyy -= 100*(1-julian/36525);
	//sprintf(Date,"%i2/%i2/%i4",mm,id,iyyy);
}