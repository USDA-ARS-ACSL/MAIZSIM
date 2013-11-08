//cpmdate.h
// Class CPMDate
/*****************************************************
//		Author/Programmer:	R.L. Olson 03/14/97
/*****************************************************
			A class that implements the date.
/*****************************************************
	char
		date[40]
			date as string
		year_string[20]
			year as string
		month_string[30]
			month as string
		day_string[20]
			day as string
	int
		day_number
			day number
		year
			year
		day
			day of month
		month
			month
		void y_to_str()
			converts year from int to string
		void m_to_str()
			converts month from int to string
		void d_to_str()
			converts day from int to string
		void set_date()
			sets date from already stored month, day, year
		void dn_to_d()
			converts day number to day of month

	public:
		CPMDate()
			constructor
		CPMDate(int dn, int y)
			constructor:  dn = day_number, y = year
		CPMDate(int m, int d, int y)
			constructor:  m = month, d = day, y = year
		int get_day_number()
			accessing for day_number
		int get_year()
			accessing for year
		int get_day()
			accessing for day
		int get_month()
			accessing for month
		char* get_year_string()
			accessing for month
		char* get_day_string()
			accessing for day_string
		char* get_month_string()
			accessing for month_string
		char* get_date()
			accessing for date
		void increment()
			increments 1 day
		int is_leap_year()
			returns 0 if false, 1 if true
/*****************************************************
/*****************************************************/
#ifndef __CPMDATE_H__
#define __CPMDATE_H__

#include <string>

class CPMDate
{
	char year_string[20];
	char month_string[30];
	char day_string[20];
	char date[40];
	char digit_string[9];
	int day_number;
	int year;
	int day;
	int month;
	void d_to_str();
	void dn_to_d();
	void m_to_str();
	void set_date();
	void y_to_str();

public:
	CPMDate();
	CPMDate(int dn, int y);
	CPMDate(int m, int d, int y);
	char* get_date() {return date;}
	int get_day() {return day;}
	int get_day_number() {return day_number;}
	char* get_day_string() {return day_string;}
	char* get_digit_string() {return digit_string;}
	int get_month() {return month;}
	char* get_month_string() {return month_string;}
	int get_year() {return year;}
	char* get_year_string() {return year_string;}
	void increment();
	int is_leap_year();
	void make_digit_string();

};

#endif