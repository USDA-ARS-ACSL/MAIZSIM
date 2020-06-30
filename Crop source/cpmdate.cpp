//cpmdate.cpp
#include "stdafx.h"
#include "cpmdate.h"
#include <cmath>
#include <cstdio>
CPMDate::CPMDate()
{
}

CPMDate::CPMDate(int dn, int y)
{
	day_number = dn;
	year = y;

	dn_to_d();	
	y_to_str();
	m_to_str();
	d_to_str();	
	set_date();
	make_digit_string();
}

void CPMDate::dn_to_d()
{
	int f = is_leap_year();
	if (day_number <=31)
	{
		month = 1;
		day = day_number;
		return;  //DCA
	}
	if ((day_number > 31) && (day_number <= (59+f)))
	{
		month = 2;
		day = (day_number - (31+f));
		return;  //DCA
	}
	if (((day_number > (59+f)) && (day_number <= (90+f))))
	{
		month = 3;
		day = (day_number - (59 + f));
		return;  //DCA
	}
	if (((day_number > (90+f)) && (day_number <= (120+f))))
	{
		month = 4;
		day = (day_number - (90+f));
		return;  //DCA
	}
	if (((day_number > (120+f)) && (day_number <= (151+f))))
	{
		month = 5;
		day = (day_number - (120+f));
		return;  //DCA
	}
	if (((day_number > (151+f)) && (day_number <= (181+f))))
	{
		month = 6;
		day = (day_number - (151+f));
		return;  //DCA
	}
	if (((day_number > (181+f)) && (day_number <= (212+f))))
	{
		month = 7;
		day = (day_number - (181+f));
		return;  //DCA
	}
	if (((day_number > (212+f)) && (day_number <= (243+f))))
	{
		month = 8;
		day = (day_number - (212+f));
		return;  //DCA
	}
	if (((day_number > (243+f)) && (day_number <= (273+f))))
	{
		month = 9;
		day = (day_number - (243+f));
		return;  //DCA
	}
	if (((day_number > (273+f)) && (day_number <= (304+f))))
	{
		month = 10;
		day = (day_number - (273+f));
		return;  //DCA
	}
	if (((day_number > (304+f)) && (day_number <= (334+f))))
	{
		month = 11;
		day = (day_number - (304+f));
		return;  //DCA
	}
	if (((day_number > (334+f)) && (day_number <= (365+f))))
	{
		month = 12;
		day = (day_number - (334+f));
		return;
	}
}


CPMDate::CPMDate(int m, int d, int y)
{
	year = y;
	month = m;
	day = d;
	y_to_str();
	m_to_str();
	d_to_str();
	set_date();
}

void CPMDate::make_digit_string()
{
	digit_string[0] = '\0';
	int yrr;
	char yy[10], yyy[10], dd[10], ddd[10], mm[10], mmm[10];
	yy[0] = '\0';
	yyy[0] = '\0';
	mm[0] = '\0';
	mmm[0] = '\0';
	dd[0] = '\0';
	ddd[0] = '\0';

	if (year > 1999)
	{ 
		yrr = year - 2000;		
	}
	else yrr = year - 1900;
	_itoa_s(yrr, yy, 10);
	if (yrr < 10)
	{
		strcat_s(yyy, "0");
		strcat_s(yyy, yy);
	}
	else strcpy_s(yyy, yy); 
	_itoa_s(month, mm, 10);
	if (month < 10)
	{
		strcat_s(mmm, "0");
		strcat_s(mmm, mm);
	}
	else strcpy_s(mmm, mm); 
	_itoa_s(day, dd, 10);
	if (day < 10)
	{
		strcat_s(ddd, "0");
		strcat_s(ddd, dd);
	}
	else strcpy_s(ddd, dd); 

	strcat_s(digit_string, mmm);
	strcat_s(digit_string, "/");
	strcat_s(digit_string, ddd);
	strcat_s(digit_string, "/");
	strcat_s(digit_string, yyy);
}



void CPMDate::y_to_str()
{
	char ybuffer[20];
	_itoa_s(year, ybuffer, 10);
	strcpy_s(year_string, ybuffer);
}

void CPMDate::m_to_str()
{
	char buffer[30];
	switch(month)
	{
		case 1:
			strcpy_s(buffer, "Jan");
			break;
		case 2:
			strcpy_s(buffer, "Feb");
			break;
		case 3:
			strcpy_s(buffer, "Mar");
			break;
		case 4:
			strcpy_s(buffer, "Apr");
			break;
		case 5:
			strcpy_s(buffer, "May");
			break;
		case 6:
			strcpy_s(buffer, "Jun");
			break;
		case 7:
			strcpy_s(buffer, "Jul");
			break;
		case 8:
			strcpy_s(buffer, "Aug");
			break;
		case 9:
			strcpy_s(buffer, "Sep");
			break;
		case 10:
			strcpy_s(buffer, "Oct");
			break;
		case 11:
			strcpy_s(buffer, "Nov");
			break;
		case 12:
			strcpy_s(buffer, "Dec");
			break;
	}
	strcpy_s(month_string, buffer);
}

void CPMDate::d_to_str()
{
	char dd[20], tbuffer[4];
	//put out in standard format, SK
	_itoa_s(day, tbuffer, 4);
	size_t origsize = strlen(dd);

	if (day < 10)
	{
		strcpy_s(dd,"0");
		strcat_s(dd,tbuffer);
	}
	else
	{
		strcpy_s(dd, tbuffer);
	}
	strcpy_s(day_string, dd);
}

void CPMDate::set_date()
{
	char dbuffer[40];
	strcpy_s(dbuffer, month_string);
	strcat_s(dbuffer, " ");
	strcat_s(dbuffer, day_string);
	strcat_s(dbuffer, " ");
	strcat_s(dbuffer, year_string);
	strcat_s(dbuffer, " ");
	strcpy_s(date, dbuffer);
}

int CPMDate::is_leap_year()
{
	if (fmod(year, 4.0) == 0.0)
	{
		return 1;
	}
	else return 0;
}

void CPMDate::increment()
{
	int j = 0;
	if (is_leap_year()) j = 1;
	if (day_number == (365+j))
	{
		day_number = 1;
		year++;
	}
	else day_number++;
	dn_to_d();	
	y_to_str();
	m_to_str();
	d_to_str();	
	set_date();
	make_digit_string();
}
