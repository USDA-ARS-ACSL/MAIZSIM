// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#include <iostream>
//#include <tchar.h>
// atlcoll is for arrays
//#include <atlcoll.h>
// atlstr is for atl string functions
//#include <atlstr.h> 

#define __min(a, b) (((a)<(b))?(a):(b))
#define __max(a, b) (((a)>(b))?(a):(b))

#include <cstring>
#ifdef _WIN32
#include <minmax.h>
#endif
//KY workaround for Windows-only functions
#ifndef _WIN32
#define strcpy_s(t, s) strcpy(t, s)
#define strcat_s(t, s) strcat(t, s)
#define strtok_s(s, d, c) strtok(s, d)
#define _itoa_s(i, a, n) snprintf(a, n, "%d", i)
#endif

// TODO: reference additional headers your program requires here
