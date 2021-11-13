/****************************************************************
 * TestKolmogorovSmirnov.cpp
 * 
 * A table with the significance levels for the KS test, together
 * with the evaluation method.
 *
 * Christian Wolf
 * chriswolf@gmx.at
 * Begin: 7.12.2004
 ****************************************************************/

// C
#include <math.h>
 
 // C++
#include <iostream>

// From the main module
#include <CIL.h>

// From this module
#include "MathFuncs.h"


using namespace std;
 
static float TableKolmogorovSmirnov[][5] = 
{
/* alpha: 0.25  0.20  0.10 0.05  0.01 */
/* 1 */  {0.9  ,0.925,0.95,0.975,0.995},	// index 0
/* 2 */  {0.684,0.726,0.776,0.842,0.929},
/* 3 */  {0.565,0.597,0.642,0.708,0.828},
/* 4 */  {0.494,0.525,0.564,0.624,0.733},
/* 5 */  {0.446,0.474,0.51,0.565,0.669},
/* 6 */  {0.41,0.436,0.47,0.521,0.618},
/* 7 */  {0.381,0.405,0.438,0.486,0.577},
/* 8 */  {0.358,0.381,0.411,0.457,0.543},
/* 9 */  {0.339,0.36,0.388,0.432,0.514},
/* 10 */ {0.322,0.342,0.368,0.41,0.49},
/* 11 */ {0.307,0.326,0.352,0.391,0.468},
/* 12 */ {0.295,0.313,0.338,0.375,0.45},
/* 13 */ {0.284,0.302,0.325,0.361,0.433},
/* 14 */ {0.274,0.292,0.314,0.349,0.418},
/* 15 */ {0.266,0.283,0.304,0.338,0.404},
/* 16 */ {0.258,0.274,0.295,0.328,0.392},
/* 17 */ {0.25,0.266,0.286,0.318,0.381},
/* 18 */ {0.244,0.259,0.278,0.309,0.371},
/* 19 */ {0.237,0.252,0.272,0.301,0.363},
/* 20 */ {0.231,0.246,0.264,0.294,0.356},	// index 19: for N in [20,22]
/* 25 */ {0.21,0.22,0.24,0.27,0.32},		// index 20: for N in [23,27]
/* 30 */ {0.19,0.2,0.22,0.24,0.29},			// index 21: for N in [28,32]
/* 35 */ {0.18,0.19,0.21,0.23,0.27},		// index 22: for N in [33,38]

/* OVER 35: multiply with 1/sqrt(N) */
	     {1.07,1.14,1.22,1.36,1.63}			// index 23: for N >= 39
};

// *********************************************************************
// return the limit of the test statistics for the given significance 
// level. If calculated ratio is greater than value shown, then reject 
// the null hypothesis at the chosen level of confidence.
// 
// alpha ... error of the first kind (significance level)
// N ....... sample size
// *********************************************************************

float getLimitKolmogorovSmirnov (float alpha, int N)
{
	int alpha_index=0;
	
	switch ((int) rint(100*alpha)) 
	{
		case 20: alpha_index=0; break;
		case 15: alpha_index=1; break;
		case 10: alpha_index=2; break;
		case 5:  alpha_index=3; break;
		case 1:  alpha_index=4; break;
		default:
			ERR_THROW ("Internal error in getLimitKolmogorovSmirnov(): don't have table for the\n"
				 << "given significance level! (alpha=" << alpha << ").\n");
	}

	if (N<=20) 
		return TableKolmogorovSmirnov[N][alpha_index-1];	
		
	if (N<=22)
		return TableKolmogorovSmirnov[19][alpha_index-1];		
		
	if (N<=27)
		return TableKolmogorovSmirnov[20][alpha_index-1];		
		
	if (N<=32)
		return TableKolmogorovSmirnov[21][alpha_index-1];		
		
	if (N<=38)
		return TableKolmogorovSmirnov[22][alpha_index-1];		
		
	return TableKolmogorovSmirnov[23][alpha_index-1]/::sqrt(N);		
}
