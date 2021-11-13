/***************************************************************************
 * TemplatesRootFinding.h
 *
 * Author: Christian Wolf
 * Begin: 27.6.2005
 ***************************************************************************/

// C
#include <math.h>

// C++
#include <iostream>

// From the main module
#include <CIL.h>

// From this module
#include "MathFuncs.h"

#define MAX_BISECTIONS		50

/***************************************************************************
 * A method which dumps the function values in a given interval.
 * This permits to plot the function values for a visual check,
 * which helps to choose the root finding method.
 ***************************************************************************/
 
template <class T>
void rootWithPlot (T functional, float a, float b, unsigned int noValues, ostream &st)
{
	float stepWidth;
	
	if (a>b) 
	{
		float tmp=a;
		a=b; 
		b=tmp;
	}
	stepWidth=(b-a)/(float)noValues;
	
	for (float i=a; i<=b; i+=stepWidth)
		st << i << "\t" << functional(i) << endl;
}

/***************************************************************************
 * Root finding
 * T is the argument which specifies the function.
 * Actually, it is a class whose operator () must be overloaded and 
 * evaluate to the function.
 ***************************************************************************/
 
template <class T>
float rootWithBisection (T functional, float a, float b, float precision)
{	
	float x_left, x_right, x_midpoint,
	      f_left, f_right, f_midpoint;

	x_left=a;
	x_right=b;
	f_left=functional(x_left);
	f_right=functional(x_right);
	if (f_left*f_right >= 0.) 
		ERR_THROW ("rootWithBisection(): interval endpoints have same sign!");
		
	while (1)
	{	
		// cerr << "Bisection: " << (x_right+x_left)/2. << " " << (x_right+x_left) << endl;				
		
		x_midpoint=(x_left+x_right)/2.;		
		f_midpoint=functional(x_midpoint);
			
		if (f_midpoint*f_left >= 0.) 
			x_left=x_midpoint;
		else
			x_right=x_midpoint;
		
		if (fabs(x_right-x_left) < precision || f_midpoint == 0.) 
			return (x_right+x_left)/2.;
	}
}
