/***************************************************************************
 * FunctionFitting.cpp
 *
 * Author: Christian Wolf
 * Begin: 2.7.2005
 ***************************************************************************/
 
// C
#include <math.h>

// C++
#include <iostream>

// From the main module
#include <CIL.h>

// From the NUMERICAL RECIPES library
#ifdef HAVE_NUM_RECIPES
#include <NumericalRecipes.h>
#endif

// From this module
#include "MathFuncs.h"
#include "DataSerie.h"

/***************************************************************************
 * Fit a generalized gaussian distribution to the data in the dataserie
 ***************************************************************************/
 
#ifdef HAVE_NUM_RECIPES 
// Functional necessary to numerically solve the 
// equation which determines the shape parameter
class ShapeFunctional
{
	public:
		ShapeFunctional (float a)				{ alpha=a; }
		float operator() (float x) 
		{
			return exp(numrec::gammln(5./x)+numrec::gammln(1./x)-2.*numrec::gammln(5./x))-alpha;	
		}
	private:
		float alpha;
}; 
 
void fitGeneralizedGaussian (DataSerie &s, float &outmean, float &outvar, float &outshape)
{
	float m,m4;
	
	// Get the necessary moments
	outmean = m = s.getMean();
	outvar  =     s.getVariance();	
	
	// Calculate the shape parameters
	m4=0.;
	for (DataSerie::iterator iter = s.begin(); iter!=s.end(); ++iter)
		m4 += pow(m-*iter,4);
	m4 /= s.size();			
	
	cerr << "fitGeneralizedGaussian(): m=" << m << " v=" << outvar << " m4=" << m4 << endl;
	
	try 
	{
		outshape = rootWithBisection (ShapeFunctional(m4/(outvar*outvar)), 0.1, 10, 0.001);
	}
	catch (EError &e)
	{
		ERR_THROW ("fitGeneralizedGaussian(): fit is numerically unstable!!");		
	}
}
#endif // HAVE_NUM_RECIPES

/***************************************************************************
 * Fit a log normal distribution to the data in the dataserie
 * s .......... the data
 * inPosition.. a location (shift) parameter
 * outMean..... the resulting mean parameter
 * outVar...... the resulting variance parameter
 ***************************************************************************/
 
void fitLogNormal (DataSerie &s, float inPosition, float &outMean, float &outVar)
{
	float m, v;
	
	// Calculate the mean
	m=0.;
	for (DataSerie::iterator iter = s.begin(); iter!=s.end(); ++iter)		
		m += log(*iter + inPosition);
	m /= s.size();			
	
	// Calculate the variance
	v=0.;
	for (DataSerie::iterator iter = s.begin(); iter!=s.end(); ++iter)
		v += pow(log(*iter + inPosition) - m, 2.);
	v /= s.size();			
		
	outMean=m;
	outVar=v;
}
 
 

