/***************************************************************************
 * MathFuncs.h
 *
 * Author: Christian Wolf
 * Begin: 27.3.2002
 ***************************************************************************/

#ifndef _WOLF_MATHFUNCS_H_
#define _WOLF_MATHFUNCS_H_

// C
#include <math.h>

// From this module
#include "Vector.h"
#include "Matrix.h"
#include "DataSerie.h"
// Template Routines are included at 
// the end of this file!!!!!

/***************************************************************************
 * Statistics
 ***************************************************************************/

// Tables with significance levels of statistical tests
float getLimitKolmogorovSmirnov (float alpha, int N);

/***************************************************************************
 * Vector distances
 ***************************************************************************/

inline double Bhattacharyya (double x1, double x2, double v1, double v2);
template <class T> double MahalanobisDiff (Vector<T> &diff_vec, Matrix<T> &invcov);
template <class T> double Mahalanobis (Vector<T> &x, Vector<T> &y, Matrix<T> &invcov);

/***************************************************************************
 * Misc
 ***************************************************************************/

inline double NeuroFuzzyMinimum (double x1, double x2, double alpha)
{
	double vadiff=alpha*(x1-x2);
	double vtanh=tanh(vadiff);
	return 0.5*(x1+x2-log(cosh(vadiff)*(vtanh*vtanh+1)));
}

/***************************************************************************
 * Root finding
 ***************************************************************************/
 
template <class T> void rootWithPlot (T functional, float a, float b, 
	unsigned int noValues,  ostream &st);
template <class T> float rootWithBisection (T functional, float a, float b, float precision);

/***************************************************************************
 * Function and distribution fitting
 ***************************************************************************/
 
#ifdef HAVE_NUM_RECIPES
void fitGeneralizedGaussian (DataSerie &s, float &outMean, float &outVar, float &outShape);
#endif
void fitLogNormal           (DataSerie &s, float inPosition, float &outMean, float &outVar);

/**************************************************************************	
 *	Include a couple of template routines defined outside
 **************************************************************************/

#include "TemplatesRootFinding.h"
 
#endif
