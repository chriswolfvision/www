#ifndef _WOLF_BICUBICPOLYNOM_H_
#define _WOLF_BICUBICPOLYNOM_H_

#include <iostream>
#include <math.h>

// *************************************************************
// The polynom function for the bicubic interpolation.
// The values are pre-computed.
// *************************************************************
		
class BicubicPolynom {

	public:

		// The instance	function of	the	singleton
		static BicubicPolynom *instance(int ifactor) {
			if (pInstance == NULL)
				pInstance =	new	BicubicPolynom (ifactor);
			return pInstance;
		}
		
		// The access operator
    	double operator() (double x) {
    		return values[(int) (int_factor*(x+2.0))];
    	}
			
	private:
	
		// The constructor
		BicubicPolynom (int ifactor);
		
		// The destructor
		~BicubicPolynom () {
			delete values;
			pInstance = NULL;
		}
		
		// The polynom function
    	double calcPolynom (double x) {
    		return (p(x+2.0,3.0) - 4.0*p(x+1.0,3.0) +
    		    6.0*p(x,3.0) - 4.0*p(x-1.0,3.0))/6.0;
    	}
    	
    	// The power function for the polynom
    	double p(double a, double b) {
    		return a>0 ? pow (a,b) : 0;
    	}
    	
    private: // STATIC MEMBERS
    	
	   	static BicubicPolynom *pInstance;
		static double *values;
		static int int_factor;
};


#endif
